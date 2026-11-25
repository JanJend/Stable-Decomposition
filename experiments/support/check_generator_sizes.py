"""Run with python3 -m experiments.support.check_generator_sizes (build pruning_torus_runtime first)."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from . import HERE, ROOT
from .runtime import digest


def main():
    sys.path.insert(0, str(HERE))
    import benchmark_generator_sizes as experiment

    worker = ROOT / "build/pruning_torus_runtime"
    fixture = ROOT / "test_presentations/test8.scc"
    with tempfile.TemporaryDirectory(prefix="generator-sizes-check-") as name:
        temporary = Path(name)
        sample = dict(input_generators=6, input_relations=8, presentation_size=14)
        for fast_i in (False, True):
            by_thread = []
            for threads in (1, 2):
                command = [str(worker), str(fixture), "10", str(threads), "512"]
                if not fast_i:
                    command.append("--standard-i")
                plain = subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
                log = temporary / "trace.tsv"
                log.write_text(subprocess.run([*command, "--generator-trace"], check=True,
                                              capture_output=True, text=True, timeout=30).stdout)
                plain_log = temporary / "plain.tsv"
                plain_log.write_text(plain.stdout)
                metrics, _ = experiment.parse_trace(log)
                plain_metrics, _ = experiment.parse_trace(plain_log)
                events = experiment.parse_generators(log)
                experiment.validate_complete(metrics, events, sample, fast_i)
                for key in ("i_iterations", "k_iterations", "i_maps", "k_maps", "output_generators", "output_relations"):
                    assert metrics[key] == plain_metrics[key], key
                by_thread.append({obj: [{key: value for key, value in event.items() if key != "elapsed_seconds"}
                                        for event in events if event["object"] == obj] for obj in ("I", "K", "J", "L")})
                # Raw columns include redundant generators; the final submodules are smaller.
                if not fast_i:
                    summary = experiment.size_summary(events, sample)
                    assert (summary["i_peak_generators"], summary["i_peak_reduced_generators"],
                            summary["i_final_generators"], summary["k_final_generators"]) == (10, 8, 8, 4)
                    assert summary["i_peak_generators_over_x_generators"] == 10 / 6
                    assert summary["i_peak_generators_over_presentation_size"] == 10 / 14
            assert by_thread[0] == by_thread[1], "Serial and parallel generator histories differ"

        output = temporary / "results"
        command = [sys.executable, "-B", str(HERE / "benchmark_generator_sizes.py"), "--input", str(fixture),
                   "--worker", str(worker), "--output-dir", str(output), "--percentages", "0", "10", "--threads", "2"]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 2 and all(row["status"] == "ok" and row["fast_i"] == "False" for row in rows)
        assert rows[0]["i_iterations"] == rows[0]["k_iterations"] == "0"  # Q=0 shortcut.
        assert rows[0]["i_final_generators"] == "6" and rows[0]["k_final_generators"] == "0"
        assert rows[1]["i_peak_generators"] == "10" and rows[1]["k_peak_reduced_generators"] == "4"
        settings = json.loads(next(output.glob("*/settings.json")).read_text())
        assert settings["worker_sha256"] == digest(worker) and settings["fast_i"] is False
        with (output / "generators.csv").open() as stream:
            assert all(float(row["over_presentation_size"]) == int(row["generators"]) / 14
                       for row in csv.DictReader(stream))

        # Empty X must retain zero counts but leave undefined ratios blank.
        empty = temporary / "empty.scc"
        empty.write_text("scc2020\n2\n0 0 0\n")
        subprocess.run([*command, "--input", str(empty), "--percentages", "0"], check=True,
                       capture_output=True, text=True, timeout=30)
        with (output / "results.csv").open() as stream:
            row = next(csv.DictReader(stream))
        assert row["i_peak_generators"] == row["k_peak_generators"] == "0"
        assert all(row[key] == "" for key in experiment.RATIO_COLUMNS)

        # Exercise manifest selection and a sweep on a tiny, deliberately relabelled fixture.
        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("family", "dataset", "target_points", "n_points",
                                                        "homology", "seed", "scc", "scc_sha256"))
            writer.writeheader()
            for family, homology, n in (("torus_gauss", 1, 640), ("gauss_clust", 0, 640), ("immuno", 1, 641)):
                writer.writerow(dict(family=family, dataset=family, target_points=640, n_points=n,
                                     homology=homology, seed=1, scc=str(fixture), scc_sha256=digest(fixture)))
        manifest_command = [sys.executable, "-B", str(HERE / "benchmark_generator_sizes.py"), "--manifest", str(manifest),
                            "--worker", str(worker), "--output-dir", str(output), "--sizes", "640", "--percentages", "10"]
        subprocess.run(manifest_command, check=True, capture_output=True, text=True, timeout=30)
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert [row["family"] for row in rows] == ["torus_gauss", "gauss_clust", "immuno"]
        assert rows[-1]["n_points"] == "641" and all(row["status"] == "ok" for row in rows)

        # A timed-out worker's complete events survive, while its truncated last event is ignored.
        partial_worker = temporary / "partial_worker.py"
        partial_worker.write_text(f"#!{sys.executable}\n"
                                  "import sys, time\n"
                                  "if '--help' in sys.argv:\n"
                                  "    print('Generator size trace v1; Column storage: compact')\n"
                                  "else:\n"
                                  "    print('GENERATORS\\tI\\tinitial\\t0\\t0\\t6\\t0.01', flush=True)\n"
                                  "    print('GENERATORS\\tI\\traw\\t1\\t1\\t20\\t0.02', flush=True)\n"
                                  "    print('GENERATORS\\tI\\treduced', end='', flush=True)\n"
                                  "    time.sleep(10)\n")
        partial_worker.chmod(0o755)
        subprocess.run([*command, "--worker", str(partial_worker), "--percentages", "10", "--timeout", "0.5"],
                       capture_output=True, text=True, timeout=30, check=False)
        with (output / "results.csv").open() as stream:
            row = next(csv.DictReader(stream))
        assert row["status"] == "timeout" and row["peaks_are_lower_bounds"] == "True"
        assert row["i_peak_generators"] == "20" and row["i_final_generators"] == row["k_peak_generators"] == ""
        assert len(list(output.glob("*/settings.json"))) == 4 and list(output.glob("*/previous_results.csv"))
    print("PASS: known counts, serial/parallel histories, standard/fast I, unchanged output, Q=0, empty X, ratios, manifest, archives and partial timeouts")


if __name__ == "__main__":
    main()

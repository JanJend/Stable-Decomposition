"""Small worker/CSV/timeout check: python3 -m experiments.support.check_torus_runtime."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from . import DATA, HERE, ROOT
from .runtime import digest


def main():
    sys.path.insert(0, str(HERE))
    import benchmark_torus_runtime as experiment

    assert experiment.select_input(DATA / "manifest.csv")["dataset"] == "torus_gauss_1280_1"
    worker = ROOT / "build/pruning_torus_runtime"
    fixture = ROOT / "test_presentations/test8.scc"
    with tempfile.TemporaryDirectory(prefix="torus-runtime-check-") as name:
        temporary = Path(name)
        results = []
        for threads in (1, 2):
            log = temporary / f"{threads}.tsv"
            with log.open("w") as stream:
                subprocess.run([str(worker), str(fixture), "10", str(threads), "512"],
                               stdout=stream, stderr=subprocess.DEVNULL, check=True, timeout=30)
            metrics, steps = experiment.parse_trace(log)
            assert steps and all(step["completed"] and step["seconds"] >= 0 for step in steps)
            total = next(step for step in steps if step["operation"] == "total")
            assert total["seconds"] == metrics["pruning_seconds"]
            for phase in ("i", "k"):
                passes = {step["pass"] for step in steps if step["operation"] == phase + "_equal"}
                assert passes == set(range(1, metrics[phase + "_iterations"] + 1))
            assert {"end0", "hom_quotient", "composition", "map_lifts", "i_phase", "k_phase"} <= {
                step["operation"] for step in steps}
            results.append(metrics)
        for key in ("end_dimension", "quotient_dimension", "i_maps", "k_maps", "i_iterations",
                    "k_iterations", "output_generators", "output_relations"):
            assert results[0][key] == results[1][key], key

        # Overlapping worker calls must retain distinct start/end IDs and valid pass numbers.
        for flags in (("--standard-i",), ("--parallel-i",), ("--parallel-i", "--minimize-i-target"),
                      ("--parallel-k",), ("--standard-i", "--parallel-k"), ("--parallel-i", "--parallel-k")):
            log = temporary / ("-".join(flags).replace("--", "") + ".tsv")
            with log.open("w") as stream:
                subprocess.run([str(worker), str(fixture), "10", "8", "512", *flags],
                               stdout=stream, stderr=subprocess.DEVNULL, check=True, timeout=30)
            metrics, steps = experiment.parse_trace(log)
            assert steps and all(step["completed"] for step in steps)
            assert metrics["parallel_i"] == ("--parallel-i" in flags)
            assert metrics["parallel_k"] == ("--parallel-k" in flags)
            assert metrics["minimize_i_target"] == ("--minimize-i-target" in flags)
            assert metrics["fast_i"] == not ("--standard-i" in flags or "--parallel-i" in flags)
            assert any(step["operation"] == "i_shift" for step in steps)
            assert any(step["operation"] == "i_target_reduce" for step in steps) == metrics["minimize_i_target"]
            for key in ("output_generators", "output_relations"):
                assert metrics[key] == results[0][key], key

        partial = temporary / "partial.tsv"
        partial.write_text('METRICS {"threads":2}\n'
                           "STEP\tstart\ttotal\t0\t1\t0\t0\n"
                           "STEP\tstart\tend0\t0\t1\t0.1\t0\n"
                           "STEP\tend\tend0\t0\t1\t1.1\t1\n"
                           "STEP\tstart\thom_quotient\t0\t1\t1.2\t0\nSTEP\tend\thom_")
        metrics, steps = experiment.parse_trace(partial)
        assert metrics == {"threads": 2} and [step["completed"] for step in steps] == [False, True, False]
        summary = experiment.aggregate_steps([dict(step, epsilon_percent=1, repeat=1) for step in steps])
        assert summary[0]["unfinished_calls"] == 1 and summary[0]["calls"] == 0
        assert summary[1]["seconds"] == 1

        # The integration manifest deliberately labels a tiny fixture as the target input.
        manifest = temporary / "manifest.csv"
        sample = dict(dataset="test-torus", family="torus_gauss", target_points=1280, n_points=1280,
                      homology=1, seed=1, generation_status="ok", scc=str(fixture), scc_sha256=digest(fixture))
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=sample)
            writer.writeheader()
            writer.writerow(sample)
        output = temporary / "results"
        command = [sys.executable, str(HERE / "benchmark_torus_runtime.py"), "--manifest", str(manifest),
                   "--worker", str(worker), "--output-dir", str(output), "--percentages", "10",
                   "--threads", "8", "--repeats", "2", "--parallel-i", "--parallel-k",
                   "--composition-certificate-visits", "65536"]
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 2 and all(row["status"] == "ok" for row in rows)
        with (output / "steps.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert {row["repeat"] for row in rows} == {"1", "2"} and all(row["completed"] == "True" for row in rows)
        settings = json.loads(next(output.glob("*/settings.json")).read_text())
        assert settings["input"]["dataset"] == "test-torus" and settings["worker_sha256"] == digest(worker)
        assert settings["generator_options"]["composition-certificate-visits"] == 65536
        assert settings["generator_options"]["composition-nilpotent-products"] == 262144
        assert settings["parallel_i"] and settings["parallel_k"] and not settings["minimize_i_target"]
        with (output / "results.csv").open() as stream:
            assert all(row["fast_i"] == "False" for row in csv.DictReader(stream))
        invalid = subprocess.run([str(worker), str(fixture), "10", "1", "512", "--composition-nilpotent-products=-1"],
                                 capture_output=True, text=True, timeout=30)
        assert invalid.returncode != 0 and "requires a nonnegative integer" in invalid.stderr
        timed_out = subprocess.run([*command, "--timeout", "0.000001"], capture_output=True, text=True, timeout=30)
        assert timed_out.returncode == 1
        with (output / "results.csv").open() as stream:
            assert all(row["status"] == "timeout" and not row["pruning_seconds"] for row in csv.DictReader(stream))
        assert len(list(output.glob("*/settings.json"))) == 2
        assert len(list(output.glob("*/*.output.scc"))) == 2
    print("PASS: exact torus selection, serial/parallel traces, pass numbers, timer reset, partial events, CSV, archives and timeouts")


if __name__ == "__main__":
    main()

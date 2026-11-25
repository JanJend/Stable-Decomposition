"""Run with python3 -m experiments.support.check_end_quotient_benchmark."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from . import DATA, HERE, ROOT


def main():
    sys.path.insert(0, str(HERE))
    import benchmark_end_quotient as experiment

    selected = experiment.select_inputs(DATA / "manifest.csv", [1280, 640], ["torus_gauss", "immuno"])
    assert [(row["dataset"], row["n_points"]) for row in selected] == [
        ("torus_gauss_640_1", "640"), ("T_E_ROI_40_locations_CD8", "640"),
        ("torus_gauss_1280_1", "1280"), ("T_A_ROI_36_locations_CD8", "1281")]
    worker = ROOT / "build/local-release/end_quotient_benchmark"
    with tempfile.TemporaryDirectory(prefix="end-quotient-comparison-") as name:
        temporary = Path(name)
        free = temporary / "free.scc"
        relation = temporary / "relation.scc"
        free.write_text("scc2020\n2\n0 2 0\n0 1 ;\n1 0 ;\n")
        relation.write_text("scc2020\n2\n1 2 0\n5 2 ; 0\n0 0 ;\n1 0 ;\n")
        for path, percentage, dimensions in [(free, 0, (2, 0, 2)), (free, 50, (4, 2, 2)),
                                              (relation, 50, (1, 0, 1))]:
            for method in experiment.METHODS:
                result = subprocess.run([str(worker), str(path), str(percentage), method, "0"],
                                        check=True, capture_output=True, text=True, timeout=30)
                metrics = json.loads(next(line[7:] for line in result.stdout.splitlines() if line.startswith("RESULT ")))
                assert tuple(metrics[key] for key in experiment.DIMENSIONS) == dimensions
                assert metrics["method"] == method and metrics["seconds"] > 0
                assert metrics["End_0_dimension"] == (-1 if method == "end-d-0" else (2 if path == free else 3))

        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("dataset", "family", "target_points", "n_points", "homology", "seed", "scc"))
            writer.writeheader()
            for n in (20, 40):
                writer.writerow(dict(dataset=f"torus_{n}", family="torus_gauss", target_points=n,
                                     n_points=n, homology=1, seed=1, scc=free.name))
        output = temporary / "results"
        command = [sys.executable, str(HERE / "benchmark_end_quotient.py"), "--manifest", str(manifest),
                   "--worker", str(worker), "--output-dir", str(output), "--sizes", "40", "20",
                   "--families", "torus_gauss", "--percentages", "0.5", "50", "--repeats", "2"]
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        with (output / "results.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert [row["target_points"] for row in rows] == ["20", "20", "40", "40"]
        assert all(row["dimensions_match"] == "yes" and row["explicit_runs"] == row["end_d_0_runs"] == "2" for row in rows)
        assert [row["explicit_quotient_dimension"] for row in rows] == ["0", "2", "0", "2"]
        assert all(float(row["speedup_explicit_over_new"]) > 0 for row in rows)
        run = next(output.glob("*/settings.json")).parent
        records = [json.loads(line) for line in (run / "observations.jsonl").read_text().splitlines()]
        assert len(records) == 16 and len({row["stdout_log"] for row in records}) == 16
        assert records[0]["method"] == "explicit" and records[2]["method"] == "end-d-0"
        assert all((run / row["stdout_log"]).is_file() for row in records)

        # A mismatch or failed process must never produce a speedup.
        sample = dict(selected[0], bounding_box_sides=[1, 2])
        pair = {method: [dict(row)] for method, row in
                (("explicit", records[0]), ("end-d-0", records[1]))}
        pair["end-d-0"][0]["quotient_dimension"] += 1
        mismatch = experiment.comparison(sample, 1, pair)
        assert mismatch["dimensions_match"] == "NO" and "speedup_explicit_over_new" not in mismatch
        assert experiment.summarize([records[1], pair["end-d-0"][0]])["status"] == "inconsistent"
        pair["end-d-0"][0]["status"] = "memory_limit"
        incomplete = experiment.comparison(sample, 1, pair)
        assert incomplete["dimensions_match"] == "incomplete" and "end_d_0_seconds" not in incomplete
        assert experiment.summarize([records[0], pair["end-d-0"][0]])["status"] == "memory_limit"

        timed_out = subprocess.run([*command, "--timeout", "0.000001"], capture_output=True, text=True, timeout=30)
        assert timed_out.returncode == 1
        with (output / "results.csv").open(newline="") as stream:
            assert all(row["dimensions_match"] == "incomplete" and not row["speedup_explicit_over_new"]
                       for row in csv.DictReader(stream))
        assert list(output.glob("*/previous_results.csv"))
    print("PASS: both quotient routes, known dimensions, increasing sizes, repetitions/order, distinct logs, CSV, archives, mismatches and timeouts")


if __name__ == "__main__":
    main()

"""Small, deterministic checks for the experiment worker and timeout handling."""
import csv
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

from . import ROOT
from .runtime import run_worker


def check(worker):
    with tempfile.TemporaryDirectory(prefix="pruning-experiment-check-") as directory:
        directory = Path(directory)
        interval = directory / "interval.scc"
        # A rectangle of side lengths 1 and 10 has diagonal lifespan 1, not 10.
        interval.write_text("scc2020\n2\n2 1 0\n1 0 ; 0\n0 10 ; 0\n0 0 ;\n")
        metrics, status, _, error = run_worker(worker, ["components", interval, 0.49, 0.5, 1], 30)
        assert status == "ok", error
        assert metrics["components"] == 1, metrics
        assert metrics["short_epsilon"] == [0, 0, 1], metrics
        assert metrics["short_2epsilon"] == [0, 1, 1], metrics
        metrics, status, _, error = run_worker(worker, ["components", interval, "0.05,1", "0.5,5"], 30)
        assert status == "ok", error
        assert metrics["short_epsilon"] == [0, 0] and metrics["short_2epsilon"] == [0, 1], metrics
        metrics, status, _, error = run_worker(worker, ["inspect", interval], 30)
        assert status == "ok" and metrics["thickness"] == 1, (metrics, error)
        output = directory / "pruned.scc"
        metrics, status, _, error = run_worker(worker, ["prune", interval, 0, output], 30)
        assert status == "ok" and output.exists(), error
        assert metrics["pruning_seconds"] >= 0, metrics
        assert metrics["i_iterations"] == metrics["k_iterations"] == 0, metrics
        free = directory / "free.scc"
        free.write_text("scc2020\n2\n0 1 0\n0 0 ;\n")
        metrics, status, _, error = run_worker(worker, ["components", free, 100], 30)
        assert status == "ok", error
        assert metrics["components"] == 1 and metrics["short_2epsilon"] == [0], metrics
        zero = directory / "zero.scc"
        zero.write_text("scc2020\n2\n0 0 0\n")
        metrics, status, _, error = run_worker(worker, ["inspect", zero], 30)
        assert status == "ok" and metrics["thickness"] == metrics["bbox_extent"] == 0, (metrics, error)
        _, status, _, _ = run_worker(worker, ["inspect", interval], 1e-9)
        assert status == "timeout", status
        # Distinct datasets with identical SCCs must each receive a row. Resume
        # must then skip those exact rows without skipping other percentages.
        twin = directory / "free-twin.scc"
        twin.write_bytes(free.read_bytes())
        manifest = directory / "manifest.csv"
        rows = [dict(dataset=path.stem, family="check", source="check", n_points=20,
                     ambient_dimension=2, homology=0, scc=path.name,
                     scc_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), generation_status="ok")
                for path in (free, twin)]
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        table = directory / "results.csv"
        command = [sys.executable, "-m", "experiments.support.aida_runtime", "--worker", str(worker),
                   "--manifest", str(manifest), "--table", str(table), "--percentages", "1"]
        for extra in ([], ["--resume"]):
            subprocess.run(command + extra, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
            with table.open(newline="") as stream:
                assert len(list(csv.DictReader(stream))) == 2
        timeout_table = directory / "timeouts.csv"
        result = subprocess.run(command + ["--table", str(timeout_table), "--pruning-timeout", "0.000000001"],
                                cwd=ROOT, stdout=subprocess.DEVNULL)
        assert result.returncode == 1
        with timeout_table.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 2 and all(row["inspection_status"] == "timeout" for row in rows)
        assert all(row["components_before"] == row["pruning_seconds"] == "" for row in rows)
    print("PASS: lifespan thresholds, thickness, zero/free modules, profiling, deadlines, CSV failures, and resume")



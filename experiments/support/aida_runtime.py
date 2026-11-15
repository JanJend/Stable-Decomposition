#!/usr/bin/env python3
"""Measure pruning and AIDA, retaining only a CSV table, never pruned modules."""
import argparse
import csv
import os
from pathlib import Path
import platform
import tempfile

from . import DATA, HERE, ROOT
from .runtime import digest, run_worker
FIELDS = [
    "dataset", "family", "source", "n_points", "ambient_dimension", "homology", "scc",
    "input_sha256", "worker_sha256", "machine", "epsilon_percent", "epsilon",
    "bbox_min_x", "bbox_min_y", "bbox_max_x", "bbox_max_y", "bbox_extent",
    "generators", "relations", "thickness", "thickness_seconds", "inspection_status",
    "pruning_seconds", "pruning_wall_seconds", "i_iterations", "k_iterations", "pruning_status",
    "components_before", "components_before_le_epsilon", "components_before_le_2epsilon",
    "components_before_gt_2epsilon", "aida_input_seconds", "aida_input_wall_seconds",
    "lifetime_classification_seconds", "aida_input_status", "components_after",
    "aida_output_seconds", "aida_output_wall_seconds", "aida_output_status",
    "pruning_timeout_seconds", "aida_timeout_seconds", "error",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/experiments/pruning_experiment")
    parser.add_argument("--table", type=Path, default=HERE / "results/aida/results.csv")
    parser.add_argument("--percentages", type=float, nargs="+", default=[1, 2, 4, 8, 16])
    parser.add_argument("--pruning-timeout", type=float, default=600)
    parser.add_argument("--aida-timeout", type=float, default=1200)
    parser.add_argument("--max-points", type=int, help="Optional small validation run")
    parser.add_argument("--families", nargs="+", help="Optional subset of manifest families")
    parser.add_argument("--resume", action="store_true", help="Skip recorded input/worker/epsilon triples, including timeouts")
    args = parser.parse_args()
    if any(not (0 < p <= 100) for p in args.percentages):
        parser.error("Percentages must lie in (0,100]")
    if not (0 < args.pruning_timeout < float("inf") and 0 < args.aida_timeout < float("inf")):
        parser.error("Timeouts must be finite and positive")
    worker = args.worker.resolve()
    if not worker.is_file() or not os.access(worker, os.X_OK):
        parser.error(f"Build the pruning_experiment target first: missing executable {worker}")
    manifest = args.manifest.resolve()
    with manifest.open(newline="") as stream:
        samples = [row for row in csv.DictReader(stream) if row["generation_status"] == "ok"]
    samples = [row for row in samples if
               (args.max_points is None or int(row["n_points"]) <= args.max_points) and
               (not args.families or row["family"] in args.families)]
    samples.sort(key=lambda row: (int(row["n_points"]), row["family"], int(row["homology"]), row["dataset"]))
    if not samples:
        parser.error("No generated presentations match; run generate_data.py first")
    table = args.table.resolve()
    table.parent.mkdir(parents=True, exist_ok=True)
    worker_hash = digest(worker)
    recorded = set()
    if table.exists() and table.stat().st_size:
        with table.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != FIELDS:
                parser.error(f"Unexpected table schema in {table}; choose a new --table")
            if args.resume:
                recorded = {(row["scc"], row["input_sha256"], row["worker_sha256"], float(row["epsilon_percent"])) for row in reader}
    needs_header = not table.exists() or table.stat().st_size == 0
    failed = 0
    with table.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        if needs_header:
            writer.writeheader()
        for sample in samples:
            scc = manifest.parent / sample["scc"]
            input_hash = digest(scc)
            if input_hash != sample["scc_sha256"]:
                parser.error(f"Input changed since generation: {scc}; regenerate its manifest entry")
            percentages = [p for p in args.percentages if (sample["scc"], input_hash, worker_hash, p) not in recorded]
            if not percentages:
                continue
            print(f"Inspecting {sample['dataset']} ({sample['n_points']} points, H{sample['homology']})", flush=True)
            info, inspection_status, _, inspection_error = run_worker(worker, ["inspect", scc], args.pruning_timeout)
            if "bbox_extent" not in info:
                epsilons = []
                before, before_status, before_wall, before_error = {}, "skipped", "", "Missing module bounds"
            else:
                epsilons = [float(info["bbox_extent"]) * p / 100 for p in percentages]
                before, before_status, before_wall, before_error = run_worker(
                    worker, ["components", scc, *epsilons], args.aida_timeout)
            for index, percentage in enumerate(percentages):
                row = {key: sample[key] for key in ["dataset", "family", "source", "n_points", "ambient_dimension", "homology", "scc"]}
                row.update(input_sha256=input_hash, worker_sha256=worker_hash, machine=platform.machine(),
                           epsilon_percent=percentage, inspection_status=inspection_status,
                           aida_input_status=before_status, aida_input_wall_seconds=before_wall,
                           pruning_timeout_seconds=args.pruning_timeout, aida_timeout_seconds=args.aida_timeout,
                           pruning_status="skipped", aida_output_status="skipped")
                row.update({key: value for key, value in info.items() if key in FIELDS})
                errors = [error for error in [inspection_error, before_error] if error]
                if epsilons:
                    epsilon = epsilons[index]
                    row["epsilon"] = epsilon
                    if before_status == "ok":
                        row.update(components_before=before["components"],
                                   components_before_le_epsilon=before["short_epsilon"][index],
                                   components_before_le_2epsilon=before["short_2epsilon"][index],
                                   components_before_gt_2epsilon=before["components"] - before["short_2epsilon"][index],
                                   aida_input_seconds=before["aida_seconds"],
                                   lifetime_classification_seconds=before["classification_seconds"])
                    print(f"  Starting {percentage:g}% pruning (epsilon={epsilon:.8g}, deadline={args.pruning_timeout:g}s)", flush=True)
                    with tempfile.TemporaryDirectory(prefix="pruning-experiment-") as temporary:
                        output = Path(temporary) / "pruned.scc"
                        pruning, status, wall, error = run_worker(worker, ["prune", scc, epsilon, output], args.pruning_timeout)
                        row.update(pruning_status=status, pruning_wall_seconds=wall)
                        row.update(pruning)
                        if error:
                            errors.append("pruning: " + error)
                        if status == "ok":
                            after, status, wall, error = run_worker(worker, ["components", output], args.aida_timeout)
                            row.update(aida_output_status=status, aida_output_wall_seconds=wall)
                            if status == "ok":
                                row.update(components_after=after["components"], aida_output_seconds=after["aida_seconds"])
                            if error:
                                errors.append("AIDA output: " + error)
                row["error"] = " | ".join(errors)
                failed += any(row[key] != "ok" for key in ["inspection_status", "pruning_status", "aida_input_status", "aida_output_status"])
                writer.writerow(row)
                stream.flush()  # Preserve completed rows if a later computation times out or is interrupted.
                recorded.add((sample["scc"], input_hash, worker_hash, percentage))
                print(f"  {percentage:g}%: pruning={row['pruning_status']} "
                      f"{row.get('pruning_seconds', '')}s; components "
                      f"{row.get('components_before', '?')} -> {row.get('components_after', '?')}", flush=True)
    print(f"Results: {table}; {failed} new rows with a failed, skipped, or timed-out phase in this run.", flush=True)
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())

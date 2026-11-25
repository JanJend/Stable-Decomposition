#!/usr/bin/env python3
"""Compare explicit and new End_d_0 quotients on precomputed inputs of increasing size."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess

from benchmark_end_sizes import EPSILON_PERCENTAGES, MEMORY_MB, TIMEOUT_SECONDS, measure
from end_quotient_size import FAMILIES
from support import DATA, HERE, ROOT, RUNTIME_SIZES, input_lock, select_samples
from support.runtime import digest
from support.scaling import presentation_info

METHODS = ("explicit", "end-d-0")
DIMENSIONS = ("End_2epsilon_dimension", "quotient_dimension", "image_dimension")
CSV_COLUMNS = ("family", "dataset", "target_points", "n_points", "homology", "epsilon_percent",
               "epsilon_x", "epsilon_y", "dimensions_match",
               "explicit_status", "end_d_0_status", "explicit_runs", "end_d_0_runs",
               "explicit_seconds", "end_d_0_seconds", "speedup_explicit_over_new",
               "explicit_end0_seconds", "explicit_with_E_seconds", "speedup_with_E_over_new",
               "explicit_wall_seconds", "end_d_0_wall_seconds", "explicit_peak_rss_mb", "end_d_0_peak_rss_mb",
               *(f"{prefix}_{key}" for prefix in ("explicit", "end_d_0") for key in DIMENSIONS))


def select_inputs(manifest, sizes, families):
    with manifest.open(newline="") as stream:
        selected = select_samples(list(csv.DictReader(stream)), sizes)
    homology = {family: degree for family, _, degree in FAMILIES}
    inputs = []
    for n in sorted(set(sizes)):
        for family in families:
            matches = [row for row in selected if row["family"] == family and int(row["target_points"]) == n
                       and row["homology"] == homology[family] and (family == "immuno" or row["seed"] == "1")]
            if len(matches) != 1:
                raise ValueError(f"Expected one {family} H{homology[family]} input at target size {n} (synthetic sample 1)")
            inputs.append(dict(matches[0], input=str(manifest.parent / matches[0]["scc"])))
    return inputs


def summarize(records):
    if not records:
        return dict(status="pending", runs=0)
    if any(row["status"] != "ok" for row in records):
        return dict(status=next(row["status"] for row in records if row["status"] != "ok"), runs=len(records))
    if len({tuple(row[key] for key in DIMENSIONS) for row in records}) != 1:
        return dict(status="inconsistent", runs=len(records))
    return dict(status="ok", runs=len(records), **{key: records[0][key] for key in DIMENSIONS},
                **{key: statistics.median(row[key] for row in records)
                   for key in ("seconds", "end0_seconds", "with_E_seconds", "wall_seconds")},
                peak_rss_mb=max(row["peak_rss_mb"] for row in records))


def comparison(sample, percentage, records):
    row = {key: sample[key] for key in CSV_COLUMNS[:5]}
    row.update(epsilon_percent=percentage, epsilon_x=sample["bounding_box_sides"][0] * percentage / 100,
               epsilon_y=sample["bounding_box_sides"][1] * percentage / 100)
    summaries = [summarize(records[method]) for method in METHODS]
    for prefix, summary in zip(("explicit", "end_d_0"), summaries):
        for key, value in summary.items():
            row[f"{prefix}_{key}"] = value
    old, new = summaries
    row["dimensions_match"] = "incomplete"
    if old["status"] == new["status"] == "ok":
        match = all(old[key] == new[key] for key in DIMENSIONS)
        row["dimensions_match"] = "yes" if match else "NO"
        if match and new["seconds"] > 0:
            row.update(speedup_explicit_over_new=old["seconds"] / new["seconds"],
                       speedup_with_E_over_new=old["with_E_seconds"] / new["seconds"])
    return row


def write_csv(path, rows):
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=RUNTIME_SIZES,
                        help="Target sizes, processed in increasing order (default: 20 through 2560, doubling)")
    parser.add_argument("--families", nargs="+", choices=[family for family, _, _ in FAMILIES],
                        default=[family for family, _, _ in FAMILIES])
    parser.add_argument("--percentages", type=float, nargs="+", default=EPSILON_PERCENTAGES)
    parser.add_argument("--repeats", type=int, default=3,
                        help="Median of this many pairs; use one pair when either first timing exceeds 5s (default: 3)")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS, help="Deadline per method/process in seconds")
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/local-release/end_quotient_benchmark")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/end_quotient_comparison")
    args = parser.parse_args()
    if (any(n <= 0 for n in args.sizes) or args.repeats < 1
            or any(not math.isfinite(p) or p < 0 for p in args.percentages)
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Sizes/repeats/deadlines must be positive; percentages/memory must be finite and nonnegative")
    args.worker, args.manifest = args.worker.resolve(), args.manifest.resolve()
    if not args.worker.is_file():
        parser.error("Configure with cmake -S . -B build/local-release, then build with "
                     "cmake --build build/local-release --target end_quotient_benchmark -j 2")
    help_text = subprocess.run([str(args.worker), "--help"], check=True, capture_output=True, text=True, timeout=10).stdout
    if "quotient_methods=explicit,end-d-0 hom_algorithm=optimised" not in help_text:
        parser.error("Rebuild end_quotient_benchmark: expected explicit and new A-1/2 quotient methods")
    if not args.manifest.is_file():
        parser.error(f"Missing shared input manifest: {args.manifest}")
    with input_lock(args.manifest.parent):
        try:
            inputs = select_inputs(args.manifest, args.sizes, list(dict.fromkeys(args.families)))
            for sample in inputs:
                path = Path(sample["input"])
                if not path.is_file() or sample.get("generation_status", "ok") != "ok":
                    raise ValueError(f"Missing or ungenerated input: {path}")
                sample["input_sha256"] = digest(path)
                if sample.get("scc_sha256", sample["input_sha256"]) != sample["input_sha256"]:
                    raise ValueError(f"SCC differs from the shared manifest: {path}")
                sample["generators"], sample["relations"], sample["bounding_box_sides"] = presentation_info(path, per_axis=True)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        directory = args.output_dir.resolve()
        run = directory / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        run.mkdir(parents=True)
        percentages = list(dict.fromkeys(args.percentages))
        worker_hash = digest(args.worker)
        (run / "settings.json").write_text(json.dumps(dict(
            inputs=inputs, sizes=sorted(set(args.sizes)), percentages=percentages, repeats=args.repeats,
            timeout_seconds=args.timeout, memory_limit_mb=args.memory_limit_mb, slow_pair_threshold_seconds=5,
            worker=str(args.worker), worker_sha256=worker_hash, manifest=str(args.manifest), manifest_sha256=digest(args.manifest),
            hom_algorithm="optimised", methods=METHODS, epsilon_convention="per-axis bounding box; shift=2*epsilon",
            timing="fresh processes; internal compute time excludes input/setup; explicit total includes End_0",
        ), indent=2) + "\n")
        if (directory / "results.csv").exists():
            shutil.copyfile(directory / "results.csv", run / "previous_results.csv")
        results, failed = [], False
        print(f"Quotient comparison: {len(inputs)} inputs x {len(percentages)} epsilons; archive: {run}", flush=True)
        try:
            with (run / "observations.jsonl").open("a") as observations:
                for sample in inputs:
                    for percentage in percentages:
                        records = {method: [] for method in METHODS}
                        print(f"RUN {sample['dataset']}: target={sample['target_points']}, n={sample['n_points']}, epsilon={percentage:g}%", flush=True)
                        for repeat in range(args.repeats):
                            order = METHODS if (len(results) + repeat) % 2 == 0 else METHODS[::-1]
                            for method in order:
                                if digest(Path(sample["input"])) != sample["input_sha256"] or digest(args.worker) != worker_hash:
                                    raise SystemExit(f"Input or worker changed; completed results are saved in {run}")
                                prefix = run / f"{sample['dataset']}_{percentage:g}pct_{method}_{repeat+1}"
                                stdout, stderr = Path(str(prefix) + ".stdout.txt"), Path(str(prefix) + ".stderr.txt")
                                command = [str(args.worker), sample["input"], str(percentage), method, str(args.memory_limit_mb)]
                                status, wall, error = measure(command, stdout, stderr, args.timeout)
                                metrics = {}
                                try:
                                    lines = [line[7:] for line in stdout.read_text().splitlines() if line.startswith("RESULT ")]
                                    if lines:
                                        metrics = json.loads(lines[-1])
                                    if status == "ok":
                                        if (metrics.get("method") != method or metrics.get("hom_algorithm") != "optimised"
                                                or any(not isinstance(metrics.get(key), int) or metrics[key] < 0 for key in DIMENSIONS)
                                                or metrics["End_2epsilon_dimension"] != metrics["quotient_dimension"] + metrics["image_dimension"]
                                                or any(not isinstance(metrics.get(key), (int, float)) or not math.isfinite(metrics[key]) or metrics[key] < 0
                                                       for key in ("seconds", "end0_seconds", "with_E_seconds", "peak_rss_mb"))):
                                            raise ValueError("Incomplete or inconsistent worker result")
                                except (ValueError, TypeError) as parse_error:
                                    status, error = "error", str(parse_error)
                                record = dict(metrics, method=method, status=status, error=error, wall_seconds=wall,
                                              dataset=sample["dataset"], epsilon_percent=percentage, repeat=repeat+1,
                                              command=command, stdout_log=stdout.name, stderr_log=stderr.name)
                                records[method].append(record)
                                observations.write(json.dumps(record) + "\n")
                                observations.flush()
                            row = comparison(sample, percentage, records)
                            write_csv(run / "results.csv", [*results, row])
                            shutil.copyfile(run / "results.csv", directory / "results.csv")
                            if (row["dimensions_match"] != "yes"
                                    or (repeat == 0 and max(records[method][-1].get("seconds", 0) for method in METHODS) > 5)):
                                break
                        results.append(row)
                        failed |= row["dimensions_match"] != "yes"
                        speed = row.get("speedup_explicit_over_new")
                        print(f"  explicit={row.get('explicit_seconds', row['explicit_status'])}; "
                              f"End_d_0={row.get('end_d_0_seconds', row['end_d_0_status'])}; "
                              f"speedup={f'{speed:.3f}x' if speed is not None else 'unavailable'}; "
                              f"dimensions_match={row['dimensions_match']}", flush=True)
        except KeyboardInterrupt:
            print(f"Stopped. Completed pairs and raw observations are saved in {run}", flush=True)
            return 130
        print(f"Results: {directory / 'results.csv'}", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

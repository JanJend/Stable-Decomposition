#!/usr/bin/env python3
"""Trace I/K generator growth relative to the presentation of X on shared large inputs."""
import argparse
from contextlib import nullcontext
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import time

from benchmark_end_quotient import select_inputs
from benchmark_end_sizes import EPSILON_PERCENTAGES, MEMORY_MB, POINT_COUNTS, TIMEOUT_SECONDS, measure
from benchmark_torus_runtime import parse_trace, write_csv
from end_quotient_size import FAMILIES
from support import DATA, HERE, ROOT, input_lock
from support.runtime import digest
from support.scaling import presentation_info

IDENTIFIERS = ("family", "dataset", "target_points", "n_points", "homology", "epsilon_percent")
SAMPLE_COLUMNS = (*IDENTIFIERS, "input_generators", "input_relations", "presentation_size")
TRACE_COLUMNS = (*SAMPLE_COLUMNS, "object", "stage", "pass", "map", "generators", "elapsed_seconds",
                 "over_x_generators", "over_presentation_size")
PEAK_COLUMNS = tuple(f"{obj.lower()}_peak_generators" for obj in ("I", "K", "J", "L"))
SIZE_COLUMNS = (*PEAK_COLUMNS, "i_peak_reduced_generators", "k_peak_reduced_generators",
                "i_final_generators", "k_final_generators")
RATIO_COLUMNS = tuple(f"{key}_{denominator}" for key in SIZE_COLUMNS
                      for denominator in ("over_x_generators", "over_presentation_size"))
RESULT_COLUMNS = (*SAMPLE_COLUMNS, "status", "error", "peaks_are_lower_bounds", "trace_events",
                  "wall_seconds", "pruning_seconds", "epsilon_x", "epsilon_y", "threads", "fast_i", "peak_rss_mb",
                  "end_dimension", "quotient_dimension", "i_maps", "k_maps", "i_iterations", "k_iterations",
                  "output_generators", "output_relations", *SIZE_COLUMNS, *RATIO_COLUMNS)


def parse_generators(path):
    events = []
    with path.open() as stream:
        for line in stream:
            if not line.startswith("GENERATORS\t") or not line.endswith("\n"):
                continue  # A killed process may leave its last line incomplete.
            _, obj, stage, iteration, map_index, count, elapsed = line.rstrip("\n").split("\t")
            iteration, map_index, count, elapsed = int(iteration), int(map_index), int(count), float(elapsed)
            if (obj not in ("I", "K", "J", "L")
                    or stage not in ("initial", "raw", "reduced", "pass_end", "final")
                    or min(iteration, map_index, count) < 0 or not math.isfinite(elapsed) or elapsed < 0):
                raise ValueError("Invalid generator event")
            events.append(dict(object=obj, stage=stage, **{"pass": iteration, "map": map_index},
                               generators=count, elapsed_seconds=elapsed))
    return events


def ratios(count, sample):
    # An empty denominator has no defined ratio, including 0/0.
    return dict(over_x_generators=count / sample["input_generators"] if sample["input_generators"] else "",
                over_presentation_size=count / sample["presentation_size"] if sample["presentation_size"] else "")


def size_summary(events, sample):
    result = {}
    for obj in ("I", "K", "J", "L"):
        rows = [row for row in events if row["object"] == obj]
        if not rows:
            continue
        result[f"{obj.lower()}_peak_generators"] = max(row["generators"] for row in rows)
        if obj in ("I", "K"):
            reduced = [row["generators"] for row in rows if row["stage"] != "raw"]
            if reduced:
                result[f"{obj.lower()}_peak_reduced_generators"] = max(reduced)
            final = [row["generators"] for row in rows if row["stage"] == "final"]
            if final:
                result[f"{obj.lower()}_final_generators"] = final[-1]
    for key, count in list(result.items()):
        result.update({f"{key}_{denominator}": ratio for denominator, ratio in ratios(count, sample).items()})
    return result


def validate_complete(metrics, events, sample, fast_i):
    if (metrics.get("input_generators") != sample["input_generators"]
            or metrics.get("input_relations") != sample["input_relations"]
            or metrics.get("fast_i") != fast_i):
        raise ValueError("Worker input dimensions or I update rule disagree with requested settings")
    for key in ("pruning_seconds", "epsilon_x", "epsilon_y", "threads", "peak_rss_mb", "end_dimension",
                "quotient_dimension", "i_maps", "k_maps", "i_iterations", "k_iterations",
                "output_generators", "output_relations"):
        value = metrics.get(key)
        if not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
            raise ValueError(f"Missing or invalid worker metric: {key}")
        if key not in ("pruning_seconds", "epsilon_x", "epsilon_y", "peak_rss_mb") and type(value) is not int:
            raise ValueError(f"Expected an integer worker metric: {key}")
    for obj, initial in (("I", sample["input_generators"]), ("K", 0)):
        rows = [row for row in events if row["object"] == obj]
        initials = [row for row in rows if row["stage"] == "initial"]
        finals = [row for row in rows if row["stage"] == "final"]
        iterations = metrics[obj.lower() + "_iterations"]
        if (len(initials) != 1 or initials[0]["generators"] != initial or initials[0]["pass"] != 0
                or len(finals) != 1 or finals[0]["pass"] != iterations or rows[-1] != finals[0]
                or [row["pass"] for row in rows if row["stage"] == "pass_end"] != list(range(1, iterations + 1))):
            raise ValueError(f"Incomplete {obj} generator trace")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=POINT_COUNTS, help="Target points (default: 640 1280)")
    parser.add_argument("--families", nargs="+", choices=[family for family, _, _ in FAMILIES],
                        default=[family for family, _, _ in FAMILIES])
    parser.add_argument("--percentages", type=float, nargs="+", default=EPSILON_PERCENTAGES)
    parser.add_argument("--input", type=Path, help="Use one SCC presentation instead of the shared manifest")
    parser.add_argument("--fast-i", action="store_true", help="Refresh I's target per map (default: fixed per pass)")
    parser.add_argument("--threads", type=int, default=0, help="Worker budget; 0 detects automatically")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS)
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    parser.add_argument("--worker", type=Path, default=ROOT / "build/pruning_torus_runtime")
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/generator_sizes")
    args = parser.parse_args()
    if (any(n <= 0 for n in args.sizes) or args.threads < 0 or args.threads > 2147483647
            or any(not math.isfinite(p) or p < 0 for p in args.percentages)
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Invalid sizes, percentages, threads or resource limits")
    args.worker, args.manifest = args.worker.resolve(), args.manifest.resolve()
    if not args.worker.is_file():
        parser.error("Build first: cmake --build build --target pruning_torus_runtime -j 2")
    help_text = subprocess.run([str(args.worker), "--help"], check=True, capture_output=True, text=True, timeout=10).stdout
    if "Generator size trace v1" not in help_text or "Column storage: compact" not in help_text:
        parser.error("Rebuild pruning_torus_runtime with compact columns and generator tracing")
    if not args.input and not args.manifest.is_file():
        parser.error(f"Missing shared manifest: {args.manifest}")
    with nullcontext() if args.input else input_lock(args.manifest.parent):
        try:
            inputs = ([dict(family="custom", dataset=args.input.stem, input=str(args.input.resolve()),
                            target_points="", n_points="", homology="")]
                      if args.input else select_inputs(args.manifest, args.sizes, list(dict.fromkeys(args.families))))
            for sample in inputs:
                path = Path(sample["input"])
                if sample.get("generation_status", "ok") != "ok" or not path.is_file():
                    raise ValueError(f"Missing or ungenerated input: {path}")
                sample["input_sha256"] = digest(path)
                if sample.get("scc_sha256", sample["input_sha256"]) != sample["input_sha256"]:
                    raise ValueError(f"SCC differs from shared manifest: {path}")
                sample["input_generators"], sample["input_relations"], _ = presentation_info(path, per_axis=True)
                sample["presentation_size"] = sample["input_generators"] + sample["input_relations"]
        except (OSError, ValueError) as error:
            parser.error(str(error))
        directory = args.output_dir.resolve()
        run = directory / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        run.mkdir(parents=True)
        percentages = list(dict.fromkeys(args.percentages))
        worker_hash = digest(args.worker)
        (run / "settings.json").write_text(json.dumps(dict(
            inputs=inputs, percentages=percentages, threads=args.threads, fast_i=args.fast_i,
            timeout_seconds=args.timeout, memory_limit_mb=args.memory_limit_mb,
            worker=str(args.worker), worker_sha256=worker_hash, platform=platform.platform(),
            manifest=None if args.input else str(args.manifest),
            manifest_sha256=None if args.input else digest(args.manifest),
            configuration="production module pruning; default optimizations; optimised A-1/2; compact columns",
            epsilon_convention="per-axis bounding box; shift=2*epsilon",
            counts="Stored columns: initial, before/after I/K minimization, pass ends and final. "
                   "J=eta(I) and L=eta(K)+sum f(K) are separate auxiliaries; L includes unreduced batches. "
                   "Peak reduced counts exclude raw events. Initial counts use X's supplied generators. "
                   "No additional minimizations/presentations; internal kernel workspace is not measured.",
            ratios="Count / generators(X), and count / (generators(X)+relations(X)); empty denominator gives blank",
        ), indent=2) + "\n")
        cache = args.worker.parent / "CMakeCache.txt"
        if cache.is_file():
            shutil.copyfile(cache, run / cache.name)
        for name in ("results.csv", "generators.csv"):
            if (directory / name).is_file():
                shutil.copyfile(directory / name, run / ("previous_" + name))
        results, traces, failed = [], [], False
        print(f"Generator sizes: {len(inputs)} inputs x {len(percentages)} epsilons; archive: {run}", flush=True)
        for index, sample in enumerate(inputs, 1):
            for epsilon_index, percentage in enumerate(percentages, 1):
                if digest(args.worker) != worker_hash or digest(Path(sample["input"])) != sample["input_sha256"]:
                    raise SystemExit(f"Input or worker changed; completed observations saved in {run}")
                prefix = f"input_{index}_epsilon_{epsilon_index}_{percentage:g}"
                stdout, stderr = run / (prefix + ".trace.tsv"), run / (prefix + ".stderr.txt")
                command = [str(args.worker), sample["input"], str(percentage), str(args.threads),
                           str(args.memory_limit_mb), "--generator-trace"]
                if not args.fast_i:
                    command.append("--standard-i")
                print(f"START {sample['dataset']}, n={sample['n_points']}, epsilon={percentage:g}%", flush=True)
                started = time.monotonic()
                try:
                    status, wall, error = measure(command, stdout, stderr, args.timeout)
                except KeyboardInterrupt:
                    status, wall, error = "interrupted", time.monotonic() - started, "Interrupted by user"
                if status == "error" and error == "Worker exited with status -9":
                    status, error = "killed_sigkill", "Worker terminated by SIGKILL; cause not reported"
                metrics, events = {}, []
                try:
                    metrics, _ = parse_trace(stdout)
                    events = parse_generators(stdout)
                    if status == "ok":
                        validate_complete(metrics, events, sample, args.fast_i)
                except (OSError, ValueError) as parse_error:
                    status, error = "error", f"Invalid worker trace: {parse_error}"
                identity = {key: sample[key] for key in SAMPLE_COLUMNS if key != "epsilon_percent"}
                identity["epsilon_percent"] = percentage
                row = dict(metrics, **size_summary(events, sample))
                row.update(identity, status=status, error=error, peaks_are_lower_bounds=status != "ok",
                           trace_events=len(events), wall_seconds=wall)
                results.append(row)
                traces.extend(dict(event, **identity, **ratios(event["generators"], sample)) for event in events)
                with (run / "observations.jsonl").open("a") as stream:
                    stream.write(json.dumps(dict(row, command=command, stdout=stdout.name, stderr=stderr.name)) + "\n")
                for name, columns, rows in (("results", RESULT_COLUMNS, results), ("generators", TRACE_COLUMNS, traces)):
                    write_csv(run / (name + ".csv"), columns, rows)
                    shutil.copyfile(run / (name + ".csv"), directory / (name + ".csv"))
                print(f"DONE {status}: peak I={row.get('i_peak_generators', '?')}, "
                      f"peak K={row.get('k_peak_generators', '?')}, "
                      f"X={sample['input_generators']} generators + {sample['input_relations']} relations", flush=True)
                failed |= status != "ok"
                if status == "interrupted":
                    return 130
        print(f"Results: {directory}; see results.csv and generators.csv", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

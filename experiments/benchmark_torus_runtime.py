#!/usr/bin/env python3
"""Time every pruning operation on the precomputed 1280-point Gaussian torus H1."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import time

from benchmark_end_sizes import measure
from benchmark_runtime import EPSILON_PERCENTAGES, MEMORY_MB
from support import DATA, HERE, ROOT, input_lock
from support.runtime import digest
from support.generator_options import add_generator_arguments, generator_flags, generator_settings

TIMEOUT_SECONDS = 1800
GENERATOR_DEFAULTS = {"composition-certificate-visits": 262144, "composition-nilpotent-products": 262144}

STEP_COLUMNS = ("epsilon_percent", "repeat", "operation", "pass", "call", "start_seconds", "end_seconds", "seconds", "completed")
TIMING_COLUMNS = ("epsilon_percent", "repeat", "operation", "calls", "unfinished_calls", "seconds")
RESULT_COLUMNS = ("dataset", "n_points", "epsilon_percent", "repeat", "status", "error", "wall_seconds",
                  "pruning_seconds", "epsilon_x", "epsilon_y", "threads", "fast_i", "parallel_i", "parallel_k", "minimize_i_target", "peak_rss_mb", "end_dimension",
                  "quotient_dimension", "i_maps", "k_maps", "i_iterations", "k_iterations",
                  "output_generators", "output_relations")


def select_input(manifest):
    with manifest.open(newline="") as stream:
        matches = [row for row in csv.DictReader(stream)
                   if row["family"] == "torus_gauss" and row["target_points"] == "1280"
                   and row["n_points"] == "1280" and row["homology"] == "1" and row["seed"] == "1"]
    if len(matches) != 1:
        raise ValueError("Expected exactly one 1280-point Gaussian torus H1, synthetic sample 1")
    sample = matches[0]
    scc = (manifest.parent / sample["scc"]).resolve()
    if sample["generation_status"] != "ok" or not scc.is_file() or digest(scc) != sample["scc_sha256"]:
        raise ValueError(f"Missing, ungenerated or modified torus input: {scc}")
    return dict(sample, input=str(scc))


def parse_trace(path):
    """Pair events by operation and call; retain unfinished timers after a kill."""
    steps, metrics = {}, {}
    with path.open() as stream:
        for line in stream:
            if not line.endswith("\n"):  # A killed process can leave a truncated last event.
                continue
            if line.startswith("RESULT "):
                metrics.update(json.loads(line[7:]))
            elif line.startswith("METRICS "):
                metrics.update(json.loads(line[8:]))
            elif line.startswith("STEP\t"):
                _, event, operation, iteration, call, timestamp, seconds = line.strip().split("\t")
                key = operation, int(call)
                timestamp, seconds = float(timestamp), float(seconds)
                if any(not math.isfinite(x) or x < 0 for x in (timestamp, seconds)):
                    raise ValueError("Invalid trace duration")
                if event == "start" and key not in steps:
                    steps[key] = dict(operation=operation, **{"pass": int(iteration)}, call=int(call),
                                      start_seconds=timestamp, end_seconds="", seconds="", completed=False)
                elif event == "end" and key in steps and not steps[key]["completed"]:
                    row = steps[key]
                    if timestamp < row["start_seconds"] or int(iteration) != row["pass"]:
                        raise ValueError("Inconsistent trace event")
                    row.update(end_seconds=timestamp, seconds=seconds, completed=True)
                else:
                    raise ValueError(f"Unpaired or duplicate trace event: {line.strip()}")
    return metrics, sorted(steps.values(), key=lambda row: row["start_seconds"])


def aggregate_steps(steps):
    groups = {}
    for step in steps:
        key = step["epsilon_percent"], step["repeat"], step["operation"]
        row = groups.setdefault(key, dict(zip(TIMING_COLUMNS[:3], key), calls=0, unfinished_calls=0, seconds=0.0))
        row["calls" if step["completed"] else "unfinished_calls"] += 1
        if step["completed"]:
            row["seconds"] += step["seconds"]
    return list(groups.values())


def write_csv(path, columns, rows):
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--percentages", type=float, nargs="+", default=EPSILON_PERCENTAGES)
    parser.add_argument("--repeats", type=int, default=1, help="Fresh processes per epsilon (default: 1; this is a slow input)")
    parser.add_argument("--threads", type=int, default=0, help="Worker budget; 0 detects automatically (default)")
    parser.add_argument("--standard-i", action="store_true", help="Compare fixed-target I with the default fast I")
    parser.add_argument("--parallel-i", action="store_true", help="Fix J per pass and process groups of constraints concurrently")
    parser.add_argument("--parallel-k", action="store_true", help="Compute and sum groups of K images concurrently")
    parser.add_argument("--minimize-i-target", action="store_true", help="Restore separate J minimization for comparison")
    add_generator_arguments(parser, defaults=GENERATOR_DEFAULTS)
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS)
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    parser.add_argument("--worker", type=Path, default=ROOT / "build/pruning_torus_runtime")
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/torus_runtime_1280")
    args = parser.parse_args()
    if args.parallel_i:
        args.standard_i = True
    generator_options = generator_settings(args)
    if (args.repeats < 1 or args.threads < 0 or args.threads > 2147483647
            or any(not math.isfinite(p) or p < 0 for p in args.percentages)
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Invalid percentages, repetitions, threads or resource limits")
    args.worker, args.manifest = args.worker.resolve(), args.manifest.resolve()
    if not args.worker.is_file():
        parser.error("Build first: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release "
                     "&& cmake --build build --target pruning_torus_runtime -j 2")
    help_text = subprocess.run([str(args.worker), "--help"], capture_output=True, text=True, check=True, timeout=10).stdout
    if "Torus runtime trace v1" not in help_text or "Column storage: compact" not in help_text:
        parser.error("Rebuild pruning_torus_runtime with compact columns")
    if not args.manifest.is_file():
        parser.error(f"Missing shared manifest: {args.manifest}; prepare the shared data first")
    with input_lock(args.manifest.parent):
        try:
            sample = select_input(args.manifest)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        directory = args.output_dir.resolve()
        run = directory / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        run.mkdir(parents=True)
        worker_hash = digest(args.worker)
        (run / "settings.json").write_text(json.dumps(dict(
            input=sample, manifest=str(args.manifest), manifest_sha256=digest(args.manifest),
            worker=str(args.worker), worker_sha256=worker_hash, percentages=args.percentages,
            repeats=args.repeats, threads=args.threads, timeout_seconds=args.timeout,
            memory_limit_mb=args.memory_limit_mb, platform=platform.platform(), machine=platform.machine(),
            configuration="all-optimizations; " + ("parallel fixed-target I" if args.parallel_i else "standard I" if args.standard_i else "fast I") + "; optimised A-1/2; compact columns; coefficient map families; batched restrictions",
            parallel_i=args.parallel_i, parallel_k=args.parallel_k, minimize_i_target=args.minimize_i_target,
            generator_options=generator_options,
            epsilon_convention="per-axis bounding box; shift=2*epsilon",
            timing="Internal wall seconds, including tracing overhead; load/sort separate. "
                   "total includes basis, I/K and finalization; basis includes end0/hom_quotient/composition. "
                   "map_lifts is coefficient-family setup; active maps are lifted inside loop phases. "
                   "I/K phases overlap with multiple workers; parallel I/K operation times sum worker durations. "
                   "Nested/overlapping rows must not be added.",
        ), indent=2) + "\n")
        cache = args.worker.parent / "CMakeCache.txt"
        if cache.is_file():
            shutil.copyfile(cache, run / "CMakeCache.txt")
        results, steps, failed = [], [], False
        print(f"1280-point torus runtime trace; archive: {run}", flush=True)
        for percentage in dict.fromkeys(args.percentages):
            for repeat in range(1, args.repeats + 1):
                if digest(args.worker) != worker_hash or digest(Path(sample["input"])) != sample["scc_sha256"]:
                    raise SystemExit(f"Worker or input changed; completed observations saved in {run}")
                prefix = f"epsilon_{percentage:g}_repeat_{repeat}"
                stdout, stderr = run / (prefix + ".trace.tsv"), run / (prefix + ".stderr.txt")
                command = [str(args.worker), sample["input"], str(percentage), str(args.threads), str(args.memory_limit_mb)]
                command.append("--output-scc=" + str(run / (prefix + ".output.scc")))
                command.extend(generator_flags(generator_options))
                if args.standard_i:
                    command.append("--standard-i")
                if args.parallel_i:
                    command.append("--parallel-i")
                if args.parallel_k:
                    command.append("--parallel-k")
                if args.minimize_i_target:
                    command.append("--minimize-i-target")
                print(f"START epsilon={percentage:g}%, repeat={repeat}", flush=True)
                started = time.monotonic()
                try:
                    status, wall, error = measure(command, stdout, stderr, args.timeout)
                except KeyboardInterrupt:
                    status, wall, error = "interrupted", time.monotonic() - started, "Interrupted"
                if status == "error" and error == "Worker exited with status -9":
                    status, error = "killed_sigkill", "Worker terminated by SIGKILL; cause not reported"
                try:
                    metrics, measured_steps = parse_trace(stdout)
                except (OSError, ValueError) as parse_error:
                    metrics, measured_steps = {}, []
                    status, error = "error", f"Invalid worker trace: {parse_error}"
                total = next((step for step in measured_steps if step["operation"] == "total"), None)
                if status == "ok" and (any(not isinstance(metrics.get(key), (int, float))
                                           or not math.isfinite(metrics[key]) or metrics[key] < 0
                                           for key in RESULT_COLUMNS[7:])
                                       or any(not step["completed"] for step in measured_steps)
                                       or not total or metrics["pruning_seconds"] != total["seconds"]):
                    status, error = "error", "Incomplete worker trace/result"
                row = dict(metrics, dataset=sample["dataset"], n_points=1280, epsilon_percent=percentage,
                           repeat=repeat, status=status, error=error, wall_seconds=wall)
                results.append(row)
                steps.extend(dict(step, epsilon_percent=percentage, repeat=repeat) for step in measured_steps)
                with (run / "observations.jsonl").open("a") as stream:
                    stream.write(json.dumps(dict(row, command=command, stdout=stdout.name, stderr=stderr.name)) + "\n")
                for name, columns, rows in (("results", RESULT_COLUMNS, results), ("steps", STEP_COLUMNS, steps),
                                            ("timings", TIMING_COLUMNS, aggregate_steps(steps))):
                    write_csv(run / (name + ".csv"), columns, rows)
                    shutil.copyfile(run / (name + ".csv"), directory / (name + ".csv"))
                duration = f"{metrics['pruning_seconds']:.3f}s" if "pruning_seconds" in metrics else "incomplete"
                print(f"DONE epsilon={percentage:g}%: {status}; pruning={duration}", flush=True)
                failed |= status != "ok"
                if status == "interrupted":
                    return 130
        print(f"Results: {directory}; see steps.csv, timings.csv and results.csv", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

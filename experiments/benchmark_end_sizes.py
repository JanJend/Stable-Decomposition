#!/usr/bin/env python3
"""Measure Hom reductions for Gaussian-torus H1; write results/End_sizes/results.csv."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import subprocess
import time

from support import DATA, HERE, ROOT
from support.runtime import digest
from support.scaling import presentation_info
from support.generator_options import add_generator_arguments, generator_settings, generator_flags

# One 640-point Gaussian torus (sample 1, RNG seed 42), with percentage epsilons.
INPUT = DATA / "presentations/torus_gauss/torus_gauss_640_1_H1.scc"
EPSILON_PERCENTAGES = [1, 2, 4, 8]
TIMEOUT_SECONDS = 10 * 60
MEMORY_MB = 12288
METHODS = (("sparse_nilpotent", "Sparse-nilpotent"),
           ("complement_orbit", "Complement-orbit"), ("combined", "Combined"))
SPACES = (("full", "Full"), ("quotient", "Quotient"))
METHOD_STAGES = tuple(f"{space}_{method}" for space, _ in SPACES for method, _ in METHODS)
STAGES = ("End_d_basis", "quotient_basis", "full_coordinates", "quotient_coordinates", *METHOD_STAGES)
# Sizes first, then shared setup and independent method timings (both sides).
COLUMNS = (
    ("Pointset + epsilon", None),
    ("Generator settings", "generator_settings"),
    ("dim End_2d", "End_d_basis_I_size"),
    ("dim End_2d/End_0", "quotient_basis_I_size"),
    *((f"{label} {name} {side}", f"{space}_{method}_{side}_size")
      for space, label in SPACES for method, name in METHODS for side in ("I", "K")),
    ("End_2d time (s)", "End_d_basis_seconds"),
    ("End_0 time (s)", "End_0_seconds"),
    ("Quotient time (s)", "quotient_basis_seconds"),
    *((f"{label} coordinates time (s)", f"{space}_coordinates_seconds") for space, label in SPACES),
    *((f"{label} {name} time (s)", f"{space}_{method}_seconds")
      for space, label in SPACES for method, name in METHODS),
)
FIELDS = tuple(label for label, _ in COLUMNS)


def parse_output(output):
    """Keep completed stages even when a later stage fails or times out."""
    result = {}
    for line in output.splitlines():
        values = dict(re.findall(r"(\w+)=([^\s]+)", line))
        if "End_0_dimension" in values:
            result.update(End_0_dimension=int(values["End_0_dimension"]), End_0_seconds=float(values["seconds"]))
        if "epsilon" in values:
            result.update(epsilon=float(values["epsilon"]), shift=float(values["shift"]))
        columns = line.split("\t")
        if len(columns) == 4 and columns[0] in STAGES:
            result.update({f"{columns[0]}_I_size": int(columns[1]),
                           f"{columns[0]}_K_size": int(columns[2]),
                           f"{columns[0]}_seconds": float(columns[3])})
    return result


def measure(command, stdout_path, stderr_path, timeout):
    start = time.monotonic()
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        try:
            while True:
                remaining = timeout - (time.monotonic() - start)
                if remaining <= 0:
                    process.kill()
                    process.wait()
                    return "timeout", time.monotonic() - start, f"Exceeded {timeout:g} seconds"
                try:
                    code = process.wait(timeout=min(30, remaining))
                    status = "ok" if code == 0 else "memory_limit" if code == 75 else "error"
                    return status, time.monotonic() - start, "" if code == 0 else f"Worker exited with status {code}"
                except subprocess.TimeoutExpired:
                    lines = stdout_path.read_text().splitlines()
                    print(f"  {time.monotonic()-start:.0f}s: {lines[-1] if lines else 'preparing input'}", flush=True)
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--percentages", type=float, nargs="+", default=EPSILON_PERCENTAGES,
                        help="Percentages of the largest generator/relation degree extent (default: 1 2 4 8)")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/local-release/hom_reduction_benchmark")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/End_sizes")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS,
                        help="Seconds per epsilon (default: 600, i.e. 10 minutes)")
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    parser.add_argument("--sparse-nilpotent-generators", action=argparse.BooleanOptionalAction, default=True,
                        help="Measure sparse-nilpotent compression (default: enabled)")
    parser.add_argument("--complement-orbit-generators", action=argparse.BooleanOptionalAction, default=True,
                        help="Measure complement-orbit compression (default: enabled)")
    add_generator_arguments(parser)
    args = parser.parse_args()
    search_settings = generator_settings(args)
    enabled = {"sparse_nilpotent": args.sparse_nilpotent_generators,
               "complement_orbit": args.complement_orbit_generators,
               "combined": args.sparse_nilpotent_generators and args.complement_orbit_generators}
    disabled_fields = {f"{space}_{method}_{suffix}"
                       for space, _ in SPACES for method, _ in METHODS if not enabled[method]
                       for suffix in ("I_size", "K_size", "seconds")}
    method_flags = [f"--{'' if active else 'no-'}{method.replace('_', '-')}-generators"
                    for method, active in enabled.items() if method != "combined"] + generator_flags(search_settings)
    if (any(not math.isfinite(p) or p < 0 for p in args.percentages)
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Percentages/memory must be finite and nonnegative; timeout must be finite and positive")
    args.input, args.worker = args.input.resolve(), args.worker.resolve()
    if not args.input.is_file():
        parser.error(f"Missing SCC input: {args.input}")
    if not args.worker.is_file():
        parser.error("Build the worker first: cmake --build build/local-release --target hom_reduction_benchmark")
    generators, relations, extent = presentation_info(args.input)
    input_hash, worker_hash = digest(args.input), digest(args.worker)
    run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    directory = args.output_dir.resolve()
    logs = directory / run
    logs.mkdir(parents=True)
    (logs / "settings.json").write_text(json.dumps(dict(
        input=str(args.input), input_sha256=input_hash, worker=str(args.worker), worker_sha256=worker_hash,
        percentages=args.percentages, bounding_box_extent=extent, adaptive=False, methods=enabled, generator_settings=search_settings,
        timeout_seconds=args.timeout, memory_limit_mb=args.memory_limit_mb), indent=2) + "\n")
    results = directory / "results.csv"
    # Preserve older tables if their column layout differs from this experiment.
    if results.exists() and results.stat().st_size:
        with results.open(newline="") as stream:
            if next(csv.reader(stream)) != list(FIELDS):
                results.rename(logs / "previous_results.csv")
    failed = False
    with results.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        if stream.tell() == 0:
            writer.writeheader()
        for percentage in dict.fromkeys(args.percentages):
            prefix = logs / f"epsilon_{percentage}pct"
            stdout_path, stderr_path = Path(str(prefix) + ".stdout.txt"), Path(str(prefix) + ".stderr.txt")
            command = [str(args.worker), str(args.input), str(percentage), "--percent",
                       "--memory-limit-mb", str(args.memory_limit_mb), *method_flags]
            row = dict(run=run, input=str(args.input), input_sha256=input_hash, worker_sha256=worker_hash,
                       epsilon_percent=percentage, epsilon=extent * (percentage / 100),
                       shift=2 * extent * (percentage / 100), bounding_box_extent=extent,
                       generators=generators, relations=relations, methods=enabled,
                       generator_settings=search_settings,
                       stdout_log=str(stdout_path.relative_to(directory)), stderr_log=str(stderr_path.relative_to(directory)))
            print(f"RUN {args.input.stem}: epsilon={percentage:g}% ({row['epsilon']:.8g})", flush=True)
            status, seconds, error = measure(command, stdout_path, stderr_path, args.timeout)
            try:
                row.update(parse_output(stdout_path.read_text()))
            except (ValueError, KeyError) as parse_error:
                status, error = "error", f"Cannot parse worker output: {parse_error}"
            if status == "ok" and any(key not in row for _, key in COLUMNS[1:] if key not in disabled_fields):
                status, error = "error", "Worker output is incomplete"
            row.update(status=status, process_seconds=seconds, error=error)
            (Path(str(prefix) + ".json")).write_text(json.dumps(dict(row, command=command), indent=2) + "\n")
            name = re.sub(r"_H\d+$", "", args.input.stem)
            table_row = {label: "disabled" if key in disabled_fields else row.get(key, status) for label, key in COLUMNS[1:]}
            table_row["Generator settings"] = json.dumps(search_settings, sort_keys=True, separators=(",", ":"))
            table_row[FIELDS[0]] = f"{name} — {percentage:g}%"
            writer.writerow(table_row)
            stream.flush()
            failed |= status != "ok"
            print(f"  {status}: {seconds:.2f}s; full Hom={row.get('End_d_basis_I_size', '?')}, "
                  f"quotient={row.get('quotient_basis_I_size', '?')}, "
                  f"combined quotient I/K={row.get('quotient_combined_I_size', 'disabled')}/"
                  f"{row.get('quotient_combined_K_size', 'disabled')}", flush=True)
    print(f"Results: {results}", flush=True)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

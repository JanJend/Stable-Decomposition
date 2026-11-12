"""Shared CLI execution, provenance, resume, and reporting for paper experiments."""
import csv
import datetime
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time

from . import HERE, ROOT
from .runtime import digest
from .scaling import SWITCHES, presentation_info
from .paper_tables import VARIANTS, write_reports


def column_storage(pruning):
    help_text = subprocess.run([str(pruning), "--help"], check=True, capture_output=True,
                               text=True, timeout=10).stdout
    match = re.search(r"^Column storage: (compact|vector)$", help_text, re.MULTILINE)
    if not match:
        raise ValueError(f"Rebuild {pruning}: its --help does not identify column storage")
    return match[1]


def execute(pruning, scc, epsilon, variant, timeout, memory_mb, expected_storage=None):
    flags = [] if variant == "baseline" else ["--" + variant]
    command = [str(pruning), str(scc), "--epsilon", repr(epsilon), "--profile", "--no-output",
               "--memory-limit-mb", str(memory_mb), *flags]
    start = time.monotonic()
    with tempfile.TemporaryFile(mode="w+b") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        status = "error"
        try:
            while True:
                remaining = timeout - (time.monotonic() - start)
                if remaining <= 0:
                    process.kill()
                    process.wait()
                    status = "timeout"
                    break
                try:
                    process.wait(timeout=min(30, remaining))
                    status = ("ok" if process.returncode == 0 else "memory_limit" if process.returncode == 75
                              else "killed_sigkill" if process.returncode == -9 else "error")
                    break
                except subprocess.TimeoutExpired:
                    print(f"RUNNING {scc.stem} {variant}: {time.monotonic()-start:.0f}/{timeout:g}s", flush=True)
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()
        elapsed = time.monotonic() - start
        log.seek(0, os.SEEK_END)
        log.seek(max(0, log.tell() - 65536))
        output = log.read().decode(errors="replace")
    result = dict(status=status, process_seconds=elapsed, returncode=process.returncode,
                  command=command, error="")
    if status == "ok":
        total = re.search(r"^Total\s+1\s+([0-9.]+)\s*$", output, re.MULTILINE)
        iterations = re.search(r"I iterations: (\d+); K iterations: (\d+)", output)
        options = re.search(r"Module pruning options \(0=off, 1=on\):([^\n]+)", output)
        storage = re.search(r"^Column storage: (compact|vector)$", output, re.MULTILINE)
        expected = {name.replace("-", "_"): int(variant == "all-optimizations" or variant == name)
                    for name in SWITCHES}
        actual = dict((name, int(value)) for name, value in re.findall(r"(\w+)=(\d)", options[1])) if options else {}
        if (not total or not iterations or actual != expected or not storage or
                (expected_storage is not None and storage[1] != expected_storage)):
            result.update(status="error", error="Missing profile or unexpected optimization switches/column storage: " + output[-2000:])
        else:
            result.update(pruning_seconds=float(total[1]), i_iterations=int(iterations[1]),
                          k_iterations=int(iterations[2]), column_storage=storage[1])
    else:
        result["error"] = output[-2000:]
    return result


def main(kind, *, points, epsilon_percents, repeats=3, timeout=600, memory_mb=12288,
         repeat_below=5, manifest=None, pruning=None, output_dir=None):
    """Run the fixed experiment; resume matching results and preserve older runs."""
    points = sorted(set(points))
    percentages = sorted(set(epsilon_percents))
    if (kind not in ("runtime", "optimizations") or not points or any(n < 6 for n in points)
            or not percentages or any(not 0 < p <= 100 for p in percentages)
            or repeats < 1 or not math.isfinite(repeat_below) or repeat_below <= 0
            or not math.isfinite(timeout) or timeout <= 0
            or not math.isfinite(memory_mb) or memory_mb < 0):
        raise ValueError("Invalid experiment settings: check sizes, epsilons, repetitions and resource limits")
    pruning = ROOT / "build/pruning" if pruning is None else Path(pruning).resolve()
    manifest = HERE / "manifest.csv" if manifest is None else Path(manifest).resolve()
    if not pruning.is_file() or not os.access(pruning, os.X_OK):
        raise SystemExit(f"Missing pruning executable: {pruning}. Build it with CMake before running the benchmark.")
    if not manifest.is_file():
        raise SystemExit(f"Prepare the datasets first: python3 {HERE / 'generate_data.py'}")
    storage = column_storage(pruning)
    if storage != "compact":
        raise SystemExit("The paper experiments use compact columns; rebuild with -DGRLINA_USE_COMPACT_COLUMNS=ON")
    binary_hash = digest(pruning)
    with manifest.open(newline="") as stream:
        samples = [row for row in csv.DictReader(stream) if int(row["target_points"]) in points]
    samples.sort(key=lambda row: (int(row["n_points"]), row["family"], int(row["homology"]), row["seed"], row["dataset"]))
    if not samples:
        raise ValueError("No manifest inputs match POINT_COUNTS")
    variants = ["all-optimizations"] if kind == "runtime" else VARIANTS
    config = dict(kind=kind, manifest=str(manifest), manifest_sha256=digest(manifest),
                  pruning=str(pruning), pruning_sha256=binary_hash, column_storage=storage,
                  epsilon_percents=percentages, repeats=repeats, repeat_below=repeat_below,
                  timeout=timeout, memory_mb=memory_mb, points=points, variants=variants,
                  machine=platform.machine(), platform=platform.platform(), python=sys.version)
    directory = Path(output_dir).resolve() if output_dir is not None else HERE / "paper_results" / kind
    metadata = directory / "run.json"
    if directory.exists() and any(directory.iterdir()):
        candidates = [directory, *sorted(directory.parent.glob(directory.name + "-*"), reverse=True)]
        for candidate in candidates:
            saved = candidate / "run.json"
            if saved.is_file() and json.loads(saved.read_text())["configuration"] == config:
                directory, metadata = candidate, saved
                print(f"Continuing saved experiment: {directory}", flush=True)
                break
        else:
            # Keep historical measurements intact when the code or settings change.
            stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            directory = directory.with_name(directory.name + "-" + stamp)
            metadata = directory / "run.json"
            print(f"Code/settings changed; preserving earlier results. New results: {directory}", flush=True)
    directory.mkdir(parents=True, exist_ok=True)
    if not metadata.exists():
        metadata.write_text(json.dumps(dict(configuration=config, created_utc=datetime.datetime.now(
            datetime.timezone.utc).isoformat(), timings="Internal pruning profile and complete process wall time; seconds",
            repetition_policy="Up to three by default; first failures and configurations above repeat-below get one observation"), indent=2) + "\n")
        cache = pruning.parent / "CMakeCache.txt"
        if cache.exists():
            shutil.copyfile(cache, directory / f"build-cache-{storage}.txt")
    raw = directory / "observations.jsonl"
    observations = [json.loads(line) for line in raw.read_text().splitlines()] if raw.exists() else []
    grouped = {}
    for row in observations:
        grouped.setdefault((row["dataset"], row["homology"], row["epsilon_percent"], row["variant"]), []).append(row)
    print(f"{kind}: {len(samples)} inputs x {len(percentages)} epsilons x {len(variants)} configurations; "
          f"up to {repeats} repeats. Tables: {directory}", flush=True)
    order = random.Random(0)
    try:
        with raw.open("a") as stream:
            for sample in samples:
                scc = manifest.parent / sample["scc"]
                base = {key: sample[key] for key in ("dataset", "family", "seed", "source")}
                base.update(n_points=int(sample["n_points"]), target_points=int(sample["target_points"]),
                            homology=int(sample["homology"]), ambient_dimension=int(sample["ambient_dimension"]),
                            scc=str(scc), input_sha256=sample["scc_sha256"],
                            timeout_seconds=timeout, memory_mb=memory_mb, column_storage=storage)
                generated = sample["generation_status"] == "ok"
                if generated:
                    if digest(scc) != sample["scc_sha256"]:
                        raise ValueError(f"Input hash differs from manifest: {scc}")
                    generators, relations, extent = presentation_info(scc)
                    base.update(generators=generators, relations=relations)
                for percentage in percentages:
                    base["epsilon_percent"] = percentage
                    if generated:
                        base["epsilon"] = extent * percentage / 100
                    for repetition in range(repeats):
                        sequence = list(variants)
                        order.shuffle(sequence)
                        if repetition == 0:
                            sequence.sort(key=lambda variant: variant != "baseline")
                        for variant in sequence:
                            key = (sample["dataset"], base["homology"], percentage, variant)
                            previous = grouped.setdefault(key, [])
                            if previous and (len(previous) >= repeats or any(row["status"] != "ok" for row in previous)
                                             or previous[0]["pruning_seconds"] > repeat_below):
                                continue
                            label = f"{sample['dataset']} H{base['homology']} epsilon={percentage:g}% {variant}"
                            print(f"START {label} ({len(previous)+1})", flush=True)
                            result = (execute(pruning, scc, base["epsilon"], variant, timeout, memory_mb, storage)
                                      if generated else dict(status="generation_" + sample["generation_status"], error=sample["error"]))
                            row = dict(base, variant=variant, repetition=len(previous)+1 if generated else 0,
                                       pruning_sha256=binary_hash)
                            row.update(result)
                            stream.write(json.dumps(row, separators=(",", ":")) + "\n")
                            stream.flush()
                            previous.append(row)
                            observations.append(row)
                            write_reports(observations, directory, kind)
                            print(f"DONE {label}: {row['status']} {row.get('pruning_seconds','')} s", flush=True)
    except KeyboardInterrupt:
        print(f"Stopped. Results saved in {directory}; rerun this script to continue.", flush=True)
        return directory
    write_reports(observations, directory, kind)
    print(f"COMPLETE: {directory}; table.csv, table.tex, summary.csv, observations.jsonl, run.json", flush=True)
    return directory

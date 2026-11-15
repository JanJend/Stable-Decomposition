"""Shared CLI execution, provenance, resume, and reporting for paper experiments."""
import csv
from collections import Counter
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

from . import DATA, HERE, ROOT, input_lock, select_samples
from .runtime import digest
from .scaling import SWITCHES, presentation_info
from .paper_tables import VARIANTS, write_reports
from .generator_options import generator_flags


def column_storage(pruning):
    help_text = subprocess.run([str(pruning), "--help"], check=True, capture_output=True,
                               text=True, timeout=10).stdout
    match = re.search(r"^Column storage: (compact|vector)$", help_text, re.MULTILINE)
    if not match:
        raise ValueError(f"Rebuild {pruning}: its --help does not identify column storage")
    return match[1]


def execute(pruning, scc, epsilon, variant, timeout, memory_mb, expected_storage=None, composition_options=None):
    # CLI defaults are fully optimized. Explicitly disable the other options for
    # one-at-a-time comparisons, keeping the unconditional improvements active.
    enabled = set(SWITCHES) - {"full-end-space"} if variant == "all-optimizations" else (
        {"sparse-nilpotent-generators", "complement-orbit-generators"} if composition_options is not None else {variant})
    flags = ["--no-" + name for name in SWITCHES if name != "full-end-space" and name not in enabled]
    if composition_options is not None:
        flags.extend(generator_flags(composition_options))
    if variant == "full-end-space":
        flags.append("--full-end-space")
    flags.append("--hom-algorithm=optimised")
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
        expected = {name.replace("-", "_"): int(name in enabled) for name in SWITCHES}
        algorithm = re.search(r"hom_algorithm=([\w-]+)", options[1]) if options else None
        expected_algorithm = "optimised"
        actual = dict((name, int(value)) for name, value in re.findall(r"(\w+)=(\d)", options[1]) if name in expected) if options else {}
        search = re.search(r"^Composition generator search:([^\n]+)", output, re.MULTILINE)
        actual_search = dict(re.findall(r"([\w-]+)=([^\s]+)", search[1])) if search else {}
        search_matches = composition_options is None or all(
            actual_search.get(name) == str(value) for name, value in composition_options.items())
        if (not search_matches or not total or not iterations or actual != expected or not algorithm or
                algorithm[1] != expected_algorithm or not storage or
                (expected_storage is not None and storage[1] != expected_storage)):
            result.update(status="error", error="Missing profile or unexpected optimization switches/column storage: " + output[-2000:])
        else:
            result.update(pruning_seconds=float(total[1]), i_iterations=int(iterations[1]),
                          k_iterations=int(iterations[2]), column_storage=storage[1],
                          generator_settings=actual_search)
            maps = re.search(r"Initial (?:Hom quotient|full Hom basis): (\d+); I maps: (\d+); K maps: (\d+)", output)
            if maps:
                result.update(initial_maps=int(maps[1]), i_maps=int(maps[2]), k_maps=int(maps[3]))
    else:
        result["error"] = output[-2000:]
    return result


def main(kind, *, points, epsilon_percents, repeats=3, timeout=600, memory_mb=12288,
         repeat_below=5, manifest=None, pruning=None, output_dir=None, families=None, generator_configs=None):
    """Run the fixed experiment; resume matching results and preserve older runs."""
    generator_configs = generator_configs or {}
    if kind != "optimizations" and generator_configs:
        raise ValueError("Generator sweeps belong to the optimization experiment")
    if any(name in VARIANTS or name == "all-optimizations" for name in generator_configs):
        raise ValueError("Generator sweep names must be distinct from the standard variants")
    points = sorted(set(points))
    percentages = sorted(set(epsilon_percents))
    families = sorted(set(families)) if families is not None else None
    if (kind not in ("runtime", "optimizations") or not points or any(n < 6 for n in points)
            or not percentages or any(not 0 < p <= 100 for p in percentages)
            or repeats < 1 or not math.isfinite(repeat_below) or repeat_below <= 0
            or not math.isfinite(timeout) or timeout <= 0
            or not math.isfinite(memory_mb) or memory_mb < 0):
        raise ValueError("Invalid experiment settings: check sizes, epsilons, repetitions and resource limits")
    pruning = ROOT / "build/pruning" if pruning is None else Path(pruning).resolve()
    shared_data = manifest is None
    manifest = DATA / "manifest.csv" if shared_data else Path(manifest).resolve()
    if not pruning.is_file() or not os.access(pruning, os.X_OK):
        raise SystemExit(f"Missing pruning executable: {pruning}. Build it with CMake before running the benchmark.")
    if shared_data and not manifest.is_file():
        subprocess.run([sys.executable, str(HERE / "generate_data.py")], check=True)
    if not manifest.is_file():
        raise SystemExit(f"Prepare the datasets first: python3 {HERE / 'generate_data.py'}")
    with input_lock(manifest.parent):
        storage = column_storage(pruning)
        if storage != "compact":
            raise SystemExit("The paper experiments use compact columns; rebuild with -DGRLINA_USE_COMPACT_COLUMNS=ON")
        binary_hash = digest(pruning)
        with manifest.open(newline="") as stream:
            samples = select_samples(list(csv.DictReader(stream)), points, include_large_immune=kind == "runtime")
        if families is not None:
            samples = [row for row in samples if row["family"] in families]
        samples.sort(key=lambda row: (int(row["n_points"]), row["family"], int(row["homology"]), row["seed"], row["dataset"]))
        if not samples:
            raise ValueError("No manifest inputs match POINT_COUNTS and FAMILIES")
        seed_counts = Counter((row["family"], int(row["n_points"]), int(row["homology"])) for row in samples)
        # Check the entire selection before creating results or timing any case.
        invalid = [str(manifest.parent / row["scc"]) for row in samples
                   if row["generation_status"] == "ok" and
                   (not (manifest.parent / row["scc"]).is_file() or
                    digest(manifest.parent / row["scc"]) != row["scc_sha256"])]
        if invalid:
            raise SystemExit("Inputs missing or different from the manifest:\n" + "\n".join(invalid) +
                             "\nFinish or rerun experiments/generate_data.py, then rerun this benchmark. "
                             "Existing results are preserved.")
        variants = ["all-optimizations"] if kind == "runtime" else [*VARIANTS, *generator_configs]
        config = dict(kind=kind, manifest=str(manifest), manifest_sha256=digest(manifest),
                      pruning=str(pruning), pruning_sha256=binary_hash, column_storage=storage,
                      epsilon_percents=percentages, repeats=repeats, repeat_below=repeat_below,
                      timeout=timeout, memory_mb=memory_mb, points=points, variants=variants,
                      machine=platform.machine(), platform=platform.platform(), python=sys.version,
                      generator_configs=generator_configs)
        if families is not None:
            config["families"] = families
        directory = Path(output_dir).resolve() if output_dir is not None else HERE / "results" / kind
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
                    base["rng_seed"] = sample.get("rng_seed", "")
                    base.update(n_points=int(sample["n_points"]), target_points=int(sample["target_points"]),
                                homology=int(sample["homology"]), ambient_dimension=int(sample["ambient_dimension"]),
                                scc=str(scc), input_sha256=sample["scc_sha256"],
                                timeout_seconds=timeout, memory_mb=memory_mb, column_storage=storage)
                    generated = sample["generation_status"] == "ok"
                    if generated:
                        if digest(scc) != sample["scc_sha256"]:
                            raise SystemExit(f"Input changed during the benchmark: {scc}. "
                                             "Existing results are preserved. Finish generating inputs, then rerun.")
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
                                result = (execute(pruning, scc, base["epsilon"], variant, timeout, memory_mb, storage,
                                                  composition_options=generator_configs.get(variant))
                                          if generated else dict(status="generation_" + sample["generation_status"], error=sample["error"]))
                                row = dict(base, variant=variant, repetition=len(previous)+1 if generated else 0,
                                           pruning_sha256=binary_hash)
                                row.update(result)
                                stream.write(json.dumps(row, separators=(",", ":")) + "\n")
                                stream.flush()
                                previous.append(row)
                                observations.append(row)
                                write_reports(observations, directory, kind, seed_counts=seed_counts, generator_configs=generator_configs)
                                print(f"DONE {label}: {row['status']} {row.get('pruning_seconds','')} s", flush=True)
        except KeyboardInterrupt:
            print(f"Stopped. Results saved in {directory}; rerun this script to continue.", flush=True)
            return directory
        write_reports(observations, directory, kind, seed_counts=seed_counts, generator_configs=generator_configs)
        reports = ("generators.csv, generators.tex, " if generator_configs else "") + "table.csv, table.tex, " + ("speedups.csv, speedups.tex, " if kind == "optimizations" else "")
        print(f"COMPLETE: {directory}; {reports}summary.csv, observations.jsonl, run.json", flush=True)
        return directory

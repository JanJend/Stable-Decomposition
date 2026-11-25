"""Optimization switches, timed worker execution, and speedup calculations."""
import json
import os
import selectors
import statistics
import subprocess
import tempfile
import time

# One-at-a-time comparisons against a quotient-based baseline. Full End is the
# unoptimized control for measuring the quotient's speedup (reverse ratio).
SWITCHES = ["projected-preimages", "amortize-reductions",
            "full-end-space", "nilpotent-reduction", "adaptive-restriction",
            "reuse-target-reductions", "complement-orbit-generators", "single-solve-quotient",
            "update-i-target"]


def presentation_info(path, *, per_axis=False):
    lines = path.read_text().splitlines()
    if lines[:2] != ["scc2020", "2"]:
        raise ValueError(f"Expected a two-parameter SCC presentation: {path}")
    sizes = [int(value) for value in lines[2].split()]
    degrees = [[float(value) for value in line.split(";")[0].split()]
               for line in lines[3:] if line.strip()]
    if len(degrees) != sum(sizes[:-1]) or any(len(degree) != 2 for degree in degrees):
        raise ValueError(f"Unexpected SCC degree table: {path}")
    sides = tuple(max(degree[a] for degree in degrees) - min(degree[a] for degree in degrees)
                  if degrees else 0 for a in range(2))
    return sizes[1], sizes[0], sides if per_axis else max(sides)



def run(worker, path, epsilon, mask, timeout, memory_mb, monitor_memory, label=None):
    start = time.monotonic()
    algorithm_start = None
    deadline = start + timeout
    last_heartbeat = start
    metrics = {}
    status = "error"
    error = "Worker returned no result"
    with tempfile.TemporaryFile(mode="w+b") as log:
        process = subprocess.Popen([str(worker), str(path), str(epsilon), str(mask),
                                    str(memory_mb if monitor_memory else 0)],
                                   stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout, selectors.EVENT_READ)
                while True:
                    now = time.monotonic()
                    if label and now - last_heartbeat >= 30:
                        last_heartbeat = now
                        elapsed = now - (algorithm_start if algorithm_start is not None else start)
                        prefix = os.pread(log.fileno(), 128 * 1024, 0).decode(errors="replace")
                        stage = "iteration" if "dim End_0 = " in prefix else "preparation"
                        print(f"RUNNING {label}: {elapsed:.0f}/{timeout:g}s ({stage})", flush=True)
                    if now >= deadline:
                        status = "timeout" if algorithm_start is not None else "load_timeout"
                        error = f"Exceeded {timeout:g} seconds"
                        if algorithm_start is not None:
                            metrics["seconds_lower_bound"] = now - algorithm_start
                        break
                    events = selector.select(min(1, deadline - now))
                    if not events:
                        continue
                    line = process.stdout.readline()
                    if line.startswith("START"):
                        algorithm_start = time.monotonic()
                        deadline = algorithm_start + timeout
                    elif line.startswith("RESULT "):
                        metrics.update(json.loads(line[7:]))
                        status, error = "ok", ""
                        break
                    elif not line:
                        process.wait()
                        error = f"Worker exited with status {process.returncode}"
                        break
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()
            process.stdout.close()
        log.seek(0)
        prefix = log.read(128 * 1024).decode(errors="replace")
        metrics["basis_completed"] = "dim End_0 = " in prefix
        if status != "ok":
            log.seek(0, os.SEEK_END)
            log.seek(max(0, log.tell() - 1000))
            error += ": " + log.read().decode(errors="replace").strip()
            if process.returncode == 75:
                status = "memory_limit"
    return dict(metrics, status=status, error=error, wall_seconds=time.monotonic() - start)


def summarize(records):
    completed = [record for record in records if record["status"] == "ok"]
    if completed:
        signatures = {(r["i_iterations"], r["k_iterations"], r["generators_out"], r["relations_out"])
                      for r in completed}
        if len(signatures) != 1:
            raise ValueError("Repeated runs disagree on their output/iteration metadata")
        return dict(status="ok", seconds=statistics.median(r["seconds"] for r in completed),
                    samples=len(completed), signature=list(next(iter(signatures))))
    if not records:
        return dict(status="pending", samples=0)
    latest = records[-1]
    result = dict(status=latest["status"], samples=0)
    if result["status"] == "error" and latest.get("error", "").startswith("Worker exited with status -9:"):
        result["status"] = "killed_sigkill"
    elif result["status"] == "error" and "Hom-space is too large (" in latest.get("error", ""):
        result["status"] = "hom_space_capacity_limit"
    if result["status"] == "timeout":
        result["lower_bound_seconds"] = latest["seconds_lower_bound"]
    return result


def speedup(baseline, variant):
    if baseline["status"] == "ok" and variant["status"] == "ok":
        return dict(kind="exact", percent=100 * (baseline["seconds"] / variant["seconds"] - 1))
    if baseline["status"] == "timeout" and variant["status"] == "ok":
        return dict(kind="lower_bound", percent=100 * (baseline["lower_bound_seconds"] / variant["seconds"] - 1))
    if baseline["status"] == "ok" and variant["status"] == "timeout":
        return dict(kind="upper_bound", percent=100 * (baseline["seconds"] / variant["lower_bound_seconds"] - 1))
    return dict(kind="unavailable", reason=f"baseline:{baseline['status']};switch:{variant['status']}")



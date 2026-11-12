"""Worker subprocess handling shared by runtime experiments and their checks."""
import hashlib
import json
import os
import subprocess
import tempfile
import time


def digest(path):
    with path.open("rb") as stream:
        result = hashlib.sha256()
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
        return result.hexdigest()


def run_worker(worker, arguments, timeout):
    start = time.perf_counter()
    metrics = {}
    status, error = "ok", ""
    # Progress can be large. Keep it in an automatically removed temporary file,
    # independently of the tiny metrics stream; include its tail only on failure.
    with tempfile.TemporaryFile(mode="w+b") as log:
        process = subprocess.Popen([str(worker), *map(str, arguments)],
                                   stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            output, _ = process.communicate(timeout=timeout)
            if process.returncode:
                status = "error"
                error = f"Worker exit status {process.returncode}"
        except subprocess.TimeoutExpired:
            process.kill()
            output, _ = process.communicate()
            status, error = "timeout", f"Exceeded {timeout:g}s"
        except BaseException:
            process.kill()
            process.wait()
            raise
        if status != "ok":
            log.seek(0, os.SEEK_END)
            log.seek(max(0, log.tell() - 1000))
            error += ": " + log.read().decode(errors="replace").strip()
    for line in output.splitlines():
        if line.startswith("RESULT "):
            metrics.update(json.loads(line[7:]))
    if status == "ok" and not metrics:
        status, error = "error", "Worker produced no metrics"
    return metrics, status, time.perf_counter() - start, error



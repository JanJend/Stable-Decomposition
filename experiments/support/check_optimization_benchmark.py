"""Small runnable check of timing, censoring, and percentage calculations."""

from . import HERE
from .scaling import presentation_info, run, speedup, summarize


def check(worker):
    path = HERE / "test_presentations/torus_gaussian/torus_gaussian_n20_seed42_h1.scc"
    generators, relations, extent = presentation_info(path)
    assert (generators, relations) == (12, 12)
    baseline = run(worker, path, extent * 0.01, 0, 10, 4096, False)
    optimized = run(worker, path, extent * 0.01, 1023, 10, 4096, False)
    assert baseline["status"] == optimized["status"] == "ok"
    assert baseline["seconds"] >= baseline["profile_seconds"] > 0
    assert all(baseline[k] == optimized[k] for k in
               ("i_iterations", "k_iterations", "generators_out", "relations_out"))
    larger = HERE / "test_presentations/torus_gaussian/torus_gaussian_n128_seed42_h1.scc"
    timeout = run(worker, larger, 0.008837960620266386, 0, 0.01, 4096, False)
    assert timeout["status"] == "timeout" and timeout["seconds_lower_bound"] >= 0.01
    memory = run(worker, larger, 0.008837960620266386, 0, 10, 1, True)
    assert memory["status"] == "memory_limit"
    assert speedup(dict(status="ok", seconds=2), dict(status="ok", seconds=1))["percent"] == 100
    assert speedup(dict(status="ok", seconds=1), dict(status="ok", seconds=2))["percent"] == -50
    assert speedup(dict(status="timeout", lower_bound_seconds=10), dict(status="ok", seconds=1)) == dict(kind="lower_bound", percent=900)
    assert speedup(dict(status="ok", seconds=1), dict(status="timeout", lower_bound_seconds=10)) == dict(kind="upper_bound", percent=-90)
    assert speedup(dict(status="timeout", lower_bound_seconds=10), dict(status="timeout", lower_bound_seconds=10))["kind"] == "unavailable"
    assert summarize([dict(status="memory_limit"), dict(status="timeout", seconds_lower_bound=600)]) == dict(status="timeout", samples=0, lower_bound_seconds=600)
    assert summarize([dict(status="error", error="Worker exited with status -9: interrupted")])["status"] == "killed_sigkill"
    assert summarize([dict(status="error", error="Hom-space is too large (113843) to convert")])["status"] == "hom_space_capacity_limit"
    print("PASS: worker timing, output metadata, timeout, memory guard, and speedup arithmetic")



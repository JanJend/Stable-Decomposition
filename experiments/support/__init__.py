"""Internal helpers for the commands in experiments/."""
import math
from contextlib import contextmanager
import fcntl
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
DATA = HERE / "data"
RUNTIME_SIZES = [20 * 2 ** power for power in range(8)]
OPTIMIZATION_SIZES = RUNTIME_SIZES[:6]
LARGE_IMMUNE = ("LargeHypoxicRegion", "LargeHypoxicRegion2")


@contextmanager
def input_lock(directory, *, writing=False):
    """Keep generated inputs fixed while benchmark processes read them."""
    with (directory / ".inputs.lock").open("a") as lock:
        try:
            fcntl.flock(lock, (fcntl.LOCK_EX if writing else fcntl.LOCK_SH) | fcntl.LOCK_NB)
        except BlockingIOError:
            operation = "regenerate inputs" if writing else "benchmark inputs"
            raise SystemExit(f"Cannot {operation}: another process is using {directory}. "
                             "Wait for it to finish, then rerun.") from None
        yield


def select_samples(rows, sizes, include_large_immune=False):
    """Share inputs, selecting complete observed clouds nearest each target size."""
    selected = [row for row in rows if row["family"] != "immuno" and int(row["target_points"]) in sizes]
    observed = {row["dataset"]: int(row["n_points"]) for row in rows if row["family"] == "immuno"}
    targets = {}
    for target in sorted(set(sizes)):
        available = observed.keys() - targets.keys()
        if available:
            name = min(available, key=lambda name: (abs(math.log(observed[name] / target)), name))
            targets[name] = target
    if include_large_immune:
        targets.update({name: observed[name] for name in LARGE_IMMUNE if name in observed})
    selected.extend(dict(row, target_points=targets[row["dataset"]]) for row in rows
                    if row["family"] == "immuno" and row["dataset"] in targets)
    return selected


#!/usr/bin/env python3
"""Generate synthetic/immune-cell inputs and function-Delaunay presentations."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from support import DATA, ROOT, RUNTIME_SIZES, OPTIMIZATION_SIZES, LARGE_IMMUNE, input_lock
SIZES = sorted(set(RUNTIME_SIZES + OPTIMIZATION_SIZES))
FAMILIES = ["torus_gauss", "torus_uni", "gauss_clust", "S2", "immuno"]
HOMOLOGY = {"torus_gauss": [1], "torus_uni": [1], "gauss_clust": [0], "S2": [2], "immuno": [1]}
# Fixed widths for the other families; tori use leave-one-out selection below.
# Synthetic widths are in original units; immune widths are in scaled coordinates.
KDE_BANDWIDTHS = {"gauss_clust": 0.2, "S2": 0.2, "immuno": 0.05}
TORUS_BANDWIDTHS = np.geomspace(0.1, 2.0, 25)
IMMUNE_URL = "https://github.com/MultiparameterTDAHistology/SpatialPatterningOfImmuneCells.git"


def torus(rng, n):
    # Same area-weighted rejection sampling as mpm_generation/noisy_torus.py.
    major, minor = 2.0, 1.0
    u = rng.uniform(0, 2 * np.pi, n)
    angles = []
    while len(angles) < n:
        v = rng.uniform(0, 2 * np.pi, 2 * n)
        keep = rng.random(len(v)) < (major + minor * np.cos(v)) / (major + minor)
        angles.extend(v[keep].tolist())
    v = np.asarray(angles[:n])
    return np.column_stack(((major + minor * np.cos(v)) * np.cos(u),
                            (major + minor * np.cos(v)) * np.sin(u), minor * np.sin(v)))


def mixed_noise(rng, points, sigma=0.27, outlier_fraction=0.03, cauchy_scale=0.6):
    """Gaussian core and clipped Cauchy displacements, as in noisy_torus.add_mixed_noise."""
    points = points + rng.normal(0, sigma, points.shape)
    count = max(1, int(len(points) * outlier_fraction))
    indices = rng.choice(len(points), count, replace=False)
    directions = rng.normal(size=(count, points.shape[1]))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    distances = np.clip(rng.standard_cauchy(count) * cauchy_scale, -3 * cauchy_scale, 3 * cauchy_scale)
    points[indices] += directions * distances[:, None]
    return points


def cluster_centers(rng, scales):
    """Randomize centers while keeping their three-sigma core balls disjoint."""
    anchors = np.asarray([[-2, 0, -0.75], [1.5, -1, 0.5], [1.5, 1.5, 1.5]])
    while True:
        centers = anchors + rng.normal(0, 0.75, anchors.shape)
        if all(np.linalg.norm(centers[i] - centers[j]) >= 3 * (scales[i] + scales[j])
               for i, j in ((0, 1), (0, 2), (1, 2))):
            return centers


def synthetic(family, n, seed):
    # Stable family seeds: adding another family does not change existing samples.
    rng = np.random.default_rng(np.random.SeedSequence([seed, FAMILIES.index(family), n]))
    if family.startswith("torus"):
        points = mixed_noise(rng, torus(rng, n))
        if family == "torus_uni":
            count = max(1, round(0.1 * n))
            points[rng.choice(n, count, replace=False)] = rng.uniform([-3, -3, -1], [3, 3, 1], (count, 3))
        return points
    if family not in ("gauss_clust", "S2"):
        raise ValueError(f"Unknown synthetic family: {family}")
    # Clusters and S2 share the exact Gaussian mixture and outlier model.
    labels = rng.choice(3, n, p=[0.4, 0.35, 0.25])
    scales = np.asarray([0.6, 0.8, 0.5])
    if family == "gauss_clust":
        centers = cluster_centers(rng, scales)
        points = centers[labels]
    else:
        points = rng.normal(size=(n, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
    points = points + rng.normal(size=(n, 3)) * scales[labels, None]
    return mixed_noise(rng, points, sigma=0, outlier_fraction=0.05, cauchy_scale=2.0)


def immune_samples(repo, sizes):
    folder = repo / "1.5mmRegions"
    if not folder.is_dir():
        raise SystemExit(f"Clone {IMMUNE_URL} to {repo}, or pass --immune-repo PATH.")
    def read_cloud(path):
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if not {"x", "y"}.issubset(reader.fieldnames or []):
                raise ValueError(f"Expected x,y columns in {path}")
            # The large CSVs contain CD8, CD68 and FoxP3 together. Keep all rows.
            points = np.asarray([[float(row["x"]), float(row["y"])] for row in reader])
        if points.ndim != 2 or points.shape[0] < 5:
            return None
        if not np.isfinite(points).all():
            raise ValueError(f"Nonfinite coordinate in {path}")
        return np.unique(points, axis=0)

    candidates = []
    # Scan every cell-type subfolder; use complete observed clouds, no subsampling.
    for path in sorted(folder.rglob("*.csv")):
        if any(part.startswith(".") for part in path.relative_to(folder).parts):
            continue
        points = read_cloud(path)
        if points is not None:
            candidates.append((path, points))
    if not candidates:
        raise ValueError(f"No usable point clouds in {folder}")
    print(f"Scanned {len(candidates)} immune-cell point clouds.", flush=True)
    used = set()
    selected = []
    for target in sorted(set(sizes)):
        if len(used) == len(candidates):
            raise ValueError(f"Not enough distinct immune-cell clouds for {len(set(sizes))} sizes")
        path, points = min((item for item in candidates if item[0] not in used),
                           key=lambda item: (abs(math.log(len(item[1]) / target)), str(item[0])))
        used.add(path)
        selected.append((target, path, points))
    for name in LARGE_IMMUNE:
        path = repo / f"{name}.csv"
        points = read_cloud(path)
        if points is None:
            raise ValueError(f"No usable immune cells in {path}")
        selected.append((len(points), path, points))
    return selected


def select_torus_bandwidth(points):
    """Maximize mean leave-one-out Gaussian log density in original spatial units."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError("Bandwidth selection needs at least two finite points in R3")
    points = points - points.mean(axis=0)
    n, dimension = points.shape
    squared = np.sum(points * points, axis=1)
    scores = np.zeros(len(TORUS_BANDWIDTHS))
    # ponytail: exact O(25*n^2) search, O(128*n) memory; use a tree if samples grow much larger.
    for start in range(0, n, 128):
        stop = min(start + 128, n)
        distances = np.maximum(squared[start:stop, None] + squared[None, :]
                               - 2 * points[start:stop] @ points.T, 0)
        distances[np.arange(stop - start), np.arange(start, stop)] = np.inf
        for index, bandwidth in enumerate(TORUS_BANDWIDTHS):
            log_kernel = -distances / (2 * bandwidth ** 2)
            maximum = log_kernel.max(axis=1)
            # Stable log-sum-exp keeps isolated points meaningful at small widths.
            scores[index] += np.sum(maximum + np.log(np.exp(log_kernel - maximum[:, None]).sum(axis=1)))
    scores = (scores / n - np.log(n - 1) - dimension / 2 * np.log(2 * np.pi)
              - dimension * np.log(TORUS_BANDWIDTHS))
    best = int(np.argmax(scores))
    return float(TORUS_BANDWIDTHS[best]), float(scores[best])


def with_codensity(points, bandwidth, *, scaled_bandwidth=False):
    # A common spatial scale makes the distance and [0,1] density axes comparable.
    points = np.asarray(points, dtype=float)
    if (points.ndim != 2 or len(points) < 2 or not np.isfinite(points).all()
            or not math.isfinite(bandwidth) or bandwidth <= 0):
        raise ValueError("Expected finite points and a finite positive KDE bandwidth")
    lower = points.min(axis=0)
    scale = float(np.ptp(points, axis=0).max())
    if scale <= 0:
        raise ValueError("Point cloud has zero extent")
    points = (points - lower) / scale
    n = len(points)
    bandwidth = float(bandwidth if scaled_bandwidth else bandwidth / scale)
    squared = np.sum(points * points, axis=1)
    density = np.empty(n)
    # ponytail: exact O(n^2) Gaussian KDE, O(128*n) memory; use a tree if larger clouds demand it.
    for start in range(0, n, 128):
        block = points[start:start + 128]
        distances = np.maximum(squared[start:start + 128, None] + squared[None, :] - 2 * block @ points.T, 0)
        density[start:start + 128] = np.exp(-distances / (2 * bandwidth ** 2)).mean(axis=1)
    extent = float(np.ptp(density))
    density = (density - density.min()) / extent if extent else np.zeros(n)
    return np.column_stack((points, -density)), scale, lower, bandwidth


def write_manifest(path, rows):
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=FAMILIES)
    parser.add_argument("--sizes", type=int, nargs="+", default=SIZES)
    parser.add_argument("--seeds", "--seed", dest="seeds", type=int, nargs="+", default=[42, 43, 44],
                        help="RNG seeds; filenames use consecutive sample numbers starting at 1")
    parser.add_argument("--homology", type=int, nargs="+", help="Override defaults: H0 clusters, H2 S2, H1 tori and immune")
    parser.add_argument("--bandwidth", type=float, help="Override automatic torus bandwidth or fixed family width (original units; scaled units for immune cells)")
    parser.add_argument("--immune-repo", type=Path, default=ROOT.parent / "SpatialPatterningOfImmuneCells")
    parser.add_argument("--function-delaunay", type=Path, default=ROOT.parent / "function_delaunay/build/main")
    parser.add_argument("--output-dir", type=Path, default=DATA)
    parser.add_argument("--delaunay-timeout", type=float, default=1200)
    parser.add_argument("--force", action="store_true", help="Recompute even when point data and homology match")
    args = parser.parse_args()
    if (any(n < 6 for n in args.sizes) or any(seed < 0 for seed in args.seeds)
            or not (0 < args.delaunay_timeout < float("inf"))):
        parser.error("Sizes must be >=6, seed nonnegative, and timeout finite and positive")
    if args.homology and any(h < 0 for h in args.homology):
        parser.error("Homology degrees must be nonnegative")
    if args.bandwidth is not None and (not math.isfinite(args.bandwidth) or args.bandwidth <= 0):
        parser.error("Bandwidth must be finite and positive")
    executable = args.function_delaunay.resolve()
    if not executable.is_file():
        parser.error(f"Build function_delaunay first: missing {executable}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with input_lock(output, writing=True):
        manifest = output / "manifest.csv"
        previous = {}
        if manifest.exists():
            with manifest.open(newline="") as stream:
                previous = {row["scc"]: row for row in csv.DictReader(stream)}
        rows = list(previous.values())
        # Keep sample numbers stable when rerunning only a subset of RNG seeds.
        seed_numbers = {int(row["rng_seed"]): int(row["seed"]) for row in rows if row["seed"] != ""}
        for rng_seed in dict.fromkeys(args.seeds):
            if rng_seed not in seed_numbers:
                seed_numbers[rng_seed] = max(seed_numbers.values(), default=0) + 1
        samples = []
        for family in args.families:
            if family == "immuno":
                for target, path, points in immune_samples(args.immune_repo.resolve(), args.sizes):
                    samples.append((family, target, path.stem, str(path), points, None, None))
            else:
                for rng_seed in dict.fromkeys(args.seeds):
                    for n in sorted(set(args.sizes)):
                        seed = seed_numbers[rng_seed]
                        samples.append((family, n, f"{family}_{n}_{seed}", "synthetic",
                                        synthetic(family, n, rng_seed), seed, rng_seed))
        for family, target, name, source, points, seed, rng_seed in samples:
            loo_score = ""
            if args.bandwidth is not None:
                bandwidth, bandwidth_method = args.bandwidth, "override"
            elif family.startswith("torus"):
                bandwidth, loo_score = select_torus_bandwidth(points)
                bandwidth_method = "loo_cv"
                print(f"{name}: selected KDE bandwidth {bandwidth:.6g} (original units)", flush=True)
            else:
                bandwidth, bandwidth_method = KDE_BANDWIDTHS[family], "fixed"
            data, scale, lower, bandwidth = with_codensity(points, bandwidth, scaled_bandwidth=family == "immuno")
            pointset = output / "point-sets" / family / f"{name}.txt"
            pointset.parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(pointset, data, fmt="%.17g")
            digest = hashlib.sha256(pointset.read_bytes()).hexdigest()
            # An override can change a cloud shared by several homology degrees.
            # Drop sibling entries computed from older codensities rather than keeping invalid hashes.
            relative_pointset = pointset.relative_to(output).as_posix()
            rows = [row for row in rows if row["pointset"] != relative_pointset or row["pointset_sha256"] == digest]
            degrees = args.homology or HOMOLOGY[family]
            for degree in degrees:
                if degree >= points.shape[1]:
                    parser.error(f"Homology degree {degree} must be smaller than ambient dimension {points.shape[1]}")
                scc = output / "presentations" / family / f"{name}_H{degree}.scc"
                scc.parent.mkdir(parents=True, exist_ok=True)
                relative = scc.relative_to(output).as_posix()
                row = dict(dataset=name, family=family, target_points=target, n_points=len(points),
                           ambient_dimension=points.shape[1], homology=degree, source=source,
                           seed=seed if source == "synthetic" else "", rng_seed=rng_seed if source == "synthetic" else "",
                           spatial_scale=scale,
                           spatial_offset=json.dumps(lower.tolist()), kde_bandwidth=bandwidth,
                           kde_bandwidth_method=bandwidth_method, kde_bandwidth_loo_score=loo_score,
                           pointset=relative_pointset, pointset_sha256=digest,
                           scc=relative, scc_sha256="", delaunay_seconds="", generation_status="ok", error="")
                old = previous.get(relative, {})
                if (not args.force and scc.exists() and old.get("pointset_sha256") == digest
                        and old.get("generation_status") == "ok"
                        and hashlib.sha256(scc.read_bytes()).hexdigest() == old.get("scc_sha256")):
                    row["scc_sha256"] = old["scc_sha256"]
                    row["delaunay_seconds"] = old["delaunay_seconds"]
                else:
                    start = time.perf_counter()
                    try:
                        result = subprocess.run([str(executable), "--minpres", str(degree), "--no-delaunay-compare",
                                                 str(pointset), str(scc)], capture_output=True, text=True,
                                                timeout=args.delaunay_timeout)
                        if result.returncode or not scc.is_file() or not scc.read_bytes().startswith(b"scc2020"):
                            raise RuntimeError((result.stderr + result.stdout)[-1500:])
                        row["scc_sha256"] = hashlib.sha256(scc.read_bytes()).hexdigest()
                    except subprocess.TimeoutExpired:
                        row["generation_status"] = "timeout"
                        row["error"] = f"Function-Delaunay exceeded {args.delaunay_timeout:g}s"
                    except (OSError, RuntimeError) as error:
                        row["generation_status"] = "error"
                        row["error"] = str(error)
                    row["delaunay_seconds"] = time.perf_counter() - start
                    if row["generation_status"] != "ok":
                        scc.unlink(missing_ok=True)
                rows = [item for item in rows if item["scc"] != relative] + [row]
                write_manifest(manifest, rows)
                print(f"{family}: n={len(points)} H{degree} {row['generation_status']} -> {relative}", flush=True)
        print(f"Manifest: {manifest}")
        return int(any(row["generation_status"] != "ok" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())

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

from support import HERE, ROOT
SIZES = [20, 32, 64, 128, 256, 512, 1024, 2048]
FAMILIES = ["torus_gaussian", "torus_uniform", "gaussian_clusters", "sphere4d", "immuno"]
IMMUNE_URL = "https://github.com/MultiparameterTDAHistology/SpatialPatterningOfImmuneCells.git"


def torus(rng, n):
    # Same area-weighted rejection sampling as mpm_generation/noisy_torus.py.
    major, minor = 3.0, 1.0
    u = rng.uniform(0, 2 * np.pi, n)
    angles = []
    while len(angles) < n:
        v = rng.uniform(0, 2 * np.pi, 2 * n)
        keep = rng.random(len(v)) < (major + minor * np.cos(v)) / (major + minor)
        angles.extend(v[keep].tolist())
    v = np.asarray(angles[:n])
    return np.column_stack(((major + minor * np.cos(v)) * np.cos(u),
                            (major + minor * np.cos(v)) * np.sin(u), minor * np.sin(v)))


def synthetic(family, n, seed):
    # Stable family seeds: adding another family does not change existing samples.
    rng = np.random.default_rng(np.random.SeedSequence([seed, FAMILIES.index(family), n]))
    if family.startswith("torus"):
        points = torus(rng, n) + rng.normal(0, 0.1, (n, 3))
        if family == "torus_uniform":
            count = max(1, round(0.1 * n))
            points[rng.choice(n, count, replace=False)] = rng.uniform([-4, -4, -1], [4, 4, 1], (count, 3))
        return points
    if family == "gaussian_clusters":
        labels = rng.choice(3, n, p=[0.4, 0.35, 0.25])
        centers = np.asarray([[-2, 0], [1.5, -1], [1.5, 1.5]])
        scales = np.asarray([0.3, 0.4, 0.25])
        return centers[labels] + rng.normal(size=(n, 2)) * scales[labels, None]
    points = rng.normal(size=(n, 4))
    points /= np.linalg.norm(points, axis=1)[:, None]
    return points + rng.normal(0, 0.05, points.shape)


def immune_samples(repo, sizes):
    folder = repo / "1.5mmRegions"
    if not folder.is_dir():
        raise SystemExit(f"Clone {IMMUNE_URL} to {repo}, or pass --immune-repo PATH.")
    candidates = []
    # Scan every cell-type subfolder; use complete observed clouds, no subsampling.
    for path in sorted(folder.rglob("*.csv")):
        if any(part.startswith(".") for part in path.relative_to(folder).parts):
            continue
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["x", "y"]:
                raise ValueError(f"Expected x,y columns in {path}")
            points = np.asarray([[float(row["x"]), float(row["y"])] for row in reader])
        if points.ndim != 2 or points.shape[0] < 5:
            continue
        if not np.isfinite(points).all():
            raise ValueError(f"Nonfinite coordinate in {path}")
        points = np.unique(points, axis=0)
        candidates.append((path, points))
    if not candidates:
        raise ValueError(f"No usable point clouds in {folder}")
    print(f"Scanned {len(candidates)} immune-cell point clouds.", flush=True)
    used = set()
    selected = []
    for target in sizes:
        path, points = min((item for item in candidates if item[0] not in used),
                           key=lambda item: (abs(math.log(len(item[1]) / target)), str(item[0])))
        used.add(path)
        selected.append((target, path, points))
    return selected


def with_codensity(points):
    # A common spatial scale makes the distance and [0,1] density axes comparable.
    lower = points.min(axis=0)
    scale = float(np.ptp(points, axis=0).max())
    if scale <= 0:
        raise ValueError("Point cloud has zero extent")
    points = (points - lower) / scale
    n, dimension = points.shape
    bandwidth = float(np.std(points, axis=0).mean() * n ** (-1 / (dimension + 4)))
    if bandwidth <= 0:
        raise ValueError("KDE bandwidth must be positive")
    squared = np.sum(points * points, axis=1)
    density = np.empty(n)
    # ponytail: exact quadratic Gaussian KDE; block the distance matrix for 2048 points.
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
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=FAMILIES)
    parser.add_argument("--sizes", type=int, nargs="+", default=SIZES)
    parser.add_argument("--seeds", "--seed", dest="seeds", type=int, nargs="+", default=[42, 43, 44])
    parser.add_argument("--homology", type=int, nargs="+", help="Override default: H0 clusters, H1 otherwise")
    parser.add_argument("--immune-repo", type=Path, default=ROOT.parent / "SpatialPatterningOfImmuneCells")
    parser.add_argument("--function-delaunay", type=Path, default=ROOT.parent / "function_delaunay/build/main")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--delaunay-timeout", type=float, default=1200)
    parser.add_argument("--force", action="store_true", help="Recompute even when point data and homology match")
    args = parser.parse_args()
    if (any(n < 6 for n in args.sizes) or any(seed < 0 for seed in args.seeds)
            or not (0 < args.delaunay_timeout < float("inf"))):
        parser.error("Sizes must be >=6, seed nonnegative, and timeout finite and positive")
    if args.homology and any(h < 0 for h in args.homology):
        parser.error("Homology degrees must be nonnegative")
    executable = args.function_delaunay.resolve()
    if not executable.is_file():
        parser.error(f"Build function_delaunay first: missing {executable}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "manifest.csv"
    previous = {}
    if manifest.exists():
        with manifest.open(newline="") as stream:
            previous = {row["scc"]: row for row in csv.DictReader(stream)}
    rows = list(previous.values())
    samples = []
    for family in args.families:
        if family == "immuno":
            for target, path, points in immune_samples(args.immune_repo.resolve(), args.sizes):
                samples.append((family, target, path.stem, str(path), points, None))
        else:
            for seed in dict.fromkeys(args.seeds):
                for n in args.sizes:
                    samples.append((family, n, f"{family}_n{n}_seed{seed}", "synthetic", synthetic(family, n, seed), seed))
    for family, target, name, source, points, seed in samples:
        data, scale, lower, bandwidth = with_codensity(points)
        pointset = output / "point-sets" / family / f"{name}.txt"
        pointset.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(pointset, data, fmt="%.17g")
        digest = hashlib.sha256(pointset.read_bytes()).hexdigest()
        degrees = args.homology or ([0] if family == "gaussian_clusters" else [1])
        for degree in degrees:
            if degree >= points.shape[1]:
                parser.error(f"Homology degree {degree} must be smaller than ambient dimension {points.shape[1]}")
            scc = output / "test_presentations" / family / f"{name}_h{degree}.scc"
            scc.parent.mkdir(parents=True, exist_ok=True)
            relative = scc.relative_to(output).as_posix()
            row = dict(dataset=name, family=family, target_points=target, n_points=len(points),
                       ambient_dimension=points.shape[1], homology=degree, source=source,
                       seed=seed if source == "synthetic" else "", spatial_scale=scale,
                       spatial_offset=json.dumps(lower.tolist()), kde_bandwidth=bandwidth,
                       pointset=pointset.relative_to(output).as_posix(), pointset_sha256=digest,
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

#!/usr/bin/env python3
"""Plot shared input clouds and Hilbert functions in data/pictures/index.html."""
import argparse
import csv
import hashlib
import inspect
import json
from html import escape
from itertools import combinations
from pathlib import Path
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from support import DATA, ROOT, RUNTIME_SIZES, OPTIMIZATION_SIZES, select_samples


def plot_pointset(path, row, output):
    data = np.loadtxt(path, ndmin=2)
    n, dimension = int(row["n_points"]), int(row["ambient_dimension"])
    if (dimension not in (2, 3, 4) or data.shape != (n, dimension + 1)
            or not np.isfinite(data).all() or np.any(data[:, -1] < -1 - 1e-12)
            or np.any(data[:, -1] > 1e-12)):
        raise ValueError(f"Expected {n} points in R{dimension} with a final codensity column: {path}")
    coordinates, density = data[:, :-1], -data[:, -1]
    figure = plt.figure(figsize=(12, 8) if dimension == 4 else (8, 7), layout="constrained")
    try:
        figure.suptitle(row.get("dataset", Path(row["scc"]).stem) + f" — {n} points in R{dimension}", fontsize=15)
        options = dict(c=density, cmap="viridis", vmin=0, vmax=1, s=14, linewidths=0)
        axes = []
        if dimension == 3:
            axis = figure.add_subplot(projection="3d")
            cloud = axis.scatter(*coordinates.T, depthshade=False, **options)
            for setter, label in ((axis.set_xlabel, "x1"), (axis.set_ylabel, "x2"), (axis.set_zlabel, "x3")):
                setter(label)
            axis.set_box_aspect(np.maximum(np.ptp(coordinates, axis=0), 1e-12), zoom=0.85)
            axes.append(axis)
        else:
            pairs = list(combinations(range(dimension), 2)) if dimension == 4 else [(0, 1)]
            for number, (i, j) in enumerate(pairs, 1):
                axis = figure.add_subplot(2, 3, number) if dimension == 4 else figure.add_subplot()
                cloud = axis.scatter(coordinates[:, i], coordinates[:, j], **options)
                axis.set(xlabel=f"x{i+1}", ylabel=f"x{j+1}", aspect="equal")
                axis.grid(alpha=0.2)
                axes.append(axis)
        figure.colorbar(cloud, ax=axes, shrink=0.8, label="Normalised KDE (higher = denser)")
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=180, facecolor="white", bbox_inches="tight")
    finally:
        plt.close(figure)


def write_gallery(output, records, experiment=None):
    fields = list(dict.fromkeys(key for row in records for key in row))
    with (output / "index.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    title = {"runtime": "Runtime experiment", "optimizations": "Optimization comparison"}.get(experiment, "Shared experiment data")
    page = [f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}: input pictures</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;color:#202a36;background:#f5f7f9;max-width:1600px;margin:auto;padding:24px}}
h1{{margin-bottom:8px}}h2{{margin-top:40px}}h3{{margin:0;overflow-wrap:anywhere}}
p{{max-width:1200px}}a{{color:#1456a0}}nav a{{display:inline-block;margin:0 24px 8px 0}}
article{{background:white;border:1px solid #d8dee6;border-radius:8px;padding:16px;margin-bottom:24px}}
.plots{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,440px),1fr));gap:16px}}
figure{{margin:12px 0 0}}img{{display:block;width:100%;height:auto}}figcaption,.meta{{color:#435164}}
</style>
<h1>{title}: input pictures</h1>
<p>{len(records)} input modules, before pruning. The same inputs are used at every epsilon and optimization setting.</p>
<p><strong>Point sets:</strong> spatial coordinates have been translated and isotropically scaled so the largest
coordinate range is one. Colour shows normalised Gaussian KDE, with yellow indicating higher density.
Immune cells are shown in 2D; clusters, tori and S2 are shown in 3D.
All points are included. Sample numbers identify RNG seeds; the original seed is recorded in each caption.</p>
<p><strong>Hilbert functions:</strong> dim H<sub>k</sub>(r,t) over F₂ for the sublevel Delaunay–Čech bifiltration.
X is minimum enclosing ball radius; Y is negative normalised KDE, with higher density entering first near −1.
White means dimension zero; positive dimensions use logarithmic blue colours. Each plot has its own radius extent
and colour maximum. Values are exact at 1000 × 700 sampled coordinates, using the complete projective resolution;
very narrow features can fall between samples.</p>
<p>Torus bandwidths are selected for each cloud by leave-one-out Gaussian likelihood
over 25 logarithmically spaced widths from 0.1 to 2.0 in original spatial units.
Clusters and S2 use 0.2 in original units; immune cells use 0.05 in scaled coordinates.
The final KDE includes self-contributions; bandwidth selection excludes them.
Captions give the actual h and hL in original spatial units (including any command-line override),
where L is the spatial scale. Click pictures for full resolution.
<a href="index.csv">Metadata and input hashes</a>.</p><nav>''']
    families = list(dict.fromkeys(row["family"] for row in records))
    page.extend(f'<a href="#{escape(family)}">{escape(family)}</a>' for family in families)
    page.append('</nav>')
    for family in families:
        page.append(f'<h2 id="{escape(family)}">{escape(family)}</h2>')
        selected = sorted((row for row in records if row["family"] == family),
                          key=lambda row: (int(row["n_points"]), row["seed"], row["dataset"]))
        for row in selected:
            name = escape(Path(row["scc"]).stem)
            sample = f" · sample {row['seed']} (RNG {row['rng_seed']})" if row["seed"] else ''
            method = {"loo_cv": "leave-one-out selection", "override": "manual override"}.get(row.get("kde_bandwidth_method"), "fixed bandwidth")
            page.append(f'''<article><h3>{name}</h3><div class="meta">n = {row['n_points']} · H{row['homology']}{sample}<br>
KDE h = {float(row['kde_bandwidth']):.8g} (scaled); hL = {float(row['kde_bandwidth'])*float(row['spatial_scale']):.8g} (original) · {method}
</div><div class="plots">''')
            if row.get("generation_status") != "ok":
                page.append(f'<p>Presentation generation: {escape(row.get("generation_status", "unknown"))}. {escape(row.get("error", ""))}</p>')
            for field, label in (("pointset_image", "Point set"), ("image", "Hilbert function")):
                if row.get(field) and (output / row[field]).is_file():
                    image = escape(row[field], quote=True)
                    page.append(f'<figure><a href="{image}"><img src="{image}" loading="lazy" alt="{label}: {name}"></a><figcaption>{label}</figcaption></figure>')
            page.append('</div></article>')
    page.append('</html>')
    (output / "index.html").write_text('\n'.join(page))


def plot_experiment(experiment=None, bundle=DATA, renderer=ROOT / "build/experiment_hilbert"):
    output = bundle / "pictures"
    if not renderer.is_file():
        raise SystemExit("Build the Hilbert renderer first: cmake --build build --target experiment_hilbert -j 2")
    renderer_hash = hashlib.sha256(renderer.read_bytes()).hexdigest()
    point_renderer_hash = hashlib.sha256(inspect.getsource(plot_pointset).encode()).hexdigest()
    with (bundle / "manifest.csv").open(newline="") as stream:
        all_inputs = list(csv.DictReader(stream))
    inputs = all_inputs
    if experiment:
        inputs = select_samples(inputs, RUNTIME_SIZES if experiment == "runtime" else OPTIMIZATION_SIZES,
                                include_large_immune=experiment == "runtime")
    output.mkdir(parents=True, exist_ok=True)
    previous = {}
    if (output / "index.csv").is_file():
        with (output / "index.csv").open(newline="") as stream:
            previous = {row["scc"]: row for row in csv.DictReader(stream)}
    records = []
    plotted = set()
    try:
        for number, row in enumerate(inputs, 1):
            source = bundle / row["pointset"]
            if hashlib.sha256(source.read_bytes()).hexdigest() != row["pointset_sha256"]:
                raise ValueError(f"Point-set hash differs from manifest: {source}")
            old = previous.get(row["scc"], {})
            stem = Path(row["scc"]).stem
            picture = Path(row["family"]) / f"{row['dataset']}_points.png"
            if source not in plotted and (not (output / picture).is_file() or old.get("pointset_sha256") != row["pointset_sha256"]
                                          or old.get("point_renderer_sha256") != point_renderer_hash):
                plot_pointset(source, row, output / picture)
            plotted.add(source)
            record = dict(row, pointset_image=picture.as_posix(), image="", renderer_sha256=renderer_hash,
                          point_renderer_sha256=point_renderer_hash)
            if row["generation_status"] == "ok":
                presentation = bundle / row["scc"]
                if hashlib.sha256(presentation.read_bytes()).hexdigest() != row["scc_sha256"]:
                    raise ValueError(f"Presentation hash differs from manifest: {presentation}")
                image_path = Path(row["family"]) / f"{stem}.png"
                if (not (output / image_path).is_file() or old.get("scc_sha256") != row["scc_sha256"]
                        or old.get("renderer_sha256") != renderer_hash):
                    result = subprocess.run([str(renderer), str(presentation), str(output / image_path)],
                                            capture_output=True, text=True, timeout=180, check=True)
                    metrics = json.loads(result.stdout)
                else:
                    metrics = {key: old[key] for key in ("radius_max", "sampled_max_dimension", "far_corner_dimension", "resolution_maps", "seconds") if key in old}
                record.update(metrics, image=image_path.as_posix())
            records.append(record)
            print(f"{number}/{len(inputs)} {stem}", flush=True)
    finally:
        # Retain valid cached pictures outside an experiment subset in the shared gallery.
        cached = {row["scc"]: previous[row["scc"]] for row in all_inputs if row["scc"] in previous
                  and all(previous[row["scc"]].get(key) == row[key] for key in ("pointset_sha256", "scc_sha256"))
                  and previous[row["scc"]].get("renderer_sha256") == renderer_hash
                  and previous[row["scc"]].get("point_renderer_sha256") == point_renderer_hash}
        cached.update({row["scc"]: row for row in records})
        if cached:
            write_gallery(output, [cached[row["scc"]] for row in all_inputs if row["scc"] in cached])
    print(f"Gallery: {output / 'index.html'}", flush=True)


def self_check():
    with tempfile.TemporaryDirectory(prefix="pointset-plot-check-") as temporary:
        folder = Path(temporary)
        for dimension in (2, 3, 4):
            row = dict(n_points=12, ambient_dimension=dimension, scc=f"test_12_1_H1.scc")
            data = np.column_stack((np.random.default_rng(42).random((12, dimension)), -np.linspace(0, 1, 12)))
            source, output = folder / "points.txt", folder / f"R{dimension}.png"
            np.savetxt(source, data)
            plot_pointset(source, row, output)
            assert plt.imread(output).shape[0] > 1000
            try:
                plot_pointset(source, dict(row, n_points=13), output)
            except ValueError:
                pass
            else:
                raise AssertionError("Invalid point-set shape was accepted")
    print("PASS: 2D, 3D, all six 4D projections, and input validation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=("runtime", "optimizations"), help="Optional experiment subset; default: all shared inputs")
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--hilbert-renderer", type=Path, default=ROOT / "build/experiment_hilbert")
    parser.add_argument("--self-check", action="store_true", help="Check the renderer on small 2D/3D/4D clouds")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        plot_experiment(args.experiment, args.data_dir.resolve(), args.hilbert_renderer.resolve())

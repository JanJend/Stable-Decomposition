"""Small data regression check: python3 -m experiments.support.check_data."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

import numpy as np

from . import HERE, select_samples


def main():
    sys.path.insert(0, str(HERE))
    from generate_data import (FAMILIES, TORUS_BANDWIDTHS, cluster_centers, immune_samples,
                               mixed_noise, select_torus_bandwidth, synthetic, torus, with_codensity)
    from plot_pointsets import plot_experiment

    for family in FAMILIES[:-1]:
        points = synthetic(family, 2000, 42)
        assert points.shape == (2000, 3) and np.isfinite(points).all()
        assert np.array_equal(points, synthetic(family, 2000, 42))
        assert not np.array_equal(points, synthetic(family, 2000, 43))
    sphere_radii = np.linalg.norm(synthetic("S2", 2000, 42), axis=1)
    assert np.median(np.abs(sphere_radii - 1)) > 0.25
    assert np.count_nonzero(sphere_radii > 4) > 10
    assert np.min(np.std(synthetic("gauss_clust", 2000, 42), axis=0)) > 0.9
    scales = np.array([0.6, 0.8, 0.5])
    assert not np.array_equal(cluster_centers(np.random.default_rng(42), scales),
                              cluster_centers(np.random.default_rng(43), scales))
    for seed in range(30):
        centers = cluster_centers(np.random.default_rng(seed), scales)
        assert all(np.linalg.norm(centers[i] - centers[j]) >= 3 * (scales[i] + scales[j])
                   for i, j in ((0, 1), (0, 2), (1, 2)))
    clean = torus(np.random.default_rng(42), 1000)
    assert np.allclose((np.linalg.norm(clean[:, :2], axis=1) - 2) ** 2 + clean[:, 2] ** 2, 1)
    displaced = mixed_noise(np.random.default_rng(42), clean, sigma=0, outlier_fraction=0.1)
    distances = np.linalg.norm(displaced - clean, axis=1)
    assert np.count_nonzero(distances) == 100 and distances.max() <= 1.8 + 1e-12

    points = np.array([[0, 0], [0.1, 0.1], [1, 0], [0, 1], [1, 1]])
    data, scale, _, h = with_codensity(points, 0.2)
    density = np.exp(-np.sum((points[:, None] - points[None, :]) ** 2, axis=2) / (2 * 0.2 ** 2)).mean(axis=1)
    assert np.allclose(data[:, -1], -(density - density.min()) / np.ptp(density))
    assert h * scale == 0.2
    assert np.allclose(data, with_codensity(points * 7 + 3, 1.4)[0])
    assert with_codensity(points * 7, 0.05, scaled_bandwidth=True)[3] == 0.05

    sparse = synthetic("torus_gauss", 20, 42)
    selected_width, selected_score = select_torus_bandwidth(sparse)
    # Independent direct probabilities check exclusion of self and Gaussian normalization.
    squared_distances = np.sum((sparse[:, None] - sparse[None, :]) ** 2, axis=2)
    np.fill_diagonal(squared_distances, np.inf)
    scores = []
    with np.errstate(divide="ignore"):
        for width in TORUS_BANDWIDTHS:
            probability = np.exp(-squared_distances / (2 * width ** 2)).sum(axis=1)
            probability /= (len(sparse) - 1) * (np.sqrt(2 * np.pi) * width) ** 3
            scores.append(np.log(probability).mean())
    assert selected_width == TORUS_BANDWIDTHS[np.argmax(scores)]
    assert np.isclose(selected_score, max(scores), rtol=0, atol=1e-12)
    assert select_torus_bandwidth(synthetic("torus_gauss", 320, 42))[0] < selected_width
    assert select_torus_bandwidth(sparse) == (selected_width, selected_score)

    with tempfile.TemporaryDirectory(prefix="experiment-data-check-") as temporary:
        temporary = Path(temporary)
        regions = temporary / "1.5mmRegions"
        regions.mkdir()
        for count in (6, 10):
            (regions / f"observed{count}.csv").write_text("x,y\n" + "".join(f"{i},{i*i}\n" for i in range(count)))
        for name in ("LargeHypoxicRegion", "LargeHypoxicRegion2"):
            (temporary / f"{name}.csv").write_text("x,y,Celltype\n" + "".join(
                f"{i},{i*i},{('CD8', 'CD68', 'FoxP3')[i % 3]}\n" for i in range(15)))
        observed = immune_samples(temporary, [6, 10])
        assert [len(points) for _, _, points in observed] == [6, 10, 15, 15]
        rows = [dict(dataset=path.stem, family="immuno", n_points=len(points), target_points=target)
                for target, path, points in observed]
        assert len(select_samples(rows, [6])) == 1
        assert len(select_samples(rows, [6], include_large_immune=True)) == 3

        bundle = temporary / "data"
        command = [sys.executable, str(HERE / "generate_data.py"), "--families", "torus_gauss", "S2", "gauss_clust",
                   "--sizes", "20", "--seeds", "42", "--output-dir", str(bundle)]
        subprocess.run(command, check=True)
        manifest = bundle / "manifest.csv"
        before = manifest.read_bytes()
        subprocess.run(command, check=True)
        assert manifest.read_bytes() == before, "Matching presentations should be reused"
        rows = list(csv.DictReader(manifest.open()))
        assert {(row["family"], row["homology"]) for row in rows} == {
            ("torus_gauss", "1"), ("S2", "2"), ("gauss_clust", "0")}
        torus_row = next(row for row in rows if row["family"] == "torus_gauss")
        assert torus_row["kde_bandwidth_method"] == "loo_cv"
        assert np.isclose(float(torus_row["kde_bandwidth"]) * float(torus_row["spatial_scale"]), selected_width)
        # An explicit S2 override exercises sharing a cloud between homology degrees.
        sphere_command = [sys.executable, str(HERE / "generate_data.py"), "--families", "S2", "--sizes", "20",
                          "--seeds", "42", "--output-dir", str(bundle)]
        subprocess.run(sphere_command + ["--homology", "1", "2"], check=True)
        plot_experiment(bundle=bundle)
        pictures = bundle / "pictures"
        assert len(list(pictures.rglob("*_points.png"))) == 3, "A cloud shared by H1/H2 is plotted once"
        assert len(list(pictures.rglob("*.png"))) == 7
        assert (pictures / "index.html").is_file()
        image_times = {path: path.stat().st_mtime_ns for path in pictures.rglob("*.png")}
        with patch("plot_pointsets.subprocess.run") as render:
            plot_experiment(bundle=bundle)
            plot_experiment("optimizations", bundle=bundle)
            render.assert_not_called()
        assert all(path.stat().st_mtime_ns == stamp for path, stamp in image_times.items())
        assert len(list(csv.DictReader((pictures / "index.csv").open()))) == 4
        # Changing one degree's bandwidth must invalidate its sibling manifest entry.
        subprocess.run(sphere_command + ["--homology", "2", "--bandwidth", "0.3"], check=True)
        rows = list(csv.DictReader(manifest.open()))
        assert [row["homology"] for row in rows if row["family"] == "S2"] == ["2"]
        subprocess.run(command[:command.index("--families")] + ["--families", "torus_gauss", "--sizes", "20",
                       "--seeds", "42", "--bandwidth", "0.4", "--output-dir", str(bundle)], check=True)
        row = next(row for row in csv.DictReader(manifest.open()) if row["family"] == "torus_gauss")
        assert row["kde_bandwidth_method"] == "override" and row["kde_bandwidth_loo_score"] == ""
        assert np.isclose(float(row["kde_bandwidth"]) * float(row["spatial_scale"]), 0.4)
    print("PASS: 3D samples, mixed outliers, KDE, adaptive torus bandwidth and override, immune types, reuse and plots")


if __name__ == "__main__":
    main()

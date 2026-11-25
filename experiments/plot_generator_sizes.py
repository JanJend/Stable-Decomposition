#!/usr/bin/env python3
"""Plot peak reduced I/K sizes at 0.5% epsilon through 640 points, on log-log axes."""
import argparse
import csv
import math
from pathlib import Path

from support import HERE

FAMILIES = (("torus_gauss", "Gaussian torus H1"), ("gauss_clust", "Gaussian clusters H0"),
            ("immuno", "Immune cells H1"))
PERCENT_COLUMNS = ("family", "dataset", "n_points", "epsilon_percent", "presentation_size",
                   "i_peak_reduced_generators", "k_peak_reduced_generators",
                   "i_percent", "k_percent", "i_complete", "k_complete", "status")


def read_percentages(path):
    measured, seen = [], set()
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if float(row["epsilon_percent"]) != 0.5 or int(row["target_points"]) > 640:
                continue
            n, size = int(row["n_points"]), int(row["input_generators"]) + int(row["input_relations"])
            key = row["family"], n
            if n <= 0 or size < 0 or key in seen:
                raise ValueError(f"Invalid or duplicate sample: {key}")
            seen.add(key)
            record = dict(family=row["family"], dataset=row["dataset"], n_points=n,
                          epsilon_percent=0.5, presentation_size=size, status=row["status"])
            for obj in ("i", "k"):
                count = row[f"{obj}_peak_reduced_generators"]
                count = int(count) if count else None
                if count is not None and count < 0:
                    raise ValueError("Generator counts must be nonnegative")
                record[f"{obj}_peak_reduced_generators"] = count
                record[f"{obj}_percent"] = 100 * count / size if count is not None and size else None
                record[f"{obj}_complete"] = bool(row[f"{obj}_final_generators"])
            measured.append(record)
    if not measured:
        raise ValueError("No 0.5% measurements through 640 points")
    return sorted(measured, key=lambda row: (row["family"], row["n_points"]))


def plot(measured):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, ScalarFormatter

    figure, axes = plt.subplots(1, 3, figsize=(12, 4.3), sharey=True)
    positive = [row[f"{obj}_percent"] for row in measured for obj in ("i", "k")
                if row[f"{obj}_percent"] is not None and row[f"{obj}_percent"] > 0]
    if not positive:
        raise ValueError("No positive percentages for a logarithmic y-axis")
    notes = set()
    for axis, (family, title) in zip(axes, FAMILIES):
        rows = [row for row in measured if row["family"] == family]
        x = [row["n_points"] for row in rows]
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        for obj, label, color, marker, style in (("i", "$I$", "#2463A6", "o", "-"),
                                                 ("k", "$K$", "#C85835", "s", "--")):
            y = [row[f"{obj}_percent"] if row[f"{obj}_percent"] is not None
                 and row[f"{obj}_percent"] > 0 and row[f"{obj}_complete"] else math.nan for row in rows]
            axis.plot(x, y, style, color=color, marker=marker, markersize=5, linewidth=1.7, label=label)
            for row in rows:
                value = row[f"{obj}_percent"]
                if value == 0:
                    notes.add("Zero percentages are omitted from the logarithmic axis.")
                elif value is None:
                    notes.add("Missing percentages appear as gaps.")
                elif not row[f"{obj}_complete"]:
                    axis.scatter(row["n_points"], value, marker="^", facecolors="none", edgecolors=color, zorder=4)
                    axis.annotate("", xy=(row["n_points"], value * 1.15), xytext=(row["n_points"], value),
                                  arrowprops=dict(arrowstyle="->", color=color, linewidth=1.2))
                    notes.add("Upward arrows mark lower bounds for unfinished phases.")
        axis.axhline(100, color="0.65", linewidth=0.8, linestyle=":", zorder=0)
        axis.set(title=title, xlabel="Number of points", xlim=(17, 750),
                 ylim=(min(positive) / 1.35, max(positive) * 1.35))
        axis.set_xticks([20, 40, 80, 160, 320, 640])
        axis.xaxis.set_major_formatter(ScalarFormatter())
        axis.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))
        axis.yaxis.set_minor_formatter(NullFormatter())
        axis.grid(True, which="major", color="0.9", linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$\frac{|b_0(-)|}{|b_0(X)|+|b_1(X)|}$")
    figure.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 1),
                  ncol=2, frameon=False)
    if notes:
        figure.text(0.5, 0.01, " ".join(sorted(notes)), ha="center", fontsize=8)
    figure.tight_layout(rect=(0, 0.05 if notes else 0, 1, 0.95))
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=HERE / "results/generator_sizes_to640_eps0p5_1h/results.csv")
    parser.add_argument("--output-dir", type=Path, help="Default: beside the input CSV")
    args = parser.parse_args()
    try:
        measured = read_percentages(args.input)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    directory = args.output_dir or args.input.parent
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "percentages.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=PERCENT_COLUMNS)
        writer.writeheader()
        writer.writerows(measured)
    figure = plot(measured)
    for suffix in ("png", "svg"):
        figure.savefig(directory / f"generator_sizes_loglog.{suffix}", dpi=220, facecolor="white")
    print(f"Saved {directory / 'generator_sizes_loglog.png'} and SVG; percentages.csv contains the plotted values")


if __name__ == "__main__":
    main()

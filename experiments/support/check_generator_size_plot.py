"""Run with a Matplotlib Python: python -m experiments.support.check_generator_size_plot."""
import csv
import math
from pathlib import Path
import sys
import tempfile

from . import HERE


def main():
    sys.path.insert(0, str(HERE))
    import plot_generator_sizes as experiment

    with tempfile.TemporaryDirectory(prefix="generator-size-plot-") as name:
        csv_path = Path(name) / "results.csv"
        rows = []
        for family, _ in experiment.FAMILIES:
            for n in (20, 40):
                unfinished = family == "gauss_clust" and n == 40
                rows.append(dict(family=family, dataset=f"{family}_{n}", target_points=n, n_points=n,
                                 epsilon_percent=0.5, input_generators=3, input_relations=7,
                                 i_peak_reduced_generators=n,
                                 k_peak_reduced_generators=0 if family == "immuno" and n == 20 else 5,
                                 i_final_generators="" if unfinished else 3, k_final_generators=1,
                                 status="timeout" if unfinished else "ok"))
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
        measured = experiment.read_percentages(csv_path)
        assert len(measured) == 6 and all(row["presentation_size"] == 10 for row in measured)
        assert [row["i_percent"] for row in measured] == [200, 400] * 3
        assert next(row for row in measured if row["family"] == "immuno")["k_percent"] == 0
        figure = experiment.plot(measured)
        figure.canvas.draw()
        assert all(axis.get_xscale() == axis.get_yscale() == "log" for axis in figure.axes)
        assert list(figure.axes[0].lines[0].get_ydata()) == [200, 400]
        assert math.isnan(figure.axes[1].lines[0].get_ydata()[1])  # Lower bound is a separate arrow.
        assert math.isnan(figure.axes[2].lines[1].get_ydata()[0])  # No invented positive value for zero.
        notes = " ".join(text.get_text() for text in figure.texts)
        assert "Zero percentages" in notes and "lower bounds" in notes
    print("PASS: percentages, actual point counts, log-log axes, phase-specific bounds and zero handling")


if __name__ == "__main__":
    main()

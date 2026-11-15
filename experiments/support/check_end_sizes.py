"""Small end-to-end check: python3 -m experiments.support.check_end_sizes."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from . import HERE


def main():
    with tempfile.TemporaryDirectory(prefix="end-sizes-check-") as temporary:
        temporary = Path(temporary)
        source = temporary / "input with spaces.scc"
        source.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
        output = temporary / "results"
        output.mkdir()
        (output / "results.csv").write_text("old_column\nold_result\n")
        command = [sys.executable, str(HERE / "benchmark_end_sizes.py"),
                   "--input", str(source), "--output-dir", str(output)]
        subprocess.run([*command, "--percentages", "0", "0.5", "50"], check=True)
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        expected = ["Pointset + epsilon", "Generator settings", "dim End_2d", "dim End_2d/End_0"]
        expected += [f"{space} {method} {side}" for space in ("Full", "Quotient")
                     for method in ("Sparse-nilpotent", "Complement-orbit", "Combined") for side in ("I", "K")]
        expected += ["End_2d time (s)", "End_0 time (s)", "Quotient time (s)",
                     "Full coordinates time (s)", "Quotient coordinates time (s)"]
        expected += [f"{space} {method} time (s)" for space in ("Full", "Quotient")
                     for method in ("Sparse-nilpotent", "Complement-orbit", "Combined")]
        assert list(rows[0]) == expected and len(expected) == 27
        assert [row["Pointset + epsilon"] for row in rows] == [
            "input with spaces — 0%", "input with spaces — 0.5%", "input with spaces — 50%"]
        assert [int(row["dim End_2d"]) for row in rows] == [2, 2, 4]
        assert [int(row["dim End_2d/End_0"]) for row in rows] == [0, 0, 2]
        for side in ("I", "K"):
            for method in ("Sparse-nilpotent", "Complement-orbit", "Combined"):
                assert [int(row[f"Quotient {method} {side}"]) for row in rows] == [0, 0, 1]
                assert [int(row[f"Full {method} {side}"]) for row in rows] == [1, 1, 2]
        for row in rows:
            assert all(float(value) >= 0 for key, value in row.items() if key.endswith("(s)"))
        metadata = [json.loads(path.read_text()) for path in sorted(output.glob("*/epsilon_*.json"))]
        assert [float(row["epsilon"]) for row in metadata] == [0, 0.01, 1]
        for row in metadata:
            assert row["status"] == "ok" and float(row["bounding_box_extent"]) == 2
            assert math.isclose(float(row["shift"]), 2 * float(row["epsilon"]))
            assert "--adaptive" not in row["command"]
            assert (output / row["stdout_log"]).is_file() and (output / row["stderr_log"]).is_file()
        settings = json.loads(next(output.glob("*/settings.json")).read_text())
        assert settings["timeout_seconds"] == 600 and settings["adaptive"] is False
        assert next(output.glob("*/previous_results.csv")).read_text() == "old_column\nold_result\n"
        assert len({row["stdout_log"] for row in metadata}) == 3
        subprocess.run([*command, "--percentages", "50"], check=True)
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 4 and rows[-1]["Full Combined I"] == "2"
        result = subprocess.run([*command, "--percentages", "50", "--timeout", "0.000001"])
        assert result.returncode == 1
        with (output / "results.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 5 and rows[-1]["dim End_2d"] == "timeout"
        for disabled, remaining in (("sparse-nilpotent", "Complement-orbit"), ("complement-orbit", "Sparse-nilpotent")):
            subprocess.run([*command, "--percentages", "50", f"--no-{disabled}-generators"], check=True)
            with (output / "results.csv").open() as stream:
                last = list(csv.DictReader(stream))[-1]
            for space in ("Full", "Quotient"):
                assert last[f"{space} Combined I"] == "disabled"
                assert last[f"{space} {disabled.capitalize()} K"] == "disabled"
                assert last[f"{space} {remaining} I"] == ("2" if space == "Full" else "1")
        subprocess.run([*command, "--percentages", "50", "--no-sparse-nilpotent-generators",
                        "--no-complement-orbit-generators"], check=True)
        with (output / "results.csv").open() as stream:
            last = list(csv.DictReader(stream))[-1]
        assert last["dim End_2d"] == "4" and last["Quotient Combined I"] == "disabled"
        subprocess.run([*command, "--percentages", "50", "--composition-orbit-rounds", "6",
                        "--composition-orbit-seed-mode", "random", "--composition-orbit-random-seed", "42"], check=True)
        with (output / "results.csv").open() as stream:
            tuned = list(csv.DictReader(stream))[-1]
        config = json.loads(tuned["Generator settings"])
        assert config["composition-orbit-rounds"] == 6 and config["composition-orbit-seed-mode"] == "random"
        assert config["composition-orbit-random-seed"] == 42
        assert int(tuned["Quotient Combined I"]) <= 2
        invalid = subprocess.run([*command, "--composition-orbit-rounds", "-1"], capture_output=True)
        assert invalid.returncode != 0
    print("PASS: 27 columns, independent methods, disabled cells, scaling, append, and timeout recording")


if __name__ == "__main__":
    main()

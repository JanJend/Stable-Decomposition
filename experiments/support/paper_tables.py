"""Absolute runtime tables and baseline-relative percentage speedup tables."""
import csv
import statistics

from .scaling import SWITCHES

VARIANTS = ["baseline", *SWITCHES]
LABELS = ["Base", "Proj", "Amort", "Full End", "Sparse", "Adaptive", "Reuse", "Orbit"]
SPEEDUP_NAMES = ["hom-quotient" if name == "full-end-space" else name for name in SWITCHES]
ROW_END = " " + chr(92) * 2


def aggregate(records):
    completed = [row for row in records if row["status"] == "ok"]
    result = {key: records[0].get(key, "") for key in
              ("dataset", "family", "seed", "n_points", "target_points", "homology",
               "generators", "relations", "epsilon_percent", "epsilon", "column_storage", "variant",
               "initial_maps", "i_maps", "k_maps", "generator_settings")}
    result.update(completed=len(completed), attempts=sum(row["repetition"] > 0 for row in records),
                  status="ok" if len(completed) == len(records) else
                  "partial" if completed else records[-1]["status"])
    for field in ("pruning_seconds", "process_seconds"):
        values = [row[field] for row in completed]
        result[field] = statistics.median(values) if values else ""
        result[field + "_min"] = min(values) if values else ""
        result[field + "_max"] = max(values) if values else ""
    result["timeout_seconds"] = records[0]["timeout_seconds"]
    return result


def cell(row, field="pruning_seconds"):
    if row["completed"]:
        text = f"{row[field]:.6f}"
        if row["status"] == "partial":
            text += f" [{row['completed']}/{row['attempts']}]"
        return text
    if row["status"] == "timeout":
        return f"timeout ({row['timeout_seconds']:g} s)"
    return row["status"]


def speedup_cell(baseline, variant):
    if baseline is None or variant is None:
        return "pending"
    if baseline["status"] == "pending" or variant["status"] == "pending":
        return "pending"
    if (not baseline["completed"] or not variant["completed"] or
            baseline["pruning_seconds"] <= 0 or variant["pruning_seconds"] <= 0):
        return "n/a"
    percent = 100 * (baseline["pruning_seconds"] / variant["pruning_seconds"] - 1)
    text = f"{percent:+.1f}%"
    if baseline["status"] == "partial" or variant["status"] == "partial":
        text += " [partial]"
    return text


def latex(text):
    return str(text).translate(str.maketrans({"\\": r"\textbackslash{}", "_": r"\_",
        "%": r"\%", "&": r"\&", "#": r"\#", "{": r"\{", "}": r"\}", "$": r"\$"}))


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def average_seeds(summaries, seed_counts):
    """Give each input seed equal weight, after taking its timing median.

    Never present a successful-seeds-only mean when a seed failed or is pending.
    Detailed per-seed measurements stay in summary.csv and observations.jsonl.
    """
    groups = {}
    for row in summaries:
        key = (row["family"], row["n_points"], row["homology"], row["epsilon_percent"],
               row["column_storage"], row["variant"])
        groups.setdefault(key, []).append(row)
    result = []
    for key, rows in groups.items():
        expected = seed_counts.get(key[:3], len(rows))
        mean = dict(rows[0], dataset=key[0], seed="", completed=0, attempts=expected,
                    initial_maps="", i_maps="", k_maps="")
        if len(rows) < expected:
            mean["status"] = "pending"
        elif any(not row["completed"] for row in rows):
            mean["status"] = "/".join(sorted({row["status"] for row in rows if not row["completed"]}))
        else:
            mean.update(completed=sum(row["completed"] for row in rows),
                        attempts=sum(row["attempts"] for row in rows),
                        status="partial" if any(row["status"] == "partial" for row in rows) else "ok")
            for field in ("pruning_seconds", "process_seconds", "generators", "relations", "epsilon", "initial_maps", "i_maps", "k_maps"):
                mean[field] = statistics.mean(row[field] for row in rows) if all(row[field] != "" for row in rows) else ""
        result.append(mean)
    return result


def write_reports(observations, directory, kind, *, seed_counts=None, generator_configs=None):
    if not observations:
        return
    groups = {}
    for row in observations:
        groups.setdefault((row["dataset"], row["homology"], row["epsilon_percent"], row["column_storage"], row["variant"]), []).append(row)
    summaries = [aggregate(rows) for rows in groups.values()]
    summaries.sort(key=lambda row: (row["family"], int(row["homology"]),
                                    int(row["n_points"]), str(row["seed"]), row["epsilon_percent"], row["column_storage"], row["variant"]))
    write_csv(directory / "summary.csv", summaries)
    if kind == "optimizations":
        summaries = average_seeds(summaries, seed_counts or {})
    cases = {}
    for row in summaries:
        cases.setdefault((row["dataset"], row["n_points"], row["homology"], row["epsilon_percent"], row["column_storage"]), {})[row["variant"]] = row
    matrix, speedups, generators = [], [], []
    for variants in cases.values():
        first = next(iter(variants.values()))
        fields = ("dataset", "family", "n_points", "column_storage", "generators", "relations", "epsilon_percent", "epsilon")
        if kind == "runtime":
            fields += ("homology", "seed")
        entry = {key: first[key] for key in fields}
        if kind == "runtime":
            row = variants["all-optimizations"]
            entry.update(pruning_seconds=cell(row), process_seconds=cell(row, "process_seconds"),
                         completed=row["completed"], attempts=row["attempts"], status=row["status"])
        else:
            comparison = dict(entry)
            comparison.update((column + "_speedup_percent", speedup_cell(variants.get(name), variants.get("baseline"))
                               if name == "full-end-space" else speedup_cell(variants.get("baseline"), variants.get(name)))
                              for name, column in zip(SWITCHES, SPEEDUP_NAMES))
            speedups.append(comparison)
            entry.update((name + "_seconds", cell(variants[name]) if name in variants else "pending")
                         for name in VARIANTS)
        matrix.append(entry)
        for name, settings in (generator_configs or {}).items():
            row = variants.get(name)
            generators.append(dict(
                {key: first[key] for key in fields}, variant=name, **settings,
                initial_maps=row.get("initial_maps", "") if row else "",
                i_maps=row.get("i_maps", "") if row else "", k_maps=row.get("k_maps", "") if row else "",
                pruning_seconds=cell(row) if row else "pending",
                speedup_percent=speedup_cell(variants.get("baseline"), row)))
    if generators:
        write_csv(directory / "generators.csv", generators)
        write_generator_latex(directory / "generators.tex", generators)
    write_csv(directory / "table.csv", matrix)
    write_latex(directory / "table.tex", matrix, kind)
    if kind == "optimizations":
        write_csv(directory / "speedups.csv", speedups)
        write_latex(directory / "speedups.tex", speedups, kind, percentages=True)


def compact_cell(value, percentages=False):
    """Short, legible cells for a normal-width optimization table."""
    value = str(value)
    token = value.split()[0]
    try:
        number = float(token.rstrip("%"))
    except ValueError:
        if value == "pending":
            return "--"
        if value == "n/a":
            return value
        if value.startswith("timeout"):
            return "timeout"
        if "memory_limit" in value:
            return "mem."
        if "generation_" in value:
            return "gen. err."
        return "error"
    return (f"{number:+.3g}%" if percentages else f"{number:.3g}") + ("*" if "[" in value else "")


def write_latex(path, matrix, kind, percentages=False):
    lines = ["% Generated automatically. All measured times are seconds.",
             "% Timeout is the process deadline, not a completed pruning time."]
    if kind == "optimizations":
        names = SPEEDUP_NAMES if percentages else VARIANTS
        labels = (["Quot" if name == "hom-quotient" else label for name, label in zip(names, LABELS[1:])]
                  if percentages else LABELS)
        caption = ("Optimization speedups (percent) from the mean runtimes. "
                   "Quot = 100*(Full End/Base - 1); other columns = 100*(Base/Variant - 1)."
                   if percentages else "Optimization runtimes (seconds).")
        caption += (" Synthetic entries average the three seed medians; immune entries use one observed input. "
                    "Epsilon is a percentage of degree extent. -- means pending; n/a means unavailable; "
                    "* marks partial timing repeats. See summary.csv for failures.")
        headings = ["Point set", "$n$", r"$\epsilon$ (\%)", *labels]
        columns = "@{}lrr" + "r" * len(names) + "@{}"
        lines += ["% Standard table/tabular; no page rotation, resizing or geometry changes.",
                  "% Equal-weight mean of per-seed timing medians. A failed or missing seed is not silently omitted.",
                  r"\begin{table}[tbp]", r"\centering", r"\begingroup\scriptsize",
                  r"\setlength{\tabcolsep}{2pt}", r"\renewcommand{\arraystretch}{1.05}",
                  r"\caption{" + latex(caption) + "}", r"\begin{tabular}{" + columns + "}", r"\hline",
                  " & ".join(headings) + ROW_END, r"\hline"]
        suffix = "_speedup_percent" if percentages else "_seconds"
        for entry in matrix:
            values = [entry["family"], entry["n_points"], f"{entry['epsilon_percent']:g}"]
            values += [compact_cell(entry[name + suffix], percentages) for name in names]
            lines.append(" & ".join(latex(value) for value in values) + ROW_END)
        lines += [r"\hline", r"\end{tabular}", r"\endgroup", r"\end{table}", "% Column key:"]
        lines += ["% " + label + ": " + name for label, name in zip(labels, names)]
        lines += ["% Base: Hom quotient, A-1/2, optional switches off. Full End omits the quotient.",
                  "% Dense pivots, fused I, deferred K minimization and throttled progress are always active."]
    else:
        # The separate runtime experiment retains its per-input detailed tables.
        lines += ["% Requires \\usepackage{longtable}.", "% Medians across timing repeats; seeds remain separate."]
        by_family = {}
        for entry in matrix:
            by_family.setdefault((entry["family"], entry["homology"], entry["epsilon_percent"], entry["column_storage"]), []).append(entry)
        for (family, homology, percentage, storage), entries in by_family.items():
            headings = ["$n$", "Seed", r"$\epsilon$", "Gens", "Rels", "Pruning (s)", "Process (s)", "Runs"]
            caption = f"{family}, H{homology}, {storage} columns. Epsilon: {percentage:g}% of degree extent. All optimizations enabled."
            column = r"p{\dimexpr(\linewidth-5em-4\tabcolsep)/6-2\tabcolsep\relax}"
            columns = r"*{2}{p{2.5em}}*{6}{" + column + "}"
            lines += [r"\begingroup\scriptsize\setlength{\tabcolsep}{1.5pt}",
                      r"\setlength{\LTcapwidth}{\linewidth}", r"\begin{longtable}[c]{" + columns + "}",
                      r"\caption{" + latex(caption) + "}" + ROW_END, r"\hline",
                      " & ".join(r"{\raggedleft " + heading + r"\par}" for heading in headings) + ROW_END,
                      r"\hline\endhead"]
            for entry in entries:
                values = [entry["n_points"], entry["seed"] if entry["seed"] != "" else "--",
                          f"{entry['epsilon']:.6g}" if entry["epsilon"] != "" else "--",
                          entry["generators"], entry["relations"], entry["pruning_seconds"],
                          entry["process_seconds"], f"{entry['completed']}/{entry['attempts']}"]
                lines.append(" & ".join(r"{\raggedleft " + latex(value) + r"\par}" for value in values) + ROW_END)
            lines += [r"\noalign{\nobreak\hrule height\arrayrulewidth}", r"\end{longtable}", r"\endgroup"]
    lines += ["% Fast configurations: requested repeats. Configurations above repeat-below: one observation.",
              "% See summary.csv for individual seeds, sample counts, min/max, failures and partial repetitions.",
              "% Hom algorithm comparisons are excluded; every configuration uses A-1/2."]
    path.write_text("\n".join(lines) + "\n")


def write_generator_latex(path, rows):
    """Tall comparison table; extra budgets add rows rather than widening the main table."""
    lines = [r"% Requires \usepackage{longtable}.",
             "% Generator configurations enable Sparse-nilpotent and Complement-orbit only; other optional switches are off.",
             "% Mean of per-seed medians. Full settings and failures are in generators.csv and observations.jsonl.",
             r"\begingroup\scriptsize\setlength{\tabcolsep}{3pt}",
             r"\begin{longtable}{@{}lrrlrrrrr@{}}",
             r"\caption{Composition-generator budget sweep. R is the maximum number of orbit attempts per side. "
             r"I/K are selected map counts; time is full pruning time in seconds. "
             r"Speedup is relative to the same baseline as the optimization table.}" + ROW_END,
             r"\hline",
             r"Point set & $n$ & $\epsilon$ (\%) & Mode & R & I & K & Time & Speedup" + ROW_END,
             r"\hline\endhead"]
    for row in rows:
        values = [row["family"], row["n_points"], f"{row['epsilon_percent']:g}",
                  row["composition-orbit-seed-mode"], row["composition-orbit-rounds"],
                  f"{row['i_maps']:.3g}" if isinstance(row["i_maps"], (int, float)) else "--",
                  f"{row['k_maps']:.3g}" if isinstance(row["k_maps"], (int, float)) else "--",
                  compact_cell(row["pruning_seconds"]), compact_cell(row["speedup_percent"], True)]
        lines.append(" & ".join(latex(value) for value in values) + ROW_END)
    lines += [r"\hline\end{longtable}", r"\endgroup"]
    path.write_text("\n".join(lines) + "\n")

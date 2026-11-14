"""Summary tables and figures."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from .classify import CLASS_ORDER
from .features import KOZAK_MODERATE, KOZAK_ORDER, KOZAK_STRONG, KOZAK_UNKNOWN, KOZAK_WEAK

CLASS_SUMMARY = "class_summary.tsv"
LENGTH_HIST = "length_hist.tsv"
KOZAK_TABLE = "kozak_by_class.tsv"

# Figure style: one blue hue throughout, on an off-white surface.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
BLUE = "#2a78d6"
# Ordered light-to-dark ramp for an ordered scale (weak < moderate < strong).
KOZAK_COLORS = {KOZAK_WEAK: "#86b6ef", KOZAK_MODERATE: "#2a78d6", KOZAK_STRONG: "#0d366b"}


def median_from_counts(counts: dict[int, int]) -> float:
    """Median of a distribution given as {value: count}."""
    total = sum(counts.values())
    if total == 0:
        return float("nan")
    lower_rank, upper_rank = (total - 1) // 2, total // 2
    seen = 0
    lower = upper = None
    for value in sorted(counts):
        seen += counts[value]
        if lower is None and seen > lower_rank:
            lower = value
        if seen > upper_rank:
            upper = value
            break
    return (lower + upper) / 2


def _ordered_classes(classes: Iterable[str]) -> list[str]:
    present = set(classes)
    return [c for c in CLASS_ORDER if c in present] + sorted(present - set(CLASS_ORDER))


def write_summaries(out_dir: Path, stats: Counter) -> None:
    """Write the per-class summary tables from the counts collected during a scan."""
    lengths: dict[str, dict[int, int]] = defaultdict(dict)
    kozak: dict[str, Counter] = defaultdict(Counter)
    starts: dict[str, Counter] = defaultdict(Counter)
    totals: dict[str, int] = {}
    for key, count in stats.items():
        if not isinstance(key, tuple):
            continue
        if key[0] == "class":
            totals[key[1]] = count
        elif key[0] == "length":
            lengths[key[1]][key[2]] = count
        elif key[0] == "kozak":
            kozak[key[1]][key[2]] = count
        elif key[0] == "start_codon":
            starts[key[1]][key[2]] = count

    classes = _ordered_classes(totals)

    with open(out_dir / CLASS_SUMMARY, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orf_class", "n_orfs", "median_aa", "atg_start", "near_cognate_start"])
        for orf_class in classes:
            atg = starts[orf_class]["ATG"]
            writer.writerow(
                [
                    orf_class,
                    totals[orf_class],
                    f"{median_from_counts(lengths[orf_class]):g}",
                    atg,
                    totals[orf_class] - atg,
                ]
            )

    with open(out_dir / LENGTH_HIST, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orf_class", "aa_length", "n_orfs"])
        for orf_class in classes:
            for aa_length in sorted(lengths[orf_class]):
                writer.writerow([orf_class, aa_length, lengths[orf_class][aa_length]])

    with open(out_dir / KOZAK_TABLE, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orf_class"] + KOZAK_ORDER)
        for orf_class in classes:
            writer.writerow([orf_class] + [kozak[orf_class][level] for level in KOZAK_ORDER])


def _read_table(path: Path) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_summaries(result_dirs: Iterable[str | Path]):
    """Read and add up the summary tables of one or more scan output folders."""
    lengths: dict[str, Counter] = defaultdict(Counter)
    kozak: dict[str, Counter] = defaultdict(Counter)
    for result_dir in result_dirs:
        result_dir = Path(result_dir)
        for row in _read_table(result_dir / LENGTH_HIST):
            lengths[row["orf_class"]][int(row["aa_length"])] += int(row["n_orfs"])
        for row in _read_table(result_dir / KOZAK_TABLE):
            for level in KOZAK_ORDER:
                kozak[row["orf_class"]][level] += int(row[level])
    return lengths, kozak


def _style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=0)


def plot_length_by_class(lengths: dict[str, Counter], path: Path) -> None:
    """Small multiples: the ORF length distribution of each class on a shared log axis."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = _ordered_classes(c for c in lengths if lengths[c])
    if not classes:
        return
    longest = max(max(lengths[c]) for c in classes)
    shortest = min(min(lengths[c]) for c in classes)
    # Log-spaced bins: ORF lengths span three orders of magnitude.
    low, high = math.log10(max(shortest, 1)), math.log10(longest + 1)
    n_bins = 40
    edges = [10 ** (low + (high - low) * i / n_bins) for i in range(n_bins + 1)]

    n_cols = min(3, len(classes))
    n_rows = math.ceil(len(classes) / n_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(3.4 * n_cols, 2.3 * n_rows), sharex=True, squeeze=False
    )
    fig.patch.set_facecolor(SURFACE)

    for ax, orf_class in zip(axes.flat, classes):
        counts = lengths[orf_class]
        values, weights = list(counts), list(counts.values())
        ax.hist(values, bins=edges, weights=weights, color=BLUE, edgecolor=SURFACE, linewidth=0.6)
        ax.set_xscale("log")
        _style_axes(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.ticklabel_format(axis="y", style="plain")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(_short_number))
        total = sum(weights)
        ax.set_title(
            f"{orf_class}\nn = {total:,} · median {median_from_counts(counts):g} aa",
            fontsize=9,
            color=INK,
            loc="left",
        )
    for ax in axes.flat[len(classes) :]:
        ax.set_visible(False)
    # Shared x: show tick labels on the lowest visible panel of each column.
    for col in range(n_cols):
        visible = [axes[row][col] for row in range(n_rows) if axes[row][col].get_visible()]
        if visible:
            visible[-1].tick_params(labelbottom=True)
            visible[-1].set_xlabel(
                "ORF length (amino acids, log scale)", fontsize=8, color=INK_SECONDARY
            )
    for row in axes:
        row[0].set_ylabel("ORFs", fontsize=8, color=INK_SECONDARY)

    fig.suptitle("ORF length by class", fontsize=12, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def _short_number(value: float, _position=None) -> str:
    for threshold, suffix in ((1e6, "M"), (1e3, "k")):
        if abs(value) >= threshold:
            return f"{value / threshold:g}{suffix}"
    return f"{value:g}"


def plot_kozak_by_class(kozak: dict[str, Counter], path: Path) -> None:
    """Stacked bars: share of strong, moderate and weak Kozak contexts in each class."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    levels = [KOZAK_STRONG, KOZAK_MODERATE, KOZAK_WEAK]
    classes = [c for c in _ordered_classes(kozak) if sum(kozak[c][lv] for lv in levels) > 0]
    if not classes:
        return

    fig, ax = plt.subplots(figsize=(7.5, 0.55 * len(classes) + 1.5))
    fig.patch.set_facecolor(SURFACE)
    _style_axes(ax)
    ax.spines["left"].set_visible(False)

    positions = list(range(len(classes)))[::-1]  # first class at the top
    left = [0.0] * len(classes)
    for level in levels:
        shares = []
        for orf_class in classes:
            scored = sum(kozak[orf_class][lv] for lv in levels)
            shares.append(100 * kozak[orf_class][level] / scored)
        ax.barh(
            positions,
            shares,
            left=left,
            height=0.62,
            color=KOZAK_COLORS[level],
            edgecolor=SURFACE,
            linewidth=2,
            label=level,
        )
        for position, start, share in zip(positions, left, shares):
            if share >= 7:
                ax.text(
                    start + share / 2,
                    position,
                    f"{share:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if level != KOZAK_WEAK else INK,
                )
        left = [a + b for a, b in zip(left, shares)]

    ax.set_yticks(positions)
    ax.set_yticklabels(classes, fontsize=9, color=INK)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of ORFs (%)", fontsize=8, color=INK_SECONDARY)
    ax.set_title("Kozak context strength by ORF class", fontsize=12, color=INK, loc="left", pad=28)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.0),
        ncol=3,
        frameon=False,
        fontsize=8,
        labelcolor=INK_SECONDARY,
        handlelength=1.2,
        columnspacing=1.4,
        borderaxespad=0.2,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def make_figures(result_dirs: Iterable[str | Path], out_dir: str | Path) -> list[Path]:
    """Draw both figures from one or more scan output folders."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lengths, kozak = load_summaries(result_dirs)
    length_path = out_dir / "orf_length_by_class.png"
    kozak_path = out_dir / "kozak_by_class.png"
    plot_length_by_class(lengths, length_path)
    plot_kozak_by_class(kozak, kozak_path)
    return [p for p in (length_path, kozak_path) if p.exists()]

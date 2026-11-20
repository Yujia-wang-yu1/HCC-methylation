"""Figure S2."""

from common import frame, save, run

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def figure_s2():
    mpl.rcdefaults()
    NAVY = "#244A73"
    GRID = "#E0E5E9"
    def _style_axis(ax: plt.Axes) -> None:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.75)
        ax.spines["bottom"].set_linewidth(0.75)
        ax.grid(axis="x", color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(
            axis="both",
            direction="out",
            length=3.5,
            width=0.75,
            labelsize=8.25,
        )


    rows = frame("Gene omission").to_dict("records")
    assert len(rows) == 30
    genes = [row["Omitted gene"] for row in rows]
    methylation = np.asarray(
        [float(row["Methylation score ρ"]) for row in rows],
        dtype=float,
    )
    rna = np.asarray(
        [float(row["RNA suppression score ρ"]) for row in rows],
        dtype=float,
    )
    if not np.isfinite(methylation).all() or not np.isfinite(rna).all():
        raise ValueError("Gene-omission correlations contain missing or infinite values.")

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Calibri"],
            "font.size": 8.5,
            "axes.labelsize": 9.25,
            "axes.titlesize": 10.0,
            "xtick.labelsize": 8.25,
            "ytick.labelsize": 8.25,
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=((15.0 / 2.54), (14.5 / 2.54)),
        facecolor="white",
    )


    axes[0].set_position([0.1457, 0.1078, 0.3004, 0.8236])
    axes[1].set_position([0.6503, 0.1122, 0.3004, 0.8192])

    y_positions = np.arange(len(rows))
    panels = (
        (axes[0], methylation, "Methylation score"),
        (axes[1], rna, "RNA suppression score"),
    )
    for ax, values, title in panels:
        ax.scatter(values, y_positions, s=18, color=NAVY, edgecolor="none", zorder=3)
        ax.set_yticks(y_positions, labels=genes)
        ax.set_xlim(0.990, 1.0003)
        ax.set_ylim(29.7, -0.7)
        ax.set_xticks(np.arange(0.990, 1.0001, 0.002))
        ax.set_xlabel("Spearman ρ with full score", labelpad=5)
        ax.set_title(title, loc="left", pad=7)
        _style_axis(ax)


    fig.text(0.112, 0.958, "a", ha="left", va="top", fontsize=10.0, weight="bold")
    fig.text(0.617, 0.958, "b", ha="left", va="top", fontsize=10.0, weight="bold")


    save(fig, "Figure_S2", dpi=300)


if __name__ == "__main__":
    run(figure_s2)

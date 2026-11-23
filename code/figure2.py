"""Figure 2."""

from common import frame, table, save, run

import math

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

import numpy as np

def figure2():

    mpl.rcdefaults()

    colors = {

        "ink": "#172033",

        "slate": "#526172",

        "grid": "#CCD5DF",

        "light": "#AEBBC7",

        "pale": "#EDF2F6",

        "navy": "#2F5B8A",

        "teal": "#148A8A",

        "blue": "#4C78A8",

        "amber": "#D59022",

        "rose": "#C85561",

    }

    tier_colors = {

        "Stable": colors["rose"],

        "Partially stable": colors["amber"],

        "Exploratory": colors["blue"],

        "Unsupported": colors["light"],

    }

    mpl.rcParams.update({

        "font.family": "DejaVu Sans",

        "font.size": 8,

        "axes.labelsize": 8,

        "xtick.labelsize": 7.5,

        "ytick.labelsize": 7.5,

        "axes.edgecolor": colors["ink"],

        "axes.linewidth": 0.8,

        "figure.facecolor": "white",

        "axes.facecolor": "white",

        "savefig.facecolor": "white",

        "pdf.fonttype": 42,

        "ps.fonttype": 42,

    })

    def to_float(row, key):

        value = row.get(key, "")

        if value is None or str(value).strip() == "":

            return np.nan

        return float(str(value).strip())

    def to_int(row, key):

        return int(round(to_float(row, key)))

    def format_axis(ax, grid_axis="both"):

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        if grid_axis in ("both", "x"):

            ax.grid(axis="x", color=colors["grid"], linewidth=0.55, zorder=0)

        if grid_axis in ("both", "y"):

            ax.grid(axis="y", color=colors["grid"], linewidth=0.55, zorder=0)

        ax.tick_params(

            direction="out",

            length=3,

            width=0.7,

            colors=colors["slate"]

        )

    def panel_mark(ax, label):

        ax.text(

            -0.12, 1.08, label,

            transform=ax.transAxes,

            fontsize=13,

            fontweight="bold",

            color="black",

            ha="left",

            va="top",

            clip_on=False

        )

    genome = frame("Genome-wide").to_dict("records")

    counts = frame("Evidence counts")

    tiers = counts[counts.Panel.eq("b")].rename(columns={"Category":"tier", "Retained (n)":"n"}).to_dict("records")

    cnv = counts[counts.Panel.eq("c")].rename(columns={"Category":"stage", "Retained (n)":"n"}).to_dict("records")

    retention = counts[counts.Panel.eq("d")].rename(columns={"Category":"evidence_unit", "Retained (n)":"retained_n", "Total (n)":"starting_n"}).to_dict("records")

    composition = counts[counts.Panel.eq("e")].rename(columns={"Category":"class", "Retained (n)":"n"}).to_dict("records")

    tier_order = ["Stable", "Partially stable", "Exploratory", "Unsupported"]

    tier_counts = {row["tier"]: to_int(row, "n") for row in tiers}

    cnv_counts = {row["stage"]: to_int(row, "n") for row in cnv}

    retention_map = {row["evidence_unit"]: row for row in retention}

    composition_counts = {row["class"]: to_int(row, "n") for row in composition}

    composition_order = [

        "Composition-robust",

        "Composition-sensitive",

        "Not fully robust",

    ]

    fig = plt.figure(figsize=(8.6, 6.0), facecolor="white", constrained_layout=True)

    outer = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0])

    top = outer[0].subgridspec(1, 2, width_ratios=[2.35, 1.0], wspace=0.28)

    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.05, 1.18, 1.30], wspace=0.25)

    ax_a = fig.add_subplot(top[0, 0])

    ax_b = fig.add_subplot(top[0, 1])

    ax_c = fig.add_subplot(bottom[0, 0])

    ax_d = fig.add_subplot(bottom[0, 1])

    ax_e = fig.add_subplot(bottom[0, 2])

    rho = np.array([to_float(row, "Methylation–RNA ρ") for row in genome])

    fdr_y = np.array([to_float(row, "−log₁₀(FDR)") for row in genome])

    gene_tiers = [row["Evidence tier"] for row in genome]

    for tier in reversed(tier_order):

        mask = np.array([x == tier for x in gene_tiers])

        ax_a.scatter(

            rho[mask],

            fdr_y[mask],

            s=8,

            color=tier_colors[tier],

            alpha=0.62,

            linewidths=0,

            rasterized=True,

            label=tier,

            zorder=2

        )

    ax_a.axvline(0, color=colors["slate"], linewidth=0.8, zorder=1)

    ax_a.axhline(-math.log10(0.05), color=colors["slate"], linewidth=0.8, linestyle="--", zorder=1)

    ax_a.set_xlim(-0.9, 0.8)

    ax_a.set_ylim(-4, 91)

    ax_a.set_xticks(np.arange(-0.8, 0.81, 0.2))

    ax_a.set_yticks([0, 20, 40, 60, 80])

    ax_a.set_xlabel("TSS200 methylation–RNA Spearman ρ")

    ax_a.set_ylabel("−log10(BH FDR)")

    reference_genes = set(table("S3")["gene"])
    if len(reference_genes) != 30:
        raise ValueError("Table S3 must contain the 30 fixed reference genes.")
    sel = np.array([row["Gene"] in reference_genes for row in genome])

    ax_a.scatter(rho[sel], fdr_y[sel], s=30, facecolors="none", edgecolors="black", linewidths=0.65, label="MISP genes", zorder=4)

    ax_a.legend(

        loc="upper right",

        frameon=False,

        fontsize=7.2,

        markerscale=1.5,

        handletextpad=0.5,

        borderaxespad=0.2

    )

    format_axis(ax_a, "both")

    panel_mark(ax_a, "a")

    y_b = np.arange(len(tier_order))[::-1]

    counts_b = [tier_counts[tier] for tier in tier_order]

    colors_b = [tier_colors[tier] for tier in tier_order]

    ax_b.barh(y_b, counts_b, color=colors_b, height=0.55, zorder=2)

    ax_b.set_yticks(y_b)

    ax_b.set_yticklabels(tier_order)

    ax_b.set_xlim(0, 4700)

    ax_b.set_xlabel("Genes")

    for y0, n0 in zip(y_b, counts_b):

        ax_b.text(

            n0 + 85,

            y0,

            f"{n0:,}",

            ha="left",

            va="center",

            fontsize=7.2,

            color=colors["ink"]

        )

    format_axis(ax_b, "x")

    panel_mark(ax_b, "b")

    cnv_x = [0.25, 1.0]

    cnv_y = [

        cnv_counts["Negative ρ + FDR"],

        cnv_counts["Retained after CNV adjustment"]

    ]

    ax_c.plot(cnv_x, cnv_y, color=colors["light"], linewidth=2.2, zorder=1)

    ax_c.scatter(

        cnv_x,

        cnv_y,

        s=48,

        color=[colors["navy"], colors["teal"]],

        edgecolors="white",

        linewidths=0.4,

        zorder=2

    )

    ax_c.set_xticks(cnv_x)

    ax_c.set_xticklabels([

        "Negative ρ\n+ FDR",

        "Retained after\nCNV adjustment"

    ])

    ax_c.set_ylabel("Genes")

    ax_c.set_xlim(0, 1.05)

    ax_c.set_ylim(0, 3700)

    ax_c.set_yticks([0, 500, 1000, 1500, 2000, 2500, 3000, 3500])

    for x0, y0 in zip(cnv_x, cnv_y):

        ax_c.text(

            x0,

            y0 + 120,

            f"{y0:,}",

            ha="center",

            va="bottom",

            fontsize=7.2,

            color=colors["ink"]

        )

    format_axis(ax_c, "y")

    panel_mark(ax_c, "c")

    stable_row = retention_map["Stable genes"]

    pathway_row = retention_map["Compact immune pathways"]

    rows_d = [

        ("Immune\npathway tests", pathway_row, colors["amber"]),

        ("Stable genes", stable_row, colors["teal"]),

    ]

    y_d = [1, 0]

    for y0, (_, row, fill_color) in zip(y_d, rows_d):

        start_n = to_int(row, "starting_n")

        keep_n = to_int(row, "retained_n")

        frac = keep_n / start_n

        ax_d.barh([y0], [1.0], color=colors["pale"], height=0.45, zorder=1)

        ax_d.barh([y0], [frac], color=fill_color, height=0.45, zorder=2)

        ax_d.text(

            frac / 2,

            y0,

            f"{keep_n:,}/{start_n:,}",

            ha="center",

            va="center",

            fontsize=7.0,

            color="black"

        )

    ax_d.set_xlim(0, 1)

    ax_d.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    ax_d.set_yticks(y_d)

    ax_d.set_yticklabels([item[0] for item in rows_d])

    ax_d.set_xlabel("Retained fraction")

    format_axis(ax_d, "x")

    panel_mark(ax_d, "d")

    e_colors = {

        "Composition-robust": colors["teal"],

        "Composition-sensitive": colors["amber"],

        "Not fully robust": colors["slate"],

    }

    draw_widths = {

        "Composition-robust": composition_counts["Composition-robust"],

        "Composition-sensitive": composition_counts["Composition-sensitive"],

        "Not fully robust": composition_counts["Not fully robust"],

    }

    left = 0

    bar_y = 0.14

    bar_h = 0.22

    for item in composition_order:

        true_n = composition_counts[item]

        draw_n = draw_widths[item]

        fill_color = e_colors[item]

        ax_e.barh(

            [bar_y],

            [draw_n],

            left=left,

            color=fill_color,

            edgecolor="white",

            linewidth=1.2,

            height=bar_h,

            zorder=2

        )

        if item == "Composition-robust":

            ax_e.text(

                left + draw_n / 2,

                bar_y,

                f"{true_n:,}",

                ha="center",

                va="center",

                fontsize=7.0,

                color="white"

            )

        elif item == "Composition-sensitive":

            ax_e.text(

                left + draw_n / 2,

                bar_y + 0.20,

                f"{true_n:,}",

                ha="center",

                va="center",

                fontsize=7.0,

                color="black"

            )

        left += draw_n

    x7_start = (

        composition_counts["Composition-robust"]

        + draw_widths["Composition-sensitive"]

    )

    x7_center = x7_start + composition_counts["Not fully robust"] / 2

    ax_e.text(

        x7_center + 95,

        bar_y,

        f"{composition_counts['Not fully robust']:,}",

        ha="left",

        va="center",

        fontsize=7.0,

        color=colors["ink"]

    )

    legend_handles = [

        Line2D(

            [0], [0],

            marker="s",

            linestyle="None",

            markerfacecolor=e_colors[item],

            markeredgecolor=e_colors[item],

            markersize=6,

            label=item

        )

        for item in composition_order

    ]

    ax_e.legend(

        handles=legend_handles,

        loc="upper right",

        frameon=False,

        fontsize=6.2,

        handletextpad=0.4,

        borderaxespad=0.1

    )

    ax_e.set_xlim(0, 2450)

    ax_e.set_ylim(0, 0.60)

    ax_e.set_yticks([])

    ax_e.set_xlabel("Stable genes")

    format_axis(ax_e, "x")

    panel_mark(ax_e, "e")

    save(fig, "Figure_2", dpi=500)

if __name__ == "__main__":

    run(figure2)


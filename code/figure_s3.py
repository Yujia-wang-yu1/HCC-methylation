"""Figure S3."""

from common import frame, save, run

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

import numpy as np

def figure_s3():

    mpl.rcdefaults()

    def _number(value: str) -> float:

        try:

            return float(value)

        except (ValueError, TypeError):

            return np.nan

    def _clean(axis) -> None:

        axis.spines[["top", "right"]].set_visible(False)

        axis.grid(axis="x", color="#E0E5E9", lw=0.5)

        axis.set_axisbelow(True)

    def _mark(axis, label: str, title: str) -> None:

        axis.set_title(f"{label}   {title}", loc="left", fontweight="bold", pad=10)

    records = frame("Immune overlap").to_dict("records")

    assert len(records) == 22

    colors = {

        "navy": "#244A73",

        "rose": "#B6445C",

        "gray": "#65717B",

        "light": "#D9E1E7",

    }

    plt.rcParams.update(

        {

            "font.family": "DejaVu Sans",

            "font.size": 8,

            "axes.labelsize": 8,

            "axes.titlesize": 9,

            "xtick.labelsize": 7,

            "ytick.labelsize": 7,

            "pdf.fonttype": 42,

            "ps.fonttype": 42,

            "svg.fonttype": "none",

            "savefig.facecolor": "white",

        }

    )

    figure = plt.figure(figsize=(190 / 25.4, 170 / 25.4))

    grid = figure.add_gridspec(

        1,

        2,

        width_ratios=[1, 1.22],

        left=0.40,

        right=0.98,

        top=0.93,

        bottom=0.17,

        wspace=0.28,

    )

    panel_a = figure.add_subplot(grid[0])

    panel_b = figure.add_subplot(grid[1])

    y_positions = np.arange(len(records))

    labels = [

        record["Immune measure"]

        .replace("_", " ")

        .replace("ESTIMATE ", "")

        .replace("MCPcounter ", "MCP ")

        for record in records

    ]

    for axis in [panel_a, panel_b]:

        axis.set(yticks=y_positions, ylim=(len(records) - 0.4, -1))

        _clean(axis)

    panel_a.set_yticklabels(labels, fontsize=7)

    panel_b.set_yticklabels([])

    panel_a.barh(

        y_positions,

        [_number(record["Original genes (n)"]) for record in records],

        color=colors["light"],

        height=0.7,

        label="Original list",

    )

    panel_a.barh(

        y_positions,

        [_number(record["Shared genes (n)"]) for record in records],

        color=colors["rose"],

        height=0.7,

        label="Shared with MISP",

    )

    panel_a.set_xlabel("Genes")

    _mark(panel_a, "a", "Gene overlap")

    for y_position, record in zip(y_positions, records):

        for column, offset, color, marker in [

            ("TCGA original-score ρ", -0.12, colors["navy"], "o"),

            ("TCGA non-overlap ρ", 0.12, colors["rose"], "D"),

        ]:

            value = _number(record[column])

            valid = np.isfinite(value) and (

                column == "TCGA original-score ρ"

                or record["Independent evaluation"] == "Evaluable"

            )

            if valid:

                panel_b.scatter(value, y_position + offset, s=18, c=color, marker=marker)

        if record["Independent evaluation"] != "Evaluable":

            panel_b.text(0.05, y_position, "NIE", va="center", fontsize=6.5)

    panel_b.set(xlim=(-1, 0.35), xlabel="TCGA Spearman ρ")

    panel_b.axvline(0, c=colors["gray"], lw=0.6)

    _mark(panel_b, "b", "TCGA associations")

    panel_a.legend(

        frameon=False,

        loc="upper left",

        bbox_to_anchor=(-1.0, -0.06),

        fontsize=7,

    )

    panel_b.legend(

        handles=[

            Line2D(

                [],

                [],

                marker="o",

                ls="",

                color=colors["navy"],

                label="Original frozen score",

            ),

            Line2D(

                [],

                [],

                marker="D",

                ls="",

                color=colors["rose"],

                label="Non-overlap reconstruction",

            ),

        ],

        frameon=False,

        loc="upper left",

        bbox_to_anchor=(-0.25, -0.06),

        fontsize=7,

    )

    panel_a.set_position([.2376, .0946, .2468, .8324])

    panel_b.set_position([.5599, .0939, .2796, .8392])

    for ax in (panel_a, panel_b):

        for text in ax.get_xticklabels()+ax.get_yticklabels():

            text.set_fontfamily("Calibri")

            text.set_fontsize(8)

        ax.xaxis.label.set_fontfamily("Calibri")

        ax._left_title.set_fontfamily("Calibri")

    panel_b.tick_params(axis="x", direction="inout")

    for j, text in enumerate(panel_b.texts):

        text.set_transform(figure.transFigure)

        text.set_position((.7774, 1-(.5713+j*.03365)))

        text.set_fontfamily("Calibri")

        text.set_fontsize(8)

    panel_b.set_xticks([-1,-.75,-.5,-.25,0,.25])

    panel_a.legend(frameon=False, loc="upper left", bbox_to_anchor=(.349,.953),

                   bbox_transform=figure.transFigure, fontsize=9)

    panel_b.legend(handles=[

        Line2D([],[],marker="o",ls="",color=colors["navy"],label="Original score"),

        Line2D([],[],marker="D",ls="",color=colors["rose"],label="Non-overlap reconstruction"),

    ], frameon=False,loc="upper left",bbox_to_anchor=(.73677,.953),

       bbox_transform=figure.transFigure,fontsize=9)

    for ax in (panel_a,panel_b):

        for text in ax.get_legend().get_texts():

            text.set_fontfamily("Calibri")

    panel_b.lines[0].set(color="#8998a4", alpha=.7, linestyle="-")

    save(figure,"Figure_S3",dpi=600)

if __name__ == "__main__":

    run(figure_s3)


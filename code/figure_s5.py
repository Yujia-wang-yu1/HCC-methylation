"""Figure S5."""

from common import frame, save, run

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap

from matplotlib.colors import Normalize

from matplotlib.colors import TwoSlopeNorm

from matplotlib.patches import Patch

import numpy as np

def figure_s5():

    mpl.rcdefaults()

    methylation = frame("Promoter methylation").set_index("Sample")

    rna = frame("RNA expression").set_index("Sample")

    tracks = frame("Heatmap tracks").set_index("Sample")

    assert methylation.shape == rna.shape == (371,30)

    assert methylation.index.equals(rna.index) and methylation.index.equals(tracks.index)

    assert methylation.isna().sum().sum() == 122

    data = {

        "genes": methylation.columns.to_numpy(),

        "methylation": methylation.to_numpy().T,

        "rna_z": rna.to_numpy().T,

        "methylation_score": tracks["Methylation score"].to_numpy(),

        "rna_suppression_score": tracks["RNA suppression score"].to_numpy(),

        "purity": tracks["ESTIMATE purity"].to_numpy(),

    }

    genes = data["genes"].tolist()

    mpl.rcParams.update(

        {

            "font.family": "Calibri",

            "font.size": 9.0,

            "axes.labelsize": 9.0,

            "axes.titlesize": 10.5,

            "xtick.labelsize": 8.25,

            "ytick.labelsize": 8.25,

            "legend.fontsize": 8.25,

            "text.color": "#243643",

            "axes.labelcolor": "#243643",

            "axes.linewidth": 0.75,

            "axes.edgecolor": "#697981",

            "xtick.direction": "in",

            "ytick.direction": "in",

            "xtick.minor.visible": False,

            "ytick.minor.visible": False,

            "legend.frameon": False,

            "pdf.fonttype": 42,

            "ps.fonttype": 42,

            "svg.fonttype": "none",

            "savefig.facecolor": "white",

        }

    )

    gray = "#697981"

    missing = "#B9BFC4"

    beta_cmap = LinearSegmentedColormap.from_list(

        "methylation", ["#F5F8F9", "#B7D4DD", "#21677F"]

    )

    beta_cmap.set_bad(missing)

    rna_cmap = LinearSegmentedColormap.from_list(

        "expression", ["#315F8A", "#F7F7F4", "#B85A51"]

    )

    rna_cmap.set_bad(missing)

    score_cmap = LinearSegmentedColormap.from_list(

        "score", ["#68439A", "#F7F7F4", "#F28C00"]

    )

    purity_cmap = LinearSegmentedColormap.from_list(

        "purity", ["#F2F4F4", "#83AEA0", "#286855"]

    )

    fig = plt.figure(

        figsize=(340 / 25.4, 180 / 25.4),

        facecolor="white",

    )

    titles = ("Promoter methylation", "RNA expression")

    subtitles = ("TSS200 gene-level β values", "Gene-wise z score of log2(TPM + 1)")

    panel_lefts = (0.1006, 0.5737)

    panel_widths = (0.3733, 0.3733)

    panel_letter_positions = ((0.0597, 0.9454), (0.5328, 0.9482))

    panel_title_positions = ((0.0997, 0.9454), (0.5728, 0.9482))

    panel_subtitle_positions = ((0.0997, 0.9154), (0.5728, 0.9182))

    track_positions = (

        ((0.1006, 0.8705), (0.1006, 0.8416), (0.1006, 0.8127)),

        ((0.5737, 0.8733), (0.5728, 0.8433), (0.5737, 0.8125)),

    )

    track_label_positions = (

        ((0.0886, 0.8795), (0.0903, 0.8521), (0.0859, 0.8217)),

        ((0.5617, 0.8823), (0.5653, 0.8490), (0.5607, 0.8215)),

    )

    panel_axes_positions = (

        (0.1006, 0.2005, 0.3733, 0.59),

        (0.5737, 0.2033, 0.3733, 0.59),

    )

    panel_xlabel_positions = ((0.2706, 0.1716), (0.7287, 0.1744))

    panel_colorbar_positions = (

        (0.3146, 0.1290, 0.16, 0.014),

        (0.8014, 0.1252, 0.1455, 0.014),

    )

    panel_matrices = (data["methylation"], data["rna_z"])

    panel_cmaps = (beta_cmap, rna_cmap)

    for panel, (left, width) in enumerate(zip(panel_lefts, panel_widths)):

        fig.text(*panel_letter_positions[panel], ("a", "b")[panel], fontsize=10.5, fontweight="bold")

        fig.text(*panel_title_positions[panel], titles[panel], fontsize=10.5, fontweight="bold")

        fig.text(*panel_subtitle_positions[panel], subtitles[panel], fontsize=8.25, color=gray)

        tracks = (

            (data["methylation_score"], "Methylation score", score_cmap, TwoSlopeNorm(0, -2, 2)),

            (

                data["rna_suppression_score"],

                "RNA suppression score",

                score_cmap,

                TwoSlopeNorm(0, -2, 2),

            ),

            (data["purity"], "Purity", purity_cmap, Normalize(0, 1)),

        )

        for row, (track_values, label, cmap, norm) in enumerate(tracks):

            x_position, y_position = track_positions[panel][row]

            track = fig.add_axes([x_position, y_position, width, 0.018])

            track.imshow(

                track_values[None, :],

                aspect="auto",

                interpolation="nearest",

                cmap=cmap,

                norm=norm,

            )

            track.set_xticks([])

            track.set_yticks([])

            for spine in track.spines.values():

                spine.set_visible(False)

            fig.text(

                *track_label_positions[panel][row],

                label,

                ha="right",

                va="center",

                fontsize=8.25,

                color=gray,

            )

        ax = fig.add_axes(panel_axes_positions[panel])

        image = ax.imshow(

            np.ma.masked_invalid(panel_matrices[panel]),

            aspect="auto",

            interpolation="nearest",

            cmap=panel_cmaps[panel],

            vmin=0 if panel == 0 else -3,

            vmax=1 if panel == 0 else 3,

        )

        ax.set(

            yticks=np.arange(30),

            yticklabels=genes,

            xticks=[0, 99, 199, 299, 370],

            xticklabels=["1", "100", "200", "300", "371"],

        )

        ax.tick_params(axis="y", length=0, labelsize=8.25, pad=4)

        ax.tick_params(axis="x", length=3, width=0.75, labelsize=8.25)

        for spine in ax.spines.values():

            spine.set_visible(False)

        ax.set_xlabel(

            "Methylation-score order",

            fontsize=8.25,

            labelpad=3,

        )

        ax.xaxis.set_label_coords(

            *panel_xlabel_positions[panel],

            transform=fig.transFigure,

        )

        colorbar_ax = fig.add_axes(panel_colorbar_positions[panel])

        colorbar = fig.colorbar(

            image,

            cax=colorbar_ax,

            orientation="horizontal",

            ticks=[0, 0.5, 1] if panel == 0 else [-3, 0, 3],

            extend="neither" if panel == 0 else "both",

        )

        colorbar.outline.set_visible(False)

        colorbar.ax.tick_params(length=2, width=0.75, labelsize=8.25)

        colorbar.set_label(

            "Promoter methylation (β)"

            if panel == 0

            else "RNA expression (z score)",

            fontsize=8.25,

            labelpad=3,

        )

    fig.legend(

        handles=[Patch(facecolor=missing, label="Missing methylation(122 cells)")],

        loc="center left",

        bbox_to_anchor=(0.309585, 0.035693),

        frameon=False,

        fontsize=8.25,

        handlelength=1.0,

    )

    panel_key_delta = panel_lefts[1] - panel_lefts[0]

    key_rows = (

        (

            0.1034,

            0.1316,

            score_cmap,

            TwoSlopeNorm(0, -2, 2),

            "Methylation score",

            [-2, 0, 2],

        ),

        (

            0.1031,

            0.0844,

            score_cmap,

            TwoSlopeNorm(0, -2, 2),

            "RNA suppression score",

            [-2, 0, 2],

        ),

        (

            0.1026,

            0.0459,

            purity_cmap,

            Normalize(0, 1),

            "Purity",

            [0, 1],

        ),

    )

    for panel_delta in (0.0, panel_key_delta):

        for x_position, y_position, cmap, norm, label, ticks in key_rows:

            key_x = x_position + panel_delta

            key = fig.add_axes([key_x, y_position, 0.085, 0.014])

            colorbar = fig.colorbar(

                plt.cm.ScalarMappable(norm=norm, cmap=cmap),

                cax=key,

                orientation="horizontal",

                ticks=ticks,

            )

            colorbar.outline.set_visible(False)

            colorbar.ax.tick_params(length=2, width=0.75, labelsize=8.25)

            fig.text(

                key_x + 0.093,

                y_position + 0.007,

                label,

                ha="left",

                va="center",

                fontsize=8.25,

            )

    save(fig,"Figure_S5",dpi=300)

if __name__ == "__main__":

    run(figure_s5)


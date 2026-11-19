"""Figure 4."""

from common import frame, save, run

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


def figure4():
    mpl.rcdefaults()
    d = frame("Gene correlations")
    genes = d["Gene"].tolist()
    values = d[["TCGA ρ", "GEO ρ"]].to_numpy()
    direction = np.sign(values[:, 0]) == np.sign(values[:, 1])
    support = d["GEO FDR"].to_numpy() < .05
    pairs = frame("Paired scores").query("Cohort == 'GEO'")
    x = pairs["Methylation score"].to_numpy()
    y = pairs["RNA suppression score"].to_numpy()
    result = frame("GEO coupling").iloc[0]
    rho, p = result["ρ"], result["P"]
    assert len(genes) == 30 and len(x) == 205
    assert sum(direction) == 28 and sum(support) == 20
    check = pd.Series(x).rank().corr(pd.Series(y).rank())
    assert np.isclose(check, rho, rtol=0, atol=1e-13)
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Calibri', 'Arial', 'DejaVu Sans'],
        'font.size': 9.5, 'axes.labelsize': 10, 'xtick.labelsize': 9.5,
        'ytick.labelsize': 9.5, 'legend.fontsize': 9.5, 'text.color': '#243643',
        'axes.labelcolor': '#243643', 'axes.edgecolor': '#697981',
        'axes.linewidth': .75, 'legend.frameon': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'savefig.facecolor': 'white',
    })
    ink, gray, teal = '#243643', '#697981', '#367F80'
    cmap = LinearSegmentedColormap.from_list('association', ['#315F8A', '#F7F7F4', '#B85A51'])
    fig = plt.figure(figsize=(183 / 25.4, 156 / 25.4), facecolor='white')
    fig.text(.028, .952, 'a', fontsize=12, fontweight='bold')
    fig.text(.082, .952, 'Program-level coupling', fontsize=11, fontweight='bold')
    fig.text(.485, .952, 'b', fontsize=12, fontweight='bold')
    fig.text(.537, .952, 'Gene-level correlations', fontsize=11, fontweight='bold')


    matrix_bottom, matrix_height = .13, .735
    scatter_width = .365
    scatter_height = scatter_width * 183 / 156
    scatter_bottom = matrix_bottom + (matrix_height - scatter_height) / 2
    ax = fig.add_axes([.092, scatter_bottom, scatter_width, scatter_height])
    points = ax.scatter(x, y, s=15, color=teal, alpha=.63, edgecolors='white', linewidths=.35, rasterized=False)
    ax.set_xlabel('MISP methylation score', labelpad=7)
    ax.set_ylabel('MISP RNA suppression score', labelpad=7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(direction='out', length=3)
    ax.margins(.07)
    fig.text(.095, .797, f'n = {len(x)} paired tumours\nSpearman ρ = {rho:.3f}; P = 2.45 × 10⁻⁶',
             va='top', fontsize=9.5, linespacing=1.5)
    fig.text(.095, .190, 'Exact GEO HCC–HCC pairs\nGSE56588 × GSE63898',
             va='top', fontsize=9.5, color=gray, linespacing=1.6)

    heat = fig.add_axes([.608, matrix_bottom, .155, matrix_height])
    im = heat.imshow(values, cmap=cmap, vmin=-1, vmax=1, aspect='auto', interpolation='nearest')
    heat.set(xticks=[0, 1], xticklabels=['TCGA', 'GEO'], yticks=np.arange(30), yticklabels=genes)
    heat.xaxis.tick_top()
    heat.tick_params(axis='both', length=0, pad=5, labelsize=9.5)
    heat.set_xticks(np.arange(-.5, 2, 1), minor=True)
    heat.set_yticks(np.arange(-.5, 30, 1), minor=True)
    heat.grid(which='minor', color='white', linewidth=.65)
    heat.tick_params(which='minor', length=0)
    for spine in heat.spines.values():
        spine.set_visible(False)
    for i in range(30):
        for j in range(2):
            value = values[i, j]
            label = f'{value:.2f}' if abs(value) >= .005 else f'{value:.1e}'
            heat.text(j, i, label, ha='center', va='center', fontsize=9,
                      color='white' if abs(value) > .6 else ink)

    flags_left, flags_width = .785, .175
    flags = fig.add_axes([flags_left, matrix_bottom, flags_width, matrix_height])
    flags.set(xlim=(-.5, 1.5), ylim=(29.5, -.5), xticks=[0, 1],
              xticklabels=['Same\ndirection', 'GEO\nFDR < 0.05'], yticks=[])
    flags.xaxis.tick_top()
    flags.tick_params(axis='x', length=0, pad=5, labelsize=9)
    for spine in flags.spines.values():
        spine.set_visible(False)
    for j, arr in enumerate([direction, support]):
        flags.scatter(np.full(30, j), np.arange(30), s=25,
                      facecolors=[teal if yes else 'white' for yes in arr],
                      edgecolors=[teal if yes else '#9AA5AD' for yes in arr], linewidth=.75, zorder=3)
    for boundary in np.arange(-.5, 30, 1):
        flags.axhline(boundary, color='#EDF0F2', lw=.45, zorder=0)
    for fraction, label in [(.25, '28 / 30'), (.75, '20 / 30')]:
        fig.text(flags_left + flags_width * fraction, .106, label,
                 ha='center', fontsize=9.5, color=teal)
    cax = fig.add_axes([.608, .071, .155, .012])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=2, labelsize=9)
    cb.set_label('Spearman ρ', fontsize=9.5, labelpad=3)
    fig.legend(handles=[
        Line2D([], [], marker='o', ls='', color=teal, markersize=5, label='Criterion met'),
        Line2D([], [], marker='o', ls='', color='#9AA5AD', markerfacecolor='white', markersize=5, label='Not met'),
    ], loc='lower left', bbox_to_anchor=(.78, .025), fontsize=9, handletextpad=.5, borderaxespad=0)


    save(fig, "Figure_4", dpi=600)


if __name__ == "__main__":
    run(figure4)

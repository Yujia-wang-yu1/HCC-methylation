"""Figure 1. Cohort design and fixed-score evaluation."""
import common
from common import frame, run, save
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.transforms import Affine2D


INK = '#293A46'
BLUE = '#3C6B8F'
TEAL = '#32827D'
PLUM = '#7D668C'
MUTED = '#697B87'
TEXTS = []


def label(ax, x, y, value, size=7.6, weight='normal', color=INK, ha='left'):
    item = ax.text(x, y, value, fontsize=size, fontweight=weight,
                   color=color, ha=ha, va='center', linespacing=1.38, zorder=5)
    TEXTS.append(item)
    return item


def panel(ax, x, y, w, h, face, edge='none'):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle='round,pad=0,rounding_size=1.6',
        linewidth=.6, edgecolor=edge, facecolor=face, zorder=0))


def arrow(ax, a, b):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>',
        mutation_scale=8.3, linewidth=.95, color='#8B9BA5',
        shrinkA=0, shrinkB=0, zorder=1))


def dna(ax, x, y, height=17, width=5.3, color=BLUE, marks=True):
    # Geometric icon only: no sequence or locus-level methylation is encoded.
    t = np.linspace(0, 1, 160)
    yy = y + height * (t - .5)
    xx = width / 2 * np.cos(3 * np.pi * t)
    ax.plot(x + xx, yy, color=color, lw=1.2)
    ax.plot(x - xx, yy, color=color, lw=1.2, alpha=.6)
    for f in np.linspace(.06, .94, 10):
        w = width / 2 * np.cos(3 * np.pi * f)
        ax.plot([x-w, x+w], [y+height*(f-.5)]*2,
                color=color, lw=.65, alpha=.7)
    if marks:
        for f, sign in [(.16, -1), (.5, 1), (.84, -1)]:
            yy0 = y + height*(f-.5)
            xx0 = x + sign * abs(width/2*np.cos(3*np.pi*f))
            xx1 = x + sign*(width/2+1.4)
            ax.plot([xx0, xx1], [yy0, yy0+.7], color=color, lw=.7)
            ax.add_patch(Circle((xx1, yy0+.7), .6, color=color))


def rna(ax, x, y, height=17, color=TEAL):
    t = np.linspace(0, 1, 140)
    xx = x + 2.4 * np.sin(2.6*np.pi*t + .3)
    yy = y + height*(t-.5)
    ax.plot(xx, yy, color=color, lw=1.35)
    for f in np.linspace(.12, .88, 8):
        xx0 = x + 2.4*np.sin(2.6*np.pi*f+.3)
        yy0 = y + height*(f-.5)
        ax.plot([xx0, xx0+1.4], [yy0, yy0+.8], color=color, lw=.8)


def gene_set(ax):
    # Exactly 30 uniform symbols; no value-dependent colour or size.
    for row in range(5):
        for col in range(6):
            ax.add_patch(FancyBboxPatch((79+col*2.7, 91+row*2.7), 1.95, 1.95,
                boxstyle='round,pad=0,rounding_size=.3',
                facecolor='#B9DCD4', edgecolor=TEAL, linewidth=.45))
    # Padlock signifies the unchanged programme definition and weights.
    ax.add_patch(Arc((99.8, 97), 4.2, 5.0, theta1=0, theta2=180,
                     color=TEAL, linewidth=1.2))
    ax.add_patch(FancyBboxPatch((96.9, 92.6), 5.8, 5.6,
        boxstyle='round,pad=0,rounding_size=.8', facecolor=TEAL, edgecolor='none'))
    ax.add_patch(Circle((99.8, 95.8), .48, color='white'))
    ax.plot([99.8, 99.8], [94.5, 95.6], color='white', lw=.9)

def gene_icon(ax, target):
    existing = set(ax.get_children())
    gene_set(ax)
    transform = (Affine2D().translate(-90.5, -96).scale(.72)
                 .translate(*target) + ax.transData)
    for item in ax.get_children():
        if item not in existing:
            item.set_transform(transform)

def figure1():
    source = frame('Cohort framework')
    numbers = source.set_index('Key')['Value'].to_dict()
    if int(numbers['fixed_genes']) != 30:
        raise ValueError('Figure 1 uses 30 schematic gene tiles.')
    for key in numbers:
        if key != 'icgc_hcc_morphology':
            numbers[key] = int(numbers[key])
    TEXTS.clear()
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 8, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'savefig.facecolor': 'white', 'text.color': INK, 'lines.solid_capstyle': 'round',
    })
    fig = plt.figure(figsize=(7.2047244094, 3.1496062992))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 183), ylim=(0, 80))
    ax.axis('off')

    for x, w, title, color in [(6, 51, 'TCGA discovery', BLUE),
                               (72, 39, f"Historical MISP\n({numbers['fixed_genes']} genes)", TEAL),
                               (122, 55, 'External evaluation', PLUM)]:
        label(ax, x, 76 if x==72 else 75.5, title, 8.0 if x==72 else 9.2, 'bold', color)
        ax.plot([x, x+w], [71.5]*2, color=color, lw=1.0)

    # Modality names plus a single count replace the repeated discovery title.
    ax.add_patch(Circle((18.5, 54), 9.0, facecolor='#EDF3F8', edgecolor='none'))
    ax.add_patch(Circle((44, 54), 9.0, facecolor='#EDF5F2', edgecolor='none'))
    dna(ax, 18.5, 54, height=14, width=4.6)
    rna(ax, 44, 54, height=14)
    label(ax, 18.5, 42.0, 'DNA methylation', 7.2, color=BLUE, ha='center')
    label(ax, 44, 42.0, 'RNA', 7.2, color=TEAL, ha='center')
    label(ax, 31.5, 35.8, f"{numbers['tcga_paired']} paired tumours", 8.2, 'bold', ha='center')
    panel(ax, 6, 6, 51, 25, '#F4F7FA')
    label(ax, 10, 25.3, f"CNV/purity complete: ≤{numbers['tcga_cnv_purity']}", 7.6)
    label(ax, 10, 17.9, f"DNA proxies: {numbers['tcga_dna']}", 7.6)
    label(ax, 10, 10.5, f"Joint DNA/ABSOLUTE\nscore analysis: {numbers['tcga_dna_absolute']}", 7.6)

    arrow(ax, (58, 54), (76, 54))
    ax.add_patch(Circle((90.5, 54), 12, facecolor='#EAF4F0', edgecolor='none'))
    gene_icon(ax, (90.5, 54))
    label(ax, 90.5, 37.7, 'Genes and weights fixed\nin primary analyses',
          8.1, 'bold', TEAL, ha='center')

    # Context is secondary: remove the five decorative icons and keep direct labels.
    label(ax, 72, 28.0, 'Supporting analyses', 7.3, 'bold', MUTED)
    ax.plot([72, 111], [25.5]*2, color='#D6DFE4', lw=.55)
    support = [
        'Composition proxies', 'Immune expression', 'Normal-cell references',
        'HCC single-cell localization', 'Post hoc clinical sensitivities',
    ]
    for y, title in zip([21.6, 17.3, 13.0, 8.7, 4.4], support):
        ax.add_patch(Circle((72.1, y), .35, color='#9CABB5'))
        label(ax, 74.2, y, title, 7.2, color=MUTED)

    ax.plot([104, 116], [54, 54], color='#8B9BA5', lw=.9, zorder=1)
    ax.plot([116, 116], [21.8, 54], color='#8B9BA5', lw=.9, zorder=1)
    arrow(ax, (116, 54), (122, 54))
    arrow(ax, (116, 21.8), (122, 21.8))

    panel(ax, 122, 41.5, 55, 27.5, '#F5F2F8')
    ax.plot([122, 122], [44, 66.5], color=PLUM, lw=1.0)
    label(ax, 127, 64.2, f"GEO HCC · {numbers['geo_paired']} exact pairs", 8.5, 'bold', PLUM)
    label(ax, 127, 58.2, 'Methylation + RNA scores', 7.6)
    label(ax, 127, 52.2, 'Fixed-score correlation', 7.6)
    label(ax, 127, 46.2, 'Matched-background comparisons', 7.6)

    panel(ax, 122, 5.5, 55, 32.5, '#F5F2F8')
    ax.plot([122, 122], [8, 35.5], color=PLUM, lw=1.0)
    label(ax, 127, 33.0, f"ICGC RNA · {numbers['icgc_hcc']} HCC donors", 8.5, 'bold', PLUM)
    label(ax, 127, 27.0, f"Morphology {numbers['icgc_hcc_morphology']}", 7.6)
    label(ax, 127, 21.5, f"{numbers['icgc_mixed']} mixed donors: sensitivity", 7.6)
    label(ax, 127, 16.0, f"{numbers['icgc_mixed']}-donor calibration retained", 7.6)
    label(ax, 127, 10.5, 'OS: exploratory', 7.6)

    common.OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(common.OUTPUT / 'Figure_1.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    save(fig, 'Figure_1', dpi=600)


if __name__ == '__main__':
    run(figure1)

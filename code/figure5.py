"""Figure 5. Methylation-score and immune-expression associations."""
from common import frame, letter, revision_style, run, save
import matplotlib.pyplot as plt
import numpy as np

def figure5():
    revision_style()
    data = frame("Methylation and immunity")
    endpoints = ["IFN-γ", "CD8 T-cell score", "NK-cell score", "CYT"]
    adjustments = ["Unadjusted", "DNA broad", "DNA hierarchical"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), layout="constrained")
    seen = {}
    for ax, cohort, panel in zip(axes, ("TCGA", "GEO"), "ab"):
        for j, adjustment in enumerate(adjustments if cohort == 'TCGA' else adjustments[:1]):
            d = data[data.Cohort.eq(cohort) & data.Adjustment.eq(adjustment)].set_index('Endpoint').loc[endpoints]
            x = d['ρ'].to_numpy()
            seen[cohort, adjustment] = x
            y = np.arange(4) + ((j-1)*.17 if cohort == 'TCGA' else 0)
            ax.scatter(x, y, s=24, color=["#2166AC", "#B35806", "#777777"][j],
                       marker=["o", "s", "^"][j], label=["Raw", "DNA broad", "DNA hierarchical"][j])
        ax.axvline(0, color="#AAAAAA", lw=.7)
        ax.set_yticks(range(4), ["IFN-γ", "CD8", "NK", "CYT"])
        ax.set_ylim(3.35, -.35)
        title = cohort
        ax.set(xlim=(-.45, .2), xlabel="Methylation score–immune expression ρ", title=title)
        letter(ax, panel)
    assert list(seen) == [('TCGA', a) for a in adjustments] + [('GEO', 'Unadjusted')]
    np.testing.assert_array_equal(np.round(seen['GEO','Unadjusted'], 3), [-.209, -.309, -.234, -.259])
    np.testing.assert_array_equal(np.round(seen['TCGA','Unadjusted'], 3), [-.336, -.372, -.339, -.324])
    np.testing.assert_array_equal(np.round(seen['TCGA','DNA broad'], 3), [-.169, -.078, -.066, -.011])
    axes[0].legend(loc='upper right', frameon=True, facecolor='white', edgecolor='white', framealpha=1, fontsize=7)
    save(fig, 'Figure_5')

if __name__ == '__main__':
    run(figure5)

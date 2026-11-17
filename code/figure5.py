"""Figure 5."""

from common import frame, save, revision_style, letter, run

import matplotlib.pyplot as plt
import numpy as np


def figure5():
    revision_style()
    data = frame("Methylation and immunity")
    endpoints = ["IFN-γ","CD8 T-cell score","NK-cell score","CYT"]
    adjustments = ["Unadjusted","DNA broad","DNA hierarchical"]
    fig, axes = plt.subplots(1,2,figsize=(7.2,3.6),layout="constrained")
    for ax, cohort, panel in zip(axes,["TCGA","GEO"],"ab"):
        for j, adjustment in enumerate(adjustments):
            d = data[data.Cohort.eq(cohort)&data.Adjustment.eq(adjustment)]
            d = d.set_index("Endpoint").loc[endpoints]
            ax.scatter(d["ρ"],np.arange(4)+(j-1)*.17,s=24,
                       color=["#2166AC","#B35806","#777777"][j],
                       marker=["o","s","^"][j],
                       label=["Raw","DNA broad","DNA hierarchical"][j])
        ax.axvline(0,color="#AAAAAA",lw=.7)
        ax.set_yticks(range(4),["IFN-γ","CD8","NK","CYT"])
        ax.invert_yaxis()
        ax.set(xlim=(-.45,.2),xlabel="Methylation score–immune expression ρ",title=cohort)
        letter(ax,panel)
    axes[1].legend(loc="lower left",frameon=False,fontsize=7)
    save(fig,"Figure_5")


if __name__ == "__main__":
    run(figure5)

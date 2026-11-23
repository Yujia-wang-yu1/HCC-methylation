"""Figure S6."""

from common import frame, save, revision_style, run

import matplotlib as mpl

import matplotlib.pyplot as plt

import numpy as np

def figure_s6():

    revision_style()

    mpl.rcParams.update({"font.size":9,"axes.titlesize":10,"axes.labelsize":9,

                         "xtick.labelsize":8,"ytick.labelsize":8})

    data = frame("Matched sets")

    statistics = frame("Matched-set statistics").set_index("Weighting")

    values=data["ρ"].dropna()

    bins=np.linspace(np.floor(values.min()*10)/10,np.ceil(values.max()*10)/10,31)

    schemes=["Equal weights","Permuted fixed weights","TCGA-derived weights"]

    fig,axes=plt.subplots(1,3,figsize=(7.2,2.7),layout="constrained")

    for ax,scheme in zip(axes,schemes):

        values=data.loc[data.Weighting.eq(scheme),"ρ"].dropna()

        row=statistics.loc[scheme]

        assert len(values)==1000

        ax.hist(values,bins=bins,color="#D4E5EE",edgecolor="white")

        ax.axvline(row["ρ"],color="#D55E00")

        ax.set(title=scheme,xlabel="GEO correlation (ρ)",ylabel="Matched sets")

        ax.text(.03,.97,f'n={len(values):,}\nP={row["Empirical tail P"]:.3f}',

                transform=ax.transAxes,va="top",fontsize=8)

    save(fig,"Figure_S6")

if __name__ == "__main__":

    run(figure_s6)


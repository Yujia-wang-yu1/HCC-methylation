"""Figure S4."""

from common import frame, save, revision_style, run

import matplotlib.pyplot as plt

import numpy as np

def figure_s4():

    revision_style()

    data = frame("Histology")

    choices = [

        ("Methylation–RNA","raw","Methylation–RNA"),

        ("IFN-γ","raw","IFN-γ raw"),("IFN-γ","ESTIMATE-adjusted","IFN-γ adjusted"),

        ("CD8 T-cell score","raw","CD8 raw"),("CD8 T-cell score","ESTIMATE-adjusted","CD8 adjusted"),

        ("NK-cell score","raw","NK raw"),("NK-cell score","ESTIMATE-adjusted","NK adjusted"),

        ("CYT","raw","CYT raw"),("CYT","ESTIMATE-adjusted","CYT adjusted"),

    ]

    cohorts = ["S0: all samples","S1: exclude AAU7","S2: exclude AAU7 and combined histology","S3: HCC NOS"]

    labels = ["S0 n=371","S1 n=370","S2 n=364","S3 n=333"]

    fig, ax = plt.subplots(figsize=(7.2,4.8),layout="constrained")

    for k,(cohort,label) in enumerate(zip(cohorts,labels)):

        values=[]

        for endpoint,adjustment,_ in choices:

            d=data[data.Cohort.eq(cohort)&data.Endpoint.eq(endpoint)&data.Adjustment.eq(adjustment)]

            assert len(d)==1

            values.append(d.iloc[0]["ρ"])

        ax.scatter(values,np.arange(9)+(k-1.5)*.12,label=label,

                   color=["#222222","#777777","#2166AC","#B35806"][k],

                   marker=["o","s","^","D"][k],s=25)

    ax.set_yticks(range(9),[v[2] for v in choices])

    ax.invert_yaxis()

    ax.set(xlim=(-1,.6),xlabel="Spearman or ESTIMATE-adjusted partial Spearman ρ")

    ax.axvline(0,color="#AAAAAA",lw=.7)

    ax.legend(loc="lower left",bbox_to_anchor=(0,1.02),ncol=4,frameon=False,

              fontsize=7.5,handletextpad=.35,columnspacing=1)

    save(fig,"Figure_S4")

if __name__ == "__main__":

    run(figure_s4)


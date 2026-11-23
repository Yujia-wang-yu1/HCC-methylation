"""Figure S7."""

from common import frame, save, revision_style, run

import matplotlib.pyplot as plt

def figure_s7():

    revision_style()

    d=frame("Normal cells").set_index("Gene")

    assert d.shape==(27,8)

    fig,ax=plt.subplots(figsize=(7.2,7.8),layout="constrained")

    im=ax.imshow(d.to_numpy(),vmin=0,vmax=1,aspect="auto",cmap="viridis")

    ax.set_yticks(range(len(d)),d.index)

    ax.set_xticks(range(8),["Hepatocyte","Monocyte","B cell","CD4 T","NK","CD8 T","Neutrophil","Endothelium"],rotation=40,ha="right")

    ax.set(title="Normal-cell reference TSS200 methylation")

    fig.colorbar(im,ax=ax,shrink=.5,label="Mean CpG β")

    save(fig,"Figure_S7")

if __name__ == "__main__":

    run(figure_s7)


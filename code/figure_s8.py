"""Figure S8."""

from common import frame, save, revision_style, run

import matplotlib.pyplot as plt
import numpy as np


def figure_s8():
    revision_style()
    d=frame("Single-cell expression").set_index("Gene")
    assert d.shape==(30,8)
    fig,ax=plt.subplots(figsize=(7.2,8.1),layout="constrained")
    im=ax.imshow(d.to_numpy(),aspect="auto",cmap="magma",vmin=0,vmax=float(np.nanmax(d)))
    ax.set_yticks(range(len(d)),d.index)
    ax.set_xticks(range(8),d.columns,rotation=40,ha="right")
    ax.set(title="HCC single-cell expression by author annotation")
    fig.colorbar(im,ax=ax,shrink=.45,label="Mean log1p(CP10K)")
    save(fig,"Figure_S8")


if __name__ == "__main__":
    run(figure_s8)

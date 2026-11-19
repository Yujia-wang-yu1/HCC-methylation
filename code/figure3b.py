"""Figure 3b."""

from common import frame, save, revision_style, run

import matplotlib as mpl
import matplotlib.pyplot as plt


def figure3b():
    revision_style()
    mpl.rcParams.update({"font.size":9,"axes.titlesize":10,"axes.labelsize":9,
                         "xtick.labelsize":8,"ytick.labelsize":8})
    q = frame("Thresholds")
    fig, ax = plt.subplots(figsize=(3.4,4.3),layout="constrained")
    ax.plot(q["Maximum attenuation (%)"],q["Immune candidates (n)"],
            "o-",color="#0072B2",lw=1.4)
    ax.scatter([50],[30],s=65,marker="s",facecolor="#D55E00",zorder=4)
    ax.annotate("Original rule\n50%: 30 genes",(50,30),xytext=(32,15),
                arrowprops={"arrowstyle":"-","color":"#777777"},fontsize=8)
    ax.set(xlabel="Maximum permitted attenuation (%)",
           ylabel="Immune candidates after all filters",ylim=(0,38),xlim=(-4,104),
           xticks=[0,25,50,75,100])
    ax.set_title("Threshold sensitivity",loc="left")
    ax.text(-.12,1.04,"b",transform=ax.transAxes,fontweight="bold",fontsize=12,va="bottom")
    save(fig,"Figure_3b")


if __name__ == "__main__":
    run(figure3b)

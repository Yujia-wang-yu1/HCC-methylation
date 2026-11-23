"""Figure S9."""
from common import frame, save, revision_style, letter, run
import matplotlib.pyplot as plt
def figure_s9():
    revision_style()
    blue,orange,gray="#2166AC","#B35806","#777777"
    data=frame("Paired scores")
    observed_by_cohort = {}
    fig,axes=plt.subplots(2,2,figsize=(7.2,6),layout="constrained")
    for ax,cohort,panel in zip(axes[0],["TCGA","GEO"],"ab"):
        d=data[data.Cohort.eq(cohort)]
        ax.scatter(d["Methylation score"],d["RNA suppression score"],s=10,alpha=.5,color=blue,edgecolors="none")
        ax.set(xlabel="Methylation score",ylabel="RNA suppression score",title=cohort+" paired tumours")
        letter(ax,panel)
        rho = d["Methylation score"].rank().corr(d["RNA suppression score"].rank())
        observed_by_cohort[cohort] = rho
        text=f"{rho:.3f}; n = {len(d)}"
        ax.text(.03,.96,("Raw Spearman ρ = " if cohort == "TCGA" else "Spearman ρ = ")+text,transform=ax.transAxes,va="top")
    ax=axes[1,0]
    bg=frame("Matched sets")
    values=bg.loc[bg.Weighting.eq("TCGA-derived weights"),"ρ"]
    assert len(values)==1000
    ax.hist(values,bins=28,color="#CCD7E3",edgecolor="white")
    statistic=frame("Matched-set statistics").set_index("Weighting").loc["TCGA-derived weights"]
    if abs(statistic["ρ"]-observed_by_cohort["GEO"]) > 1e-12:
        raise ValueError("Matched-background reference does not agree with the GEO paired scores.")
    ax.axvline(statistic["ρ"],color=orange,lw=1.5)
    ax.set(xlabel="GEO paired-score ρ",ylabel="Matched sets",title="Conditional matched background")
    ax.text(.03,.96,f"Upper-tail P = {statistic['Empirical tail P']:.3f}",transform=ax.transAxes,va="top")
    letter(ax,"c")
    ax=axes[1,1]
    d=frame("DNA-adjusted coupling")
    d=d.loc[d["Cohort"].eq("TCGA")]
    if len(d) != 2 or set(d["n"]) != {354}:
        raise ValueError("Fig. S9d requires the two TCGA estimates on 354 matched cases.")
    positions,labels=[],[]
    for j,(_,row) in enumerate(d.iterrows()):
        y=1-j
        raw=row["Adjustment"]=="Unadjusted, DNA-eligible samples"
        ax.errorbar(row["ρ"],y,xerr=[[row["ρ"]-row["Bootstrap 95% CI lower"]],
                                   [row["Bootstrap 95% CI upper"]-row["ρ"]]],
                    fmt="o",color=blue if raw else orange,capsize=2,ms=4)
        positions.append(y)
        labels.append("Unadjusted" if raw else "DNA adjusted")
    ax.axvline(0,color=gray,lw=.7)
    ax.set_yticks(positions,labels)
    ax.set(xlim=(-.17,.62),xlabel="ρ with bootstrap 95% CI",title="TCGA same-case comparison",ylim=(-.55,1.55))
    letter(ax,"d")
    save(fig,"Figure_S9")
if __name__ == "__main__":
    run(figure_s9)

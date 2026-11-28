"""Exploratory clinical and TMB families using fixed scores; no outcome-based refitting."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from analysis_core import bh,corr

def run(input_dir,output_dir):
    D=Path(input_dir);O=Path(output_dir);O.mkdir(parents=True,exist_ok=True)
    t=pd.read_csv(D/'TCGA_scores_covariates.tsv.gz',sep='\t');grade=pd.to_numeric(t.grade.astype(str).str.extract(r'([1-4])',expand=False),errors='coerce');stage=t.stage_numeric
    clinical=[]
    for score in ['MISP_RNA_suppression_score','MISP_methylation_score']:
        for name,y in [('age',t.age_years),('ordinal_stage_I_IV',stage),('ordinal_grade_1_4',grade)]:clinical.append({'score':score,'clinical_variable':name,'test':'Spearman',**corr(t[score],y)})
        for name,y in [('stage_group_I_IV',stage),('grade_group_1_4',grade),('sex_at_birth',t.male)]:
            ok=t[score].notna()&y.notna();groups=[t.loc[ok&y.eq(v),score].to_numpy() for v in sorted(y[ok].unique())]
            result=stats.mannwhitneyu(*groups,alternative='two-sided') if name=='sex_at_birth' else stats.kruskal(*groups)
            clinical.append({'score':score,'clinical_variable':name,'test':'Mann_Whitney_U' if name=='sex_at_birth' else 'Kruskal_Wallis_H','n':int(ok.sum()),'statistic':float(result.statistic),'p':float(result.pvalue)})
    c=pd.DataFrame(clinical);c['FDR_12_clinical_tests']=bh(c.p);c.to_csv(O/'TCGA_clinical_associations.tsv',sep='\t',index=False)
    # The MAF-derived TMB proxy uses the same accepted GDC flags as the CTNNB1
    # comparison. Counts were filtered before deduplication; absent MAF cases
    # remain missing rather than being assigned zero variants. The fixed 38-Mb
    # denominator is a proxy convention, not an individual callable-exome size.
    tc=pd.read_csv(D/'TCGA_TMB_variant_filter_counts.tsv',sep='\t').set_index('patient_id').reindex(t.patient_id)
    expected_count=tc.accepted_unique_nonsyn_count.to_numpy()
    if not np.allclose(t.TMB_nonsyn_variant_count,expected_count,equal_nan=True):
        raise ValueError('TMB input does not match accepted-variant count metadata')
    if not np.allclose(t.TMB_nonsyn_per_Mb,expected_count/38.0,equal_nan=True):
        raise ValueError('TMB proxy input does not match the fixed 38-Mb convention')
    rows=[]
    for score in ['MISP_RNA_suppression_score','MISP_methylation_score']:
        for adjustment,cov in [('raw',None),('ESTIMATE_rank',t.ESTIMATE_TumorPurity.to_numpy())]:rows.append({'score':score,'adjustment':adjustment,**corr(t[score],t.TMB_nonsyn_per_Mb,cov)})
    b=pd.DataFrame(rows);b['FDR_four_TMB_tests']=bh(b.p);b.to_csv(O/'TCGA_TMB_associations.tsv',sep='\t',index=False)
    subsets={'S0_all':np.ones(len(t),bool),'S1_exclude_AAU7':t.patient_id.ne('TCGA-K7-AAU7').to_numpy(),'S2_exclude_AAU7_combined':(t.patient_id.ne('TCGA-K7-AAU7')&t.primary_diagnosis.ne('Combined hepatocellular carcinoma and cholangiocarcinoma')).to_numpy(),'S3_HCC_NOS':t.primary_diagnosis.eq('Hepatocellular carcinoma, NOS').to_numpy()}
    hist=[]
    for cohort,mask in subsets.items():
        for ep in ['IFN_gamma_response_score','CD8_T_cell_score','NK_cell_score','cytolytic_activity','ESTIMATE_ImmuneScore','immune_checkpoint_composite','TMB_nonsyn_per_Mb']:
            for adjustment,cv in [('raw',None),('ESTIMATE_rank',t.loc[mask,'ESTIMATE_TumorPurity'].to_numpy())]:hist.append({'cohort':cohort,'endpoint':ep,'adjustment':adjustment,**corr(t.loc[mask,'MISP_RNA_suppression_score'],t.loc[mask,ep],cv)})
        hist.append({'cohort':cohort,'endpoint':'paired_methylation_RNA','adjustment':'raw',**corr(t.loc[mask,'MISP_methylation_score'],t.loc[mask,'MISP_RNA_suppression_score'])})
        for label,y in [('age',t.age_years),('stage_I_IV',stage),('grade_1_4',grade)]:hist.append({'cohort':cohort,'endpoint':label,'adjustment':'raw',**corr(t.loc[mask,'MISP_RNA_suppression_score'],y[mask])})
        ok=mask&t.male.notna().to_numpy();groups=[t.loc[ok&t.male.eq(v),'MISP_RNA_suppression_score'].to_numpy() for v in [0,1]];u=stats.mannwhitneyu(*groups,alternative='two-sided');hist.append({'cohort':cohort,'endpoint':'sex_at_birth','adjustment':'raw','n':int(ok.sum()),'statistic':float(u.statistic),'p':float(u.pvalue)})
    h=pd.DataFrame(hist);h['FDR_19_tests_within_subset']=h.groupby('cohort').p.transform(bh);h.to_csv(O/'TCGA_histology_sensitivity.tsv',sep='\t',index=False)
    return c,b

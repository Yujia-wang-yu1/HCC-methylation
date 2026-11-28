"""Post hoc, fully specified score-level covariate sensitivities; v4 inputs remain read-only."""
from pathlib import Path
import sys,json,hashlib,os
import numpy as np,pandas as pd
from scipy import stats
R = Path(__file__).resolve().parents[1]
V4 = R.parent / 'MISP_Code'
D = Path(os.environ.get('MISP_REPRO_INPUTS', R.parent)) / 'MISP_Code/data'
A = Path(os.environ.get('MISP_REPRO_OUTPUT', R.parent)) / 'clinical_covariates'
A.mkdir(parents=True, exist_ok=True)
sys.path.insert(0,str(V4));from analysis_core import bh
SEED=20261114;REPS=2000
models={
 'DNA_broad':{'continuous':['IC','Fib'],'binary':[],'categorical':[]},
 'DNA_ABSOLUTE':{'continuous':['IC','Fib','ABSOLUTE_purity'],'binary':[],'categorical':[]},
 'DNA_ABSOLUTE_age_sex':{'continuous':['IC','Fib','ABSOLUTE_purity','age_years'],'binary':['male'],'categorical':[]},
 'DNA_ABSOLUTE_age_sex_stage_grade':{'continuous':['IC','Fib','ABSOLUTE_purity','age_years'],'binary':['male'],'categorical':['stage_numeric','grade_numeric']}
}
spec={'objective':'Condition the original fixed molecular scores on measured age/sex/stage/grade and composition/purity; do not refit genes/weights or construct a clinical predictor.', 'models':models,'stage_grade':'Categorical I–IV/1–4, dummy coded after selecting complete cases; no imputation.','rank':'Spearman ranks for continuous scores/covariates; binary/dummy covariates remain 0/1, equivalent residual spaces to their average ranks.','matching':'Each raw/adjusted comparison uses exactly the same model-specific complete-case patients.','families':'Four methylation–immune outcomes per model, and four RNA–immune outcomes per model separately; paired-score sensitivities corrected across four model specifications.','intervals':'2000 patient bootstrap replicates with ranks recomputed for all nine full-clinical-model associations; percentile95%. Intervals condition on original scores, weights and cohort calibration.','seed':SEED,'repetitions':REPS,'no_predictive_claim':'No outcome-driven refitting; no cross-validated predictive performance or clinical treatment endpoint.'}
(A/'analysis_specification.json').write_text(json.dumps(spec,indent=2),encoding='utf8')
t=pd.read_csv(D/'TCGA_scores_covariates.tsv.gz',sep='\t').set_index('sample_id');t['grade_numeric']=pd.to_numeric(t.grade.astype(str).str.extract(r'([1-4])',expand=False),errors='coerce')
assert len(t)==371 and t.index.is_unique
endpoints=['IFN_gamma_response_score','CD8_T_cell_score','NK_cell_score','cytolytic_activity']
tests=[('paired_scores','paired_methylation_RNA','MISP_methylation_score','MISP_RNA_suppression_score')]
tests += [(kind,ep,score,ep) for kind,score in [('methylation_immune','MISP_methylation_score'),('RNA_immune','MISP_RNA_suppression_score')] for ep in endpoints]
def design(frame,spec):
 parts=[];names=[]
 for c in spec['continuous']:
  parts.append(stats.rankdata(frame[c]));names.append(c)
 for c in spec['binary']:parts.append(frame[c].to_numpy(float));names.append(c)
 for c in spec['categorical']:
  levels=sorted(frame[c].unique());
  for level in levels[1:]:parts.append(frame[c].eq(level).to_numpy(float));names.append(f'{c}_{int(level)}')
 return np.column_stack([np.ones(len(frame)),*parts]),names
def calculate(x,y,C):
 rx=stats.rankdata(x);ry=stats.rankdata(y);rank=np.linalg.matrix_rank(C)
 ex=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];ey=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0]
 rho=float(np.corrcoef(ex,ey)[0,1]);df=len(x)-rank-1
 p=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-16,1-rho*rho)),df))
 return rho,p,df,rank
rows=[];availability=[];cases=[];checks=[]
for label,spec0 in models.items():
 needed=spec0['continuous']+spec0['binary']+spec0['categorical']
 eligible=t[needed].notna().all(axis=1)
 availability.append({'model':label,'eligible_covariates':int(eligible.sum()),'missing_covariates':int((~eligible).sum()),'stage_counts':json.dumps(t.loc[eligible,'stage_numeric'].value_counts().sort_index().to_dict()),'grade_counts':json.dumps(t.loc[eligible,'grade_numeric'].value_counts().sort_index().to_dict())})
 for patient in t.loc[eligible,'patient_id']:cases.append({'model':label,'patient_id':patient})
 for family,endpoint,xcol,ycol in tests:
  q=t.loc[eligible&np.isfinite(t[xcol])&np.isfinite(t[ycol])];C,names=design(q,spec0);assert np.linalg.matrix_rank(C)==C.shape[1],(label,names)
  x=q[xcol].to_numpy();y=q[ycol].to_numpy();rho,p,df,rank=calculate(x,y,C);raw=stats.spearmanr(x,y)
  # Exact invariance under centering/scaling of nuisance covariates verifies the numerical rank-residual implementation.
  Cs=C.copy();Cs[:,1:]=(Cs[:,1:]-Cs[:,1:].mean(0))/Cs[:,1:].std(0);check_rho=calculate(x,y,Cs)[0];assert abs(check_rho-rho)<1e-10
  # Project both ranked endpoints together as an independent algebraic residual check.
  Q=np.linalg.qr(Cs,mode='reduced')[0];rr=np.column_stack([stats.rankdata(x),stats.rankdata(y)]);res=rr-Q@(Q.T@rr);qr_rho=float(np.corrcoef(res.T)[0,1]);assert abs(qr_rho-rho)<1e-10
  out={'model':label,'family':family,'endpoint':endpoint,'n':len(q),'covariate_rank':rank-1,'df':df,'raw_same_case_rho':float(raw.statistic),'raw_same_case_p':float(raw.pvalue),'adjusted_rho':rho,'p':p,'covariates':';'.join(names),'bootstrap_reps':0}
  if label=='DNA_ABSOLUTE_age_sex_stage_grade':
   rng=np.random.default_rng(SEED+len(rows));boot=[];deficient=0
   for _ in range(REPS):
    ix=rng.integers(0,len(q),len(q));qb=q.iloc[ix];Cb,_=design(qb,spec0);deficient+=np.linalg.matrix_rank(Cb)<Cb.shape[1]
    val=calculate(qb[xcol].to_numpy(),qb[ycol].to_numpy(),Cb)[0]
    if np.isfinite(val):boot.append(val)
   lo,hi=np.quantile(boot,[.025,.975]);out.update(CI_lower=float(lo),CI_upper=float(hi),bootstrap_reps=len(boot),bootstrap_rank_deficient_reps=deficient)
  rows.append(out);checks.append({'model':label,'endpoint':endpoint,'family':family,'scaled_rho_error':abs(check_rho-rho),'QR_rho_error':abs(qr_rho-rho)})
  print(label,family,endpoint,'n',len(q),'rho',round(rho,6),flush=True)
result=pd.DataFrame(rows);result['BH_adjusted_P']=np.nan
for (model,family),ix in result[result.family.ne('paired_scores')].groupby(['model','family']).groups.items():result.loc[ix,'BH_adjusted_P']=bh(result.loc[ix,'p'])
ix=result.family.eq('paired_scores');result.loc[ix,'BH_adjusted_P']=bh(result.loc[ix,'p'])
result.to_csv(A/'fixed_MISP_fully_adjusted_score_associations.tsv',sep='\t',index=False)
pd.DataFrame(availability).to_csv(A/'covariate_case_availability.tsv',sep='\t',index=False);pd.DataFrame(cases).to_csv(A/'eligible_patient_registry.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(A/'independent_algebra_checks.tsv',sep='\t',index=False)
(A/'input_provenance.json').write_text(json.dumps({'input_file':'TCGA_scores_covariates.tsv.gz','source_version':V4.name,'sha256':hashlib.sha256((D/'TCGA_scores_covariates.tsv.gz').read_bytes()).hexdigest(),'371_original_fixed_score_calibration':True,'no_genes_or_weights_refitted':True},indent=2),encoding='utf8')
print('COMPLETE',len(result),'association results',flush=True)

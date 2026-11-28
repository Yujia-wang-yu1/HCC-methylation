"""Run from any directory. Uses only explicit packaged inputs; no network or installation."""
from pathlib import Path
import argparse,json,subprocess,sys,os
import numpy as np,pandas as pd
from scipy import stats
from analysis_core import bh,corr,ols,zrows,frozen_score

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input-dir',type=Path,default=Path(__file__).resolve().parent/'data');ap.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'results');ap.add_argument('--stage',choices=['scores','genes','background','all'],default='scores');ap.add_argument('--rscript',default=None);args=ap.parse_args();D=args.input_dir.resolve();O=args.output_dir.resolve();O.mkdir(parents=True,exist_ok=True)
    def read(n,index=None):return pd.read_csv(D/n,sep='\t',index_col=index)
    def save(x,n):x.to_csv(O/n,sep='\t',index=False)
    t=read('TCGA_scores_covariates.tsv.gz',0);g=read('GEO_scores.tsv.gz',0);ic=read('ICGC_scores.tsv.gz',0);hcc=ic.tumour_histological_type.eq('8170/3').to_numpy();weights=read('fixed_MISP_weights.tsv',0).MISP_weight.sort_index();genes=list(weights.index);w=weights.to_numpy();w=w/w.sum();tx=np.load(D/'TCGA_RNA.npz');td=np.load(D/'TCGA_DNA.npz');gx=np.load(D/'GEO_standardized.npz');ix=np.load(D/'ICGC_standardized.npz');tg=list(tx['genes']);cg=list(gx['genes']);tidx={v:i for i,v in enumerate(tg)};cidx={v:i for i,v in enumerate(cg)}
    assert list(tx['samples'])==list(t.index) and int(hcc.sum())==203
    xr=pd.DataFrame(tx['RNA'].astype(float),index=tg,columns=t.index);xm=pd.DataFrame(td['methylation'].astype(float),index=tg,columns=t.index);cnv_array=td['CNV'].astype(float)
    tcga_z=zrows(xr).to_numpy();arrays={'TCGA':tcga_z,'ICGC':ix['RNA_z'].astype(float),'GEO_methylation':gx['methylation_z'].astype(float),'GEO_RNA':gx['RNA_z'].astype(float)}
    def weighted(key,gs,v):
        index=tidx if key=='TCGA' else cidx
        values = arrays[key][[index[x] for x in gs]]
        scores = np.nansum(values * np.asarray(v)[:, None], axis=0)
        scores[~np.isfinite(values).any(axis=0)] = np.nan
        return scores
    checks={}
    checks['TCGA_methylation_max_abs_error']=float(np.max(np.abs(frozen_score(xm,weights).to_numpy()-t.MISP_methylation_score)))
    checks['TCGA_RNA_max_abs_error']=float(np.max(np.abs(-weighted('TCGA',genes,w)-t.MISP_RNA_suppression_score)))
    checks['GEO_methylation_max_abs_error']=float(np.max(np.abs(weighted('GEO_methylation',genes,w)-g.MISP_GEO_methylation_suppression_score)))
    checks['GEO_RNA_max_abs_error']=float(np.max(np.abs(-weighted('GEO_RNA',genes,w)-g.MISP_GEO_RNA_suppression_score)))
    checks['ICGC_RNA_max_abs_error']=float(np.max(np.abs(-weighted('ICGC',genes,w)-ic.MISP_RNA_suppression_score)))
    assert max(checks.values())<5e-6,checks
    defs=json.loads((D/'immune_definitions.json').read_text());eps=list(defs);ex=np.load(D/'immune_and_Inflamed_expression.npz');expression={co:pd.DataFrame(ex[co+'_expression'],index=ex[co+'_genes']) for co in ['TCGA','GEO','ICGC']}
    def endpoint(co,ep,exclude=()):
        gs=[v for v in defs[ep] if v in expression[co].index and v not in exclude];x=expression[co].loc[gs]
        return (x.mean(axis=0) if ep=='cytolytic_activity' else zrows(x).mean(axis=0)).to_numpy()
    rows=[]
    for co,tab,met,rna in [('TCGA',t,t.MISP_methylation_score,t.MISP_RNA_suppression_score),('GEO',g,g.MISP_GEO_methylation_suppression_score,g.MISP_GEO_RNA_suppression_score)]:
        adjustments = {'raw': np.empty((len(tab), 0))}
        if co == 'TCGA':
            dna = read('TCGA_DNA_composition_proxies.tsv', 0).reindex(t.index)
            adjustments['DNA_broad'] = dna[['IC', 'Fib']].to_numpy()
            adjustments['DNA_ABSOLUTE'] = t[['IC', 'Fib', 'ABSOLUTE_purity']].to_numpy()
            adjustments['DNA_ABSOLUTE_age_sex'] = t[['IC', 'Fib', 'ABSOLUTE_purity', 'age_years', 'male']].to_numpy()
        for model,cv in adjustments.items():
            rows.append({'cohort':co,'model':model,'endpoint':'paired_methylation_RNA',**corr(met,rna,cv)})
            for ep in eps:rows.append({'cohort':co,'model':model,'endpoint':ep,**corr(met,endpoint(co,ep),cv),**{'OLS_'+k:v for k,v in ols(endpoint(co,ep),met,cv).items()}})
    associations=pd.DataFrame(rows);associations['FDR_four_immune']=np.nan
    for (co,mo),z in associations.groupby(['cohort','model']):
        ids=z.index[z.endpoint.ne('paired_methylation_RNA')];associations.loc[ids,'FDR_four_immune']=bh(associations.loc[ids,'p'])
    save(associations,'reproduced_score_associations.tsv')
    if args.stage in ('scores', 'all'):
        from reported_score_results import run as reported_score_run
        reported_score_run(D, O, t, g, ic, hcc, expression, defs, genes)
    from clinical_endpoints import run as clinical_run
    clinical_run(D,O)
    if args.stage in ['genes','all']:
        models={'CNV_DNA_broad':['IC','Fib'],'CNV_DNA_broad_age_sex':['IC','Fib','age_years','male'],'CNV_DNA_ABSOLUTE':['IC','Fib','ABSOLUTE_purity'],'CNV_DNA_ABSOLUTE_age_sex':['IC','Fib','ABSOLUTE_purity','age_years','male']};rows=[]
        for j,gene in enumerate(tg):
            for model,cols in models.items():
                cv=np.column_stack([cnv_array[j],t[cols].to_numpy()])
                try:fit=ols(xr.iloc[j],xm.iloc[j],cv)
                except (ValueError,np.linalg.LinAlgError):fit={'n':int((np.isfinite(cv).all(1)&np.isfinite(xr.iloc[j])&np.isfinite(xm.iloc[j])).sum()),'beta':np.nan,'p':np.nan}
                rows.append({'gene':gene,'model':model,**fit})
            if (j+1)%3000==0:print('Gene models',j+1,flush=True)
        result=pd.DataFrame(rows);result['FDR_genomewide']=result.groupby('model').p.transform(bh);save(result,'reproduced_genomewide_models.tsv.gz')
        summary=[];audit=read('gene_selection_audit.tsv',0);selected=[]
        for model,z in result.groupby('model'):
            f=z[z.gene.isin(genes)];summary.append({'model':model,'genomewide_tests':int(z.p.notna().sum()),'frozen_negative':int(f.beta.lt(0).sum()),'frozen_FDR05':int((f.beta.lt(0)&f.FDR_genomewide.lt(.05)).sum()),'n_min':int(f.n.min()),'n_max':int(f.n.max())})
            az=audit.join(z.set_index('gene')[['beta','FDR_genomewide','std_beta']]);atten=100*(1-az.beta.abs()/az.adjusted_methylation_coefficient.abs());keep=az.genomewide_evidence_tier.eq('Stable')&az.beta.lt(0)&az.FDR_genomewide.lt(.05)&atten.le(50)&az.high_purity_methylation_beta.lt(0)&az.leading_edge_eligible&az.promoter_direction_supported&az.RAB27A_excluded&az.ICGC_coverage&az.GPL13667_coverage
            gs=az.index[keep].tolist();vw=az.loc[gs,'std_beta'].abs().to_numpy();vw=vw/vw.sum()
            for gene,weight in zip(gs,vw):selected.append({'model':model,'gene':gene,'weight':weight,'in_original_MISP':gene in genes})
        save(pd.DataFrame(summary),'reproduced_gene_model_summary.tsv');save(pd.DataFrame(selected),'reproduced_criterion_sensitivity_sets.tsv')
        # Evaluate post hoc selected sets without replacing primary fixed 30-gene scores.
        rows=[]
        for model,z in pd.DataFrame(selected).groupby('model',sort=False):
            gs=z.gene.tolist();v=z.weight.to_numpy();ms=weighted('GEO_methylation',gs,v);rs=-weighted('GEO_RNA',gs,v)
            rows.append({'selection_model':model,'cohort':'GEO','endpoint':'paired_methylation_RNA','adjustment':'raw',**corr(ms,rs)})
            for ep in eps:
                rows.append({'selection_model':model,'cohort':'GEO','endpoint':ep,'adjustment':'raw',**corr(ms,endpoint('GEO',ep,gs))})
            for co in ['TCGA','ICGC']:
                sr=-weighted(co,gs,v);mask=np.ones(len(t),bool) if co=='TCGA' else hcc
                for ep in eps:
                    y=endpoint(co,ep,gs);rows.append({'selection_model':model,'cohort':co if co=='TCGA' else 'ICGC_HCC_203','endpoint':ep,'adjustment':'raw_nonoverlap',**corr(sr[mask],y[mask])})
        se=pd.DataFrame(rows);se['FDR_four_immune']=np.nan
        for _,z in se.groupby(['selection_model','cohort','adjustment']):
            ids=z.index[z.endpoint.ne('paired_methylation_RNA')];se.loc[ids,'FDR_four_immune']=bh(se.loc[ids,'p'])
        save(se,'reproduced_selected_set_evaluations.tsv')
    if args.stage in ['background','all']:
        from matched_registry import regenerate
        sets=regenerate(D,genes);save(sets,'reproduced_matched_registry.tsv');coefs=read('TCGA_background_weight_coefficients.tsv',0);rng=np.random.default_rng(20261113);obs=[];random=[]
        def effects(gs,v):
            rs=-weighted('TCGA',gs,v);ri=-weighted('ICGC',gs,v);result=[{'endpoint':'GEO_methylation_RNA','cohort':'GEO',**corr(weighted('GEO_methylation',gs,v),-weighted('GEO_RNA',gs,v))}]
            for ep in eps:
                yt=endpoint('TCGA',ep,gs);yi=endpoint('ICGC',ep,gs)
                for co,x,y,cv,rank_cv in [('TCGA_raw',rs,yt,None,True),('TCGA_ESTIMATE',rs,yt,t.ESTIMATE_TumorPurity.to_numpy(),True),('TCGA_ABSOLUTE',rs,yt,t.ABSOLUTE_purity.to_numpy(),False),('ICGC_HCC_203',ri[hcc],yi[hcc],None,True)]:result.append({'endpoint':ep,'cohort':co,**corr(x,y,cv,rank_covariates=rank_cv)})
            return result
        for scheme in ['equal','permuted_fixed','refit_TCGA']:
            ow=np.ones(30)/30 if scheme=='equal' else w
            obs.extend([{'scheme':scheme,**q} for q in effects(genes,ow)])
        for row in sets.itertuples(index=False,name=None):
            sid=row[0];gs=list(row[1:])
            for scheme in ['equal','permuted_fixed','refit_TCGA']:
                if scheme=='equal':v=np.ones(30)/30
                elif scheme=='permuted_fixed':v=rng.permutation(w)
                else:v=coefs.reindex(gs).purity_methylation_standardized_beta.abs().to_numpy();v=v/v.sum()
                random.extend([{'set_id':sid,'scheme':scheme,**q} for q in effects(gs,v)])
            if int(sid)%200==0:print('Matched sets',sid,flush=True)
        rr=pd.DataFrame(random);oo=pd.DataFrame(obs);save(rr,'reproduced_background_effects.tsv.gz');save(oo,'reproduced_background_observed.tsv');summary=[]
        for o in oo.itertuples(index=False):
            v=rr.loc[rr.scheme.eq(o.scheme)&rr.endpoint.eq(o.endpoint)&rr.cohort.eq(o.cohort),'rho'].dropna().to_numpy();upper=o.endpoint=='GEO_methylation_RNA';ext=int((v>=o.rho).sum() if upper else (v<=o.rho).sum());lo,hi=np.quantile(v,[.025,.975]);summary.append({'scheme':o.scheme,'endpoint':o.endpoint,'cohort':o.cohort,'observed_rho':o.rho,'n':o.n,'finite_random_n':len(v),'random_q025':lo,'random_q975':hi,'empirical_tail_p':(ext+1)/(len(v)+1)})
        ss=pd.DataFrame(summary);ss['FDR_within_scheme_cohort']=ss.groupby(['scheme','cohort']).empirical_tail_p.transform(bh);save(ss,'reproduced_background_summary.tsv')
    if args.rscript:
        code_dir=Path(__file__).resolve().parent
        for rfile in ['reproduce_DNA_proxies.R','icgc_survival.R']:
            subprocess.run([args.rscript, rfile, str(D), str(O)], cwd=code_dir, check=True)
    (O/'scoring_reproduction_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf8');print(json.dumps(checks,indent=2),flush=True)

if __name__=='__main__':main()

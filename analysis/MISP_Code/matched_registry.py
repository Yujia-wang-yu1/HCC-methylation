"""Recreate the frozen technical-matched registry from its published input features."""
from pathlib import Path
import numpy as np,pandas as pd

def regenerate(input_dir,genes):
    D=Path(input_dir);meta=pd.read_csv(D/'random_background_matching_features.tsv',sep='\t',index_col='gene');coverage=pd.read_csv(D/'random_background_technical_coverage.tsv',sep='\t',index_col='gene')
    keep=coverage.eligible_complete_background.reindex(meta.index).fillna(False).astype(bool)
    pool=meta.loc[keep].copy();pool['cpg_proxy']=pool.cpg_proxy.fillna(pool.cpg_proxy.median())
    features=['expr_mean','expr_sd','probe_n','cpg_proxy','gene_length','missingness']
    med=pool[features].median();scale=(pool[features]-med).abs().median().replace(0,1);xx=pool[features].to_numpy(float);nearest={};names=pool.index.to_numpy();rng=np.random.default_rng(20261114)
    for g in genes:
        target=meta.loc[g];dist=np.abs((xx-target[features].fillna(med).to_numpy(float))/scale.to_numpy(float)).sum(axis=1);same=(pool.cnv_available.to_numpy()==bool(target.cnv_available))&(pool.immune_gene.to_numpy()==bool(target.immune_gene));ix=np.flatnonzero(same);assert len(ix)>=100;nearest[g]=ix[np.argsort(dist[ix],kind='stable')]
    sets=[]
    for _ in range(1000):
        used=set();row=[]
        for g in genes:
            candidates=[]
            for j in nearest[g]:
                if int(j) not in used:candidates.append(j)
                if len(candidates)==100:break
            pick=int(rng.choice(candidates));used.add(pick);row.append(str(names[pick]))
        sets.append(row)
    result=pd.DataFrame({'set_id':np.arange(1,1001),**{f'gene_{j+1}':[row[j] for row in sets] for j in range(30)}})
    expected=pd.read_csv(D/'coverage_matched_random_sets.tsv',sep='\t')
    if not result.equals(expected):raise AssertionError('Regenerated matching registry differs from the frozen reference')
    return result

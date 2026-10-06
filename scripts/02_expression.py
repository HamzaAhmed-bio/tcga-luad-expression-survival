"""Create verified matrices, fit paired PyDESeq2, and freeze candidates."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='1'
import hashlib,json,pickle
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT=Path(__file__).resolve().parents[1]; WORK=Path(os.environ.get('TCGA_WORKDIR',str(ROOT.parents[1]/'work/tcga-luad')))
OUT=ROOT/'results'; FIG=OUT/'figures'
def save(name):
    plt.tight_layout(); plt.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight'); plt.savefig(FIG/(name+'.svg'),bbox_inches='tight'); plt.close()
def main():
    OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
    manifest=pd.read_csv(ROOT/'metadata/download_manifest.csv').sort_values(['patient','sample_type'])
    assert manifest.file_id.is_unique and manifest['sample'].is_unique
    cache=WORK/'matrices.pkl'
    input_hash=hashlib.sha256((ROOT/'metadata/download_manifest.csv').read_bytes()).hexdigest()
    cache_stamp=WORK/'matrix_manifest.sha256'
    if cache.exists() and cache_stamp.exists() and cache_stamp.read_text()==input_hash:
        with cache.open('rb') as f: counts,tpm,annotation=pickle.load(f)
    else:
        cs={}; ts={}; annotation=None; qc=[]; headers=set()
        for i,r in enumerate(manifest.itertuples(),1):
            path=WORK/'counts'/(r.file_id+'.tsv'); data=path.read_bytes()
            assert hashlib.md5(data).hexdigest()==r.md5sum
            headers.add(data.splitlines()[0].decode())
            x=pd.read_csv(path,sep='\t',comment='#')
            special=x[~x.gene_id.str.startswith('ENSG')]
            x=x[x.gene_id.str.startswith('ENSG')].set_index('gene_id')
            assert x.index.is_unique and np.isfinite(x.unstranded).all() and (x.unstranded>=0).all() and (x.unstranded%1==0).all()
            if annotation is None: annotation=x[['gene_name','gene_type']]
            else: assert annotation.equals(x[['gene_name','gene_type']])
            cs[r.sample]=x.unstranded.astype('int64'); ts[r.sample]=x.tpm_unstranded
            q=dict(sample=r.sample,patient=r.patient,condition=r.sample_type,library_size=int(x.unstranded.sum()),detected_genes=int((x.unstranded>0).sum()),tpm_sum=float(x.tpm_unstranded.sum()))
            q.update({row.gene_id:int(row.unstranded) for row in special.itertuples()}); qc.append(q)
            if i%100==0: print(f'Parsed {i} verified files',flush=True)
        counts=pd.DataFrame(cs).T; tpm=pd.DataFrame(ts).T
        pd.DataFrame(qc).to_csv(OUT/'sample_qc.csv',index=False)
        (OUT/'reference_headers.json').write_text(json.dumps(sorted(headers)))
        with cache.open('wb') as f: pickle.dump((counts,tpm,annotation),f)
        cache_stamp.write_text(input_hash)
        annotation.to_csv(OUT/'gene_annotation.csv')
    paired=manifest[manifest.paired_analysis].set_index('sample')
    paired['condition']=np.where(paired.sample_type=='Primary Tumor','Tumor','Normal')
    paired['condition']=pd.Categorical(paired.condition,categories=['Normal','Tumor'])
    assert paired.groupby('patient').size().eq(2).all()
    c=counts.loc[paired.index]; keep=(c>=10).sum(axis=0)>=len(paired)//2
    annotation.assign(passed_expression_filter=keep).to_csv(OUT/'gene_filter_audit.csv')
    c=c.loc[:,keep]; assert c.index.equals(paired.index)
    paired.to_csv(OUT/'paired_metadata.csv'); c.to_csv(OUT/'paired_raw_counts.csv.gz',compression='gzip')
    print(f'Fitting {c.shape[1]} genes across {c.shape[0]} paired samples',flush=True)
    dds=DeseqDataSet(counts=c,metadata=paired[['patient','condition']],design='~patient + condition',refit_cooks=False,n_cpus=4)
    dm=np.asarray(dds.obsm['design_matrix']); assert np.linalg.matrix_rank(dm)==dm.shape[1]
    model_hash=hashlib.sha256((input_hash+'|pydeseq2-0.5.4|patient+condition|counts10-in50|no-refit').encode()).hexdigest()
    model_stamp=WORK/'model_inputs.sha256'
    if not ((WORK/'dds.pkl').exists() and model_stamp.exists() and model_stamp.read_text()==model_hash):
        dds.deseq2()
        with (WORK/'dds.pkl').open('wb') as f: pickle.dump(dds,f)
        model_stamp.write_text(model_hash)
    else:
        with (WORK/'dds.pkl').open('rb') as f: dds=pickle.load(f)
    stats=DeseqStats(dds,contrast=['condition','Tumor','Normal'],alpha=0.05,n_cpus=4,cooks_filter=True,independent_filter=True)
    stats.summary(); result=annotation.join(stats.results_df)
    convergence=dds.var[[k for k in dds.var.columns if 'converged' in k]].fillna(False).all(axis=1)
    result['fit_converged']=convergence.reindex(result.index,fill_value=False)
    result['padj_before_convergence_gate']=result.padj
    result.loc[~result.fit_converged,'padj']=np.nan
    result['passed_expression_filter']=keep
    result['log2FC_CI_low']=result.log2FoldChange-1.96*result.lfcSE
    result['log2FC_CI_high']=result.log2FoldChange+1.96*result.lfcSE
    result.to_csv(OUT/'differential_expression_all_genes.csv.gz',compression='gzip')
    unique=~result.gene_name.duplicated(keep=False)
    candidates=result[(result.gene_type=='protein_coding') & unique & result.gene_name.notna() & (result.baseMean>=100)&(result.padj<.05)&(result.log2FoldChange.abs()>=1)].copy()
    candidates['absolute_stat']=candidates.stat.abs()
    candidates=candidates.reset_index().sort_values(['absolute_stat','gene_id'],ascending=[False,True]).head(10)
    candidates.to_csv(OUT/'frozen_candidates.csv',index=False)
    (OUT/'candidate_freeze.json').write_text(json.dumps(dict(sha256=hashlib.sha256((OUT/'frozen_candidates.csv').read_bytes()).hexdigest(),rule='v2 expression-only rule; frozen before survival and external tests'),indent=2))
    sf=dds.obs['size_factors'].to_numpy()
    norm=c.div(sf,axis=0); log=np.log2(norm+1)
    pd.DataFrame({'sample':c.index,'size_factor':sf}).to_csv(OUT/'normalization_factors.csv',index=False)
    top=log.var().nlargest(1000).index; pca=PCA(n_components=2); coords=pca.fit_transform(log[top])
    plt.figure(figsize=(8,6))
    for label,color in [('Normal','#158a8a'),('Tumor','#cb5264')]:
        ix=paired.condition==label; plt.scatter(coords[ix,0],coords[ix,1],label=label,c=color,s=38,alpha=.8)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})');plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})');plt.title('Matched samples: expression overview\nTop 1,000 variable genes; log2(normalized counts + 1)');plt.legend();save('01_pca')
    rr=result.dropna(subset=['padj','log2FoldChange']);colors=np.where(rr.padj<.05,np.where(rr.log2FoldChange>0,'#cb5264','#158a8a'),'#bbc2cc')
    plt.figure(figsize=(8,6));plt.scatter(rr.log2FoldChange,-np.log10(rr.padj.clip(lower=1e-300)),c=colors,s=5,alpha=.55);plt.axhline(-np.log10(.05),c='grey',ls='--');plt.xlabel('log2 fold change: tumor / normal');plt.ylabel('-log10 adjusted p-value');plt.title('Paired tumor versus non-tumor expression');save('02_volcano')
    hm=log[candidates.gene_id]; hm=(hm-hm.mean())/hm.std();order=paired.sort_values(['condition','patient']).index
    plt.figure(figsize=(12,5));plt.imshow(hm.loc[order].T,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5);plt.yticks(range(len(candidates)),candidates.gene_name);plt.xticks([24.5,74.5],['50 normal samples','50 tumors']);plt.axvline(49.5,c='black',lw=1);plt.colorbar(label='Within-gene z score');plt.title('Ten expression-selected candidates');save('03_heatmap')
    plt.figure(figsize=(8,5));plt.scatter(np.arange(len(c)),c.sum(axis=1)/1e6,c=np.where(paired.condition=='Tumor','#cb5264','#158a8a'));plt.ylabel('Gene-assigned reads (millions)');plt.xlabel('Paired sample index');plt.title('Library sizes: all matched samples retained');save('04_library_sizes')
    diagnostics={key:dds.var[key].value_counts(dropna=False).astype(int).to_dict() for key in dds.var.columns if 'converged' in key}
    # Boolean keys must be serialized as strings for stable JSON.
    diagnostics={k:{str(a):int(b) for a,b in v.items()} for k,v in diagnostics.items()}
    summary=dict(pairs=len(paired)//2,genes_total=len(result),genes_filtered_in=int(keep.sum()),genes_with_padj=int(result.padj.notna().sum()),significant_fdr05=int((result.padj<.05).sum()),significant_fdr05_lfc1=int(((result.padj<.05)&(result.log2FoldChange.abs()>=1)).sum()),up=int(((result.padj<.05)&(result.log2FoldChange>0)).sum()),down=int(((result.padj<.05)&(result.log2FoldChange<0)).sum()),design_rank=int(np.linalg.matrix_rank(dm)),design_columns=dm.shape[1],fit_diagnostics=diagnostics,candidates=candidates.gene_name.tolist())
    (OUT/'expression_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

if __name__=='__main__': main()

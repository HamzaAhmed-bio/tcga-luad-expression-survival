"""Independent identity checks and prespecified-candidate sensitivity analyses."""
import json,os,pickle
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from lifelines import CoxPHFitter
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';WORK=Path(os.environ.get('TCGA_WORKDIR',str(ROOT.parents[1]/'work/tcga-luad')))
def bh(p):
 p=np.asarray(p,float);i=np.argsort(p);q=np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];a=np.empty(len(p));a[i]=np.minimum(q,1);return a
def main():
 manifest=pd.read_csv(ROOT/'metadata/download_manifest.csv');pairs=manifest[manifest.paired_analysis].copy();pairs['plate']=pairs.aliquot.str.split('-').str[-2];pairs['center']=pairs.aliquot.str.split('-').str[-1]
 assert manifest.file_id.is_unique and manifest['sample'].is_unique
 assert pairs.groupby('patient').size().eq(2).all()
 assert pairs.groupby('patient').sample_type.nunique().eq(2).all()
 assert manifest[manifest.sample_type=='Primary Tumor'].patient.is_unique
 candidates=pd.read_csv(OUT/'frozen_candidates.csv');de=pd.read_csv(OUT/'differential_expression_all_genes.csv.gz');assert candidates.fit_converged.all()
 assert de.loc[~de.fit_converged,'padj'].isna().all()
 with (WORK/'matrices.pkl').open('rb') as f:counts,tpm,ann=pickle.load(f)
 sf=pd.read_csv(OUT/'normalization_factors.csv').set_index('sample').size_factor
 log=np.log2(counts.loc[pairs['sample'],candidates.gene_id].div(sf,axis=0)+1)
 log.columns=candidates.gene_name
 identical=set(pairs.groupby('patient').plate.nunique().loc[lambda x:x==1].index)
 matched=pairs.pivot(index='patient',columns='sample_type',values='sample')
 effects=[]
 for group,pp in [('all_pairs',matched),('same_plate_pairs',matched.loc[sorted(identical)])]:
  for gene in candidates.gene_name:
   tumor=log.loc[pp['Primary Tumor'],gene].to_numpy();normal=log.loc[pp['Solid Tissue Normal'],gene].to_numpy();diff=tumor-normal
   w=wilcoxon(diff,alternative='two-sided',method='auto')
   expected=float(candidates.set_index('gene_name').loc[gene,'log2FoldChange'])
   effects.append(dict(subset=group,gene=gene,pairs=len(pp),median_log2_normalized_difference=float(np.median(diff)),p=float(w.pvalue),direction_agrees=bool(np.sign(np.median(diff))==np.sign(expected))))
 sensitivity=pd.DataFrame(effects);sensitivity['q']=sensitivity.groupby('subset').p.transform(lambda a:bh(a));sensitivity.to_csv(OUT/'paired_candidate_sensitivity.csv',index=False)
 tcga=pd.read_csv(OUT/'TCGA_survival.csv');geo=pd.read_csv(OUT/'GSE31210_survival.csv'); assert len(tcga)==len(geo)==10
 assert np.allclose(tcga.q,bh(tcga.p.fillna(1))) and np.allclose(geo.q,bh(geo.p.fillna(1)))
 # Stage PH diagnostics motivate a clearly labeled, post-fit sensitivity model.
 clinical=pd.read_csv(OUT/'tcga_survival_eligibility.csv',index_col=0)
 tx=pd.read_csv(OUT/'tcga_candidate_log2_tpm.csv',index_col=0)
 stratified=[]
 for gene in candidates.gene_name:
  z=clinical[clinical.eligible].join(tx[gene]).copy();z['z']=(z[gene]-z[gene].mean())/z[gene].std();z['age10']=z.age/10
  fit=CoxPHFitter().fit(z[['time','event','z','age10','stage']],'time','event',strata=['stage'],formula='z + age10');s=fit.summary.loc['z']
  stratified.append(dict(gene=gene,n=len(z),events=int(z.event.sum()),HR=s['exp(coef)'],CI_low=s['exp(coef) lower 95%'],CI_high=s['exp(coef) upper 95%'],p=s.p))
 st=pd.DataFrame(stratified);st['q']=bh(st.p);st.to_csv(OUT/'TCGA_stage_stratified_sensitivity.csv',index=False)
 receipt=dict(identity_checks_passed=True,paired_patients=matched.shape[0],same_plate_pairs=len(identical),different_plate_pairs=matched.shape[0]-len(identical),source_centers=pairs.center.unique().tolist(),workflow_versions=manifest.workflow_version.unique().tolist(),nonconverged_tested_genes=int((de.passed_expression_filter&~de.fit_converged).sum()),candidate_sensitivity_direction_agreements=int(sensitivity.direction_agrees.sum()),candidate_sensitivity_tests=len(sensitivity),candidate_sensitivity_fdr05=int((sensitivity.q<.05).sum()),survival_bh_independently_verified=True,tcga_fit_status=tcga.status.value_counts().to_dict(),geo_fit_status=geo.status.value_counts().to_dict())
 (OUT/'validation_receipt.json').write_text(json.dumps(receipt,indent=2)); print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()

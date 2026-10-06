"""Frozen-candidate survival, external replication, and pathway analysis."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[key]='1'
import csv,gzip,hashlib,io,json,pickle,re,warnings,zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter,KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.plotting import add_at_risk_counts
import gseapy

ROOT=Path(__file__).resolve().parents[1]; WORK=Path(os.environ.get('TCGA_WORKDIR',str(ROOT.parents[1]/'work/tcga-luad'))); EXT=WORK/'external'; OUT=ROOT/'results'; FIG=OUT/'figures'
def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]; result=np.empty(len(p));result[order]=np.minimum(q,1);return result
def save(name):
    plt.tight_layout();plt.savefig(FIG/(name+'.png'),dpi=180,bbox_inches='tight');plt.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close()
def read_xlsx(path):
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        shared=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
        rows=[]
        for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row',ns):
            d={}
            for c in row:
                col=re.sub('[0-9]','',c.get('r'));v=c.find('s:v',ns)
                d[col]=shared[int(v.text)] if c.get('t')=='s' and v is not None else v.text if v is not None else ''
            rows.append(d)
    headers=rows[0];return pd.DataFrame([{name:r.get(col,'') for col,name in headers.items() if name} for r in rows[1:]])
def read_geo():
    lines=gzip.open(EXT/'GSE31210_series_matrix.txt.gz','rt').readlines()
    meta={}; chars=[]
    for line in lines:
        row=next(csv.reader([line],delimiter='\t'))
        if row and row[0].startswith('!Sample_'):
            if row[0]=='!Sample_characteristics_ch1': chars.append(row[1:])
            else: meta[row[0]]=row[1:]
    ids=meta['!Sample_geo_accession'];records=[]
    for i,sample in enumerate(ids):
        r={'sample':sample,'title':meta['!Sample_title'][i]}
        for a in chars:
            if ': ' in a[i]:
                key,value=a[i].split(': ',1)
                # GEO contains a duplicate mislabeled months field; use documented days fields only.
                if key not in r:r[key]=value
        records.append(r)
    clinical=pd.DataFrame(records).set_index('sample');assert clinical.index.is_unique and clinical.title.is_unique
    start=next(i for i,x in enumerate(lines) if x.startswith('!series_matrix_table_begin'))+1
    end=next(i for i,x in enumerate(lines) if x.startswith('!series_matrix_table_end'))
    matrix=pd.read_csv(io.StringIO(''.join(lines[start:end])),sep='\t',index_col=0)
    assert matrix.columns.tolist()==ids and matrix.index.is_unique
    raw=matrix.to_numpy(float); assert np.isfinite(raw).all() and (raw>0).all()
    # Source states MAS5-normalized values; confirm linear intensity range.
    assert np.quantile(raw,.99)>100
    matrix=np.log2(matrix)
    ann_lines=gzip.open(EXT/'GPL570.annot.gz','rt').readlines();a=next(i for i,x in enumerate(ann_lines) if x.startswith('ID\t'))
    ann=pd.read_csv(io.StringIO(''.join(x for x in ann_lines[a:] if not x.startswith('!'))),sep='\t',dtype=str).set_index('ID')
    sym=ann['Gene symbol'].dropna();sym=sym[~sym.str.contains('///',regex=False)&sym.ne('')]
    common=matrix.index.intersection(sym.index);expr=matrix.loc[common].assign(symbol=sym.loc[common]).groupby('symbol').mean().T
    clinical.to_csv(OUT/'geo_metadata.csv');sym.loc[common].to_csv(OUT/'geo_probe_mapping.csv')
    return expr,clinical
def fit_cohort(expression,clinical,label,candidates,notification=None):
    rows=[];diagnostics=[];sens=[];fits={}
    for gene in candidates:
        base={'cohort':label,'gene':gene,'p':1.0,'status':'unavailable'}
        if gene not in expression.columns: rows.append(base);continue
        x=clinical.join(expression[gene].rename('expression'),how='inner')
        x=x[x.time.gt(0)&x.event.isin([0,1])].copy()
        x=x.dropna(subset=['age','stage','expression'])
        # Prespecified complete-case cohort, even when event count limits covariate adjustment.
        x['z']=(x.expression-x.expression.mean())/x.expression.std(ddof=1);x['age10']=x.age/10
        stages=sorted(x.stage.unique());ncoef=2+len(stages)-1;events=int(x.event.sum())
        adjusted=events>=10*ncoef
        formula='z + age10 + C(stage)' if adjusted else 'z'
        d=x[['time','event','z','age10','stage']].copy()
        base.update(n=len(d),events=events,model=formula,adjusted=adjusted,expression_mean=float(x.expression.mean()),expression_sd=float(x.expression.std(ddof=1)))
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always');fit=CoxPHFitter().fit(d,'time','event',formula=formula)
            warning_text='; '.join(str(w.message) for w in caught)
            s=fit.summary.loc['z']; ph=proportional_hazard_test(fit,d,time_transform='rank').summary
            for cov,row in ph.iterrows():diagnostics.append(dict(cohort=label,gene=gene,covariate=cov,p=row.p,test_statistic=row.test_statistic))
            base.update(status='ok' if not warning_text else 'fit_warning',warning=warning_text,HR=s['exp(coef)'],CI_low=s['exp(coef) lower 95%'],CI_high=s['exp(coef) upper 95%'],p=s.p,PH_expression_p=float(ph.loc['z','p']),PH_any_min_p=float(ph.p.min()))
            q=d.copy();q['z2']=q.z**2
            nonlinear=CoxPHFitter().fit(q,'time','event',formula=formula+' + z2')
            base['nonlinearity_p']=float(stats.chi2.sf(max(0,2*(nonlinear.log_likelihood_-fit.log_likelihood_)),1))
            fits[gene]=(x,fit)
            if notification is not None:
                clean=d.loc[~d.index.isin(notification)]
                if clean.event.sum()>=10*ncoef:
                    cf=CoxPHFitter().fit(clean,'time','event',formula=formula);ss=cf.summary.loc['z'];sens.append(dict(cohort=label,gene=gene,n=len(clean),events=int(clean.event.sum()),HR=ss['exp(coef)'],p=ss.p,analysis='Exclude all annotation-positive patients; original expression scaling'))
        except Exception as exc: base.update(status='failed',error=str(exc),p=1.0)
        rows.append(base)
    result=pd.DataFrame(rows);result['q']=bh(result.p.fillna(1));result['PH_expression_q']=bh(result.get('PH_expression_p',pd.Series(1,index=result.index)).fillna(1));result['nonlinearity_q']=bh(result.get('nonlinearity_p',pd.Series(1,index=result.index)).fillna(1))
    result.to_csv(OUT/(label+'_survival.csv'),index=False)
    pd.DataFrame(diagnostics).to_csv(OUT/(label+'_PH_diagnostics.csv'),index=False)
    if sens:pd.DataFrame(sens).to_csv(OUT/(label+'_sensitivity.csv'),index=False)
    return result,fits
def main():
    candidates=pd.read_csv(OUT/'frozen_candidates.csv');names=candidates.gene_name.tolist();assert len(names)==10
    freeze=json.loads((OUT/'candidate_freeze.json').read_text());assert hashlib.sha256((OUT/'frozen_candidates.csv').read_bytes()).hexdigest()==freeze['sha256']
    with (WORK/'matrices.pkl').open('rb') as f:counts,tpm,annotation=pickle.load(f)
    manifest=pd.read_csv(ROOT/'metadata/download_manifest.csv');tumor=manifest[manifest.sample_type=='Primary Tumor'].set_index('sample')
    tx=np.log2(tpm.loc[tumor.index,candidates.gene_id]+1);tx.columns=names;tx.index=tumor.patient;assert tx.index.is_unique
    tx.to_csv(OUT/'tcga_candidate_log2_tpm.csv')
    clinical=read_xlsx(EXT/'TCGA-CDR-SupplementalTableS1.xlsx');clinical=clinical[clinical.type=='LUAD'].set_index('bcr_patient_barcode');assert clinical.index.is_unique
    clinical.to_csv(OUT/'tcga_cdr_luad.csv')
    c=pd.DataFrame(index=clinical.index);c['time']=pd.to_numeric(clinical['OS.time'],errors='coerce');c['event']=pd.to_numeric(clinical.OS,errors='coerce');c['age']=pd.to_numeric(clinical.age_at_initial_pathologic_diagnosis,errors='coerce');c['stage']=clinical.ajcc_pathologic_tumor_stage.str.extract(r'^Stage (IV|III|II|I)')[0]
    c['redaction']=clinical.Redaction.fillna('');redacted=c.redaction.ne(''); c.loc[redacted,['time','event']]=np.nan
    c['has_expression']=c.index.isin(tx.index);c['eligible']=c.has_expression&c.time.gt(0)&c.event.isin([0,1])&c.age.notna()&c.stage.notna()&~redacted
    c.to_csv(OUT/'tcga_survival_eligibility.csv')
    notified=set(tumor.loc[tumor.annotations.notna(),'patient'])
    tcga,tfits=fit_cohort(tx,c,'TCGA',names,notified)
    gx,gmeta=read_geo();gx.reindex(columns=names).to_csv(OUT/'geo_candidate_log2_intensity.csv');g=pd.DataFrame(index=gmeta.index);g['time']=pd.to_numeric(gmeta['days before death/censor'],errors='coerce');g['event']=gmeta.death.map({'dead':1,'alive':0});g['age']=pd.to_numeric(gmeta['age (years)'],errors='coerce');g['stage']=gmeta['pathological stage'].str.extract(r'^(II|I)')[0]
    exclusion='exclude for prognosis analysis due to incomplete resection or adjuvant therapy'
    g['excluded_by_submitter']=gmeta[exclusion].ne('none');g['eligible']=g.time.gt(0)&g.event.notna()&g.age.notna()&g.stage.notna()&~g.excluded_by_submitter;g.to_csv(OUT/'geo_survival_eligibility.csv')
    geo,gfits=fit_cohort(gx,g[g.eligible],'GSE31210',names)
    er=[]
    for name in names:
        if name not in gx:er.append(dict(gene=name,p=1,status='not_mapped'));continue
        a=gx.loc[gmeta.tissue.eq('primary lung tumor'),name];b=gx.loc[gmeta.tissue.eq('normal lung'),name];test=stats.ttest_ind(a,b,equal_var=False)
        er.append(dict(gene=name,n_tumor=len(a),n_normal=len(b),mean_log2_difference=float(a.mean()-b.mean()),p=float(test.pvalue),status='exploratory_unpaired'))
    external=pd.DataFrame(er);external['q']=bh(external.p);external=external.merge(candidates[['gene_name','log2FoldChange']],left_on='gene',right_on='gene_name');external['direction_agrees']=np.sign(external.mean_log2_difference)==np.sign(external.log2FoldChange);external.to_csv(OUT/'external_expression_validation.csv',index=False)
    combined=tcga.merge(geo,on='gene',suffixes=('_TCGA','_GEO'));combined['direction_agrees']=np.sign(np.log(combined.HR_TCGA))==np.sign(np.log(combined.HR_GEO));combined['statistical_replication']=(combined.q_TCGA<.05)&(combined.q_GEO<.05)&combined.direction_agrees&combined.adjusted_TCGA&combined.adjusted_GEO
    combined['diagnostic_flag']=(combined.PH_expression_q_TCGA<.05)|(combined.PH_expression_q_GEO<.05)|(combined.nonlinearity_q_TCGA<.05)|(combined.nonlinearity_q_GEO<.05)|combined.status_TCGA.ne('ok')|combined.status_GEO.ne('ok');combined.to_csv(OUT/'survival_replication.csv',index=False)
    plt.figure(figsize=(9,7));ax=plt.gca()
    for offset,table,color,label in [(-.13,tcga,'#158a8a','TCGA'),(.13,geo,'#cb5264','Independent GSE31210')]:
        table=table.set_index('gene').reindex(names);y=np.arange(len(names))+offset
        ax.errorbar(table.HR,y,xerr=[table.HR-table.CI_low,table.CI_high-table.HR],fmt='o',color=color,label=label,capsize=3)
    ax.set_yticks(range(len(names)),names);ax.set_xscale('log');ax.axvline(1,c='grey',ls='--');ax.set_xlabel('Overall-survival hazard ratio per SD expression (95% CI)');ax.set_title('Expression-selected candidates: survival associations');ax.legend();save('05_survival_forest')
    for label,fits in [('TCGA',tfits),('GSE31210',gfits)]:
        gene=names[0]
        if gene not in fits:continue
        x,_=fits[gene];high=x.expression>x.expression.median();fig,ax=plt.subplots(figsize=(8,6));models=[]
        for mask,title,color in [(~high,'At/below median','#158a8a'),(high,'Above median','#cb5264')]:
            km=KaplanMeierFitter().fit(x.loc[mask,'time']/365.25,x.loc[mask,'event'],label=title);km.plot_survival_function(ax=ax,color=color);models.append(km)
        add_at_risk_counts(*models,ax=ax);ax.set_title(f'{label}: {gene} (descriptive median split)');ax.set_xlabel('Years');ax.set_ylabel('Overall survival');ax.set_ylim(0,1.02);save('06_KM_'+label)
    de=pd.read_csv(OUT/'differential_expression_all_genes.csv.gz').dropna(subset=['stat','gene_name']);de=de[de.fit_converged & np.isfinite(de.stat)].sort_values('baseMean',ascending=False).drop_duplicates('gene_name');rank=de[['gene_name','stat']].sort_values(['stat','gene_name'],ascending=[False,True]);rank.to_csv(OUT/'pathway_rank.csv',index=False)
    gs=gseapy.prerank(rnk=rank,gene_sets=str(EXT/'MSigDB_Hallmark_2020.gmt'),threads=4,min_size=15,max_size=500,permutation_num=1000,outdir=str(OUT/'pathways'),seed=20261005,no_plot=True,verbose=True)
    pathways=gs.res2d;pathways.to_csv(OUT/'pathway_enrichment.csv',index=False)
    top=pathways.assign(NES=pd.to_numeric(pathways.NES),q=pd.to_numeric(pathways['FDR q-val'])).sort_values('q').head(12).sort_values('NES')
    plt.figure(figsize=(10,6));plt.barh(top.Term,top.NES,color=np.where(top.NES>0,'#cb5264','#158a8a'));plt.xlabel('Normalized enrichment score (positive: tumor)');plt.title('Hallmark pathways: 12 lowest enrichment FDR values');save('07_pathways')
    summary=dict(tcga_cdr_cases=len(c),tcga_expression_patients=len(tx),tcga_complete_cases=int(c.eligible.sum()),tcga_events=int(c.loc[c.eligible,'event'].sum()),geo_total_samples=len(gmeta),geo_survival_eligible=int(g.eligible.sum()),geo_events=int(g.loc[g.eligible,'event'].sum()),tcga_survival_fdr05=int((tcga.q<.05).sum()),geo_survival_fdr05=int((geo.q<.05).sum()),statistical_survival_replications=combined.loc[combined.statistical_replication,'gene'].tolist(),diagnostic_flags=combined.loc[combined.diagnostic_flag,'gene'].tolist(),external_expression_direction_agrees=int(external.direction_agrees.sum()),external_expression_fdr05=int((external.q<.05).sum()),hallmark_fdr05=int((pd.to_numeric(pathways['FDR q-val'])<.05).sum()))
    (OUT/'followup_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    provenance=[]
    for p in EXT.iterdir():
        if p.is_file():provenance.append(dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (ROOT/'metadata/external_checksums.json').write_text(json.dumps(provenance,indent=2))

if __name__=='__main__':main()

"""Audit GDC sample identity and download a frozen, checksum-verified cohort."""
import collections, concurrent.futures, csv, hashlib, json, os, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get('TCGA_WORKDIR',str(ROOT.parents[1] / 'work' / 'tcga-luad'))) / 'counts'

def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def main():
    obj = json.loads((ROOT/'metadata/gdc-files-expanded.json').read_text(encoding='utf-8-sig'))
    assert len(obj['data']['hits']) == obj['data']['pagination']['total']
    rows=[]
    for f in obj['data']['hits']:
        assert len(f['cases']) == 1
        c=f['cases'][0]; assert len(c['samples']) == 1
        s=c['samples'][0]
        aliquots=[a for p in s.get('portions',[]) for n in p.get('analytes',[]) for a in n.get('aliquots',[])]
        assert len(aliquots)==1
        categories=sorted({a['category'] for a in c.get('annotations',[])})
        rows.append(dict(file_id=f['file_id'],file_name=f['file_name'],md5sum=f['md5sum'],bytes=f['file_size'],patient=c['submitter_id'],sample=s['submitter_id'],sample_id=s['sample_id'],sample_type=s['sample_type'],aliquot=aliquots[0]['submitter_id'],workflow_version=f['analysis']['workflow_version'],annotations='; '.join(categories),status=''))
    multiplicity=collections.Counter((r['patient'],r['sample_type']) for r in rows)
    for r in rows:
        if r['sample_type'] not in ['Primary Tumor','Solid Tissue Normal']: r['status']='exclude_recurrent'
        elif multiplicity[r['patient'],r['sample_type']] != 1: r['status']='exclude_multiple_files_or_samples_per_patient_condition'
        elif any(term.lower() in r['annotations'].lower() for term in ['neoadjuvant','unacceptable','qualification metrics','may not meet']): r['status']='exclude_quality_or_prior_treatment_annotation'
        else: r['status']='eligible'
    eligible=[r for r in rows if r['status']=='eligible']
    groups=collections.defaultdict(set)
    for r in eligible: groups[r['patient']].add(r['sample_type'])
    paired={p for p,g in groups.items() if len(g)==2}
    for r in rows: r['paired_analysis']=r['status']=='eligible' and r['patient'] in paired
    selected=[r for r in rows if r['status']=='eligible' and (r['sample_type']=='Primary Tumor' or r['paired_analysis'])]
    write_csv(ROOT/'metadata/cohort_audit.csv',rows)
    write_csv(ROOT/'metadata/download_manifest.csv',selected)
    summary=dict(inventory_files=len(rows),selected_files=len(selected),paired_patients=len(paired),primary_tumors=sum(r['sample_type']=='Primary Tumor' for r in selected),exclusions=dict(collections.Counter(r['status'] for r in rows)),download_bytes=sum(r['bytes'] for r in selected),selection_policy='Exclude every ambiguous patient-condition group; no arbitrary aliquot selection; exclude flagged neoadjuvant, unacceptable treatment, qualification changes, or protocol deviation.')
    (ROOT/'metadata/cohort_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)
    RAW.mkdir(parents=True,exist_ok=True)
    def download(r):
        dest=RAW/(r['file_id']+'.tsv')
        if dest.exists() and hashlib.md5(dest.read_bytes()).hexdigest()==r['md5sum']: return r['file_id']
        for attempt in range(5):
            try:
                with urllib.request.urlopen('https://api.gdc.cancer.gov/data/'+r['file_id'],timeout=120) as response: data=response.read()
                assert len(data)==r['bytes'] and hashlib.md5(data).hexdigest()==r['md5sum'], 'Checksum mismatch'
                dest.write_bytes(data); return r['file_id']
            except Exception:
                if attempt==4: raise
                time.sleep(2**attempt)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        for i,_ in enumerate(pool.map(download,selected),1):
            if i%20==0 or i==len(selected): print(f'Verified {i}/{len(selected)} expression files',flush=True)
    (ROOT/'metadata/download_receipt.json').write_text(json.dumps(dict(files=len(selected),all_md5_verified=True,completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2))

if __name__=='__main__': main()

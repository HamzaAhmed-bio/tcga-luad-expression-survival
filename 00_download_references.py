"""Download source references and verify them against the delivered snapshot."""
import hashlib,json,os,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXT=Path(os.environ.get('TCGA_WORKDIR',str(ROOT.parents[1]/'work/tcga-luad')))/'external'
SOURCES={
 'TCGA-CDR-SupplementalTableS1.xlsx':'https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81',
 'GSE31210_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz',
 'GPL570.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz',
 'MSigDB_Hallmark_2020.gmt':'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020',
}
def main():
 EXT.mkdir(parents=True,exist_ok=True)
 checksum_file=ROOT/'metadata/external_checksums.json'
 expected={r['file']:r['sha256'] for r in json.loads(checksum_file.read_text())} if checksum_file.exists() else {}
 for name,url in SOURCES.items():
  dest=EXT/name
  if not dest.exists():
   with urllib.request.urlopen(url,timeout=180) as response:data=response.read()
   if name in expected:assert hashlib.sha256(data).hexdigest()==expected[name],f'Source changed: {name}'
   dest.write_bytes(data)
  if name in expected:assert hashlib.sha256(dest.read_bytes()).hexdigest()==expected[name],f'Checksum failed: {name}'
  print('Verified',name,flush=True)
if __name__=='__main__':main()

"""Hash/config acceptance only: no R, package installation, SSH or sbatch."""
from pathlib import Path
import argparse,csv,datetime,getpass,hashlib,json,os,platform,subprocess,sys
REQUIRED={'TSR_ACCOUNT':'hbulut','TSR_QOS':'normal','TSR_PARTITION':'debug',
 'TSR_GCC_MODULE':'comp/gcc/14.1.0','TSR_R_MODULE':'apps/R/4.5.2-gcc-14.1.0','TSR_ENV_VERIFIED':'1'}
def validate_files(project,manifest_path,revision):
 manifest=json.loads(manifest_path.read_text());assert manifest['schema']=='tsr-transfer-v1'
 listed=set((revision/'transfer/files_from_project.txt').read_text().splitlines())
 assert listed=={x['path'] for x in manifest['files']}|{'revision_scientific_reports/transfer/payload_manifest.json'},'Whitelist/manifest mismatch'
 for item in manifest['files']:
  relative=Path(item['path']);assert not relative.is_absolute() and '..' not in relative.parts,'Unsafe path'
  p=project/relative
  assert p.is_file() and not p.is_symlink(),'Missing/non-regular file: '+str(relative)
  assert os.path.commonpath([str(project),str(p.resolve())])==str(project),'File escapes bundle'
  assert p.stat().st_size==item['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],'SHA256 mismatch: '+str(relative)
 for pin in csv.DictReader((revision/'config/runtime_md5.csv').open()):
  p=(revision/pin['path']).resolve()
  assert os.path.commonpath([str(project),str(p)])==str(project),'Runtime pin outside bundle'
  assert hashlib.md5(p.read_bytes()).hexdigest()==pin['md5'],'Frozen runtime mismatch: '+pin['path']
 return manifest

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--audit-only',action='store_true');parser.add_argument('--verify-only',action='store_true');a=parser.parse_args()
 assert not(a.audit_only and a.verify_only)
 revision=Path(__file__).resolve().parents[1];project=revision.parent
 for key,value in REQUIRED.items():assert os.environ.get(key)==value,'Configuration mismatch/unset: '+key+' (source config/truba.env)'
 assert Path(os.environ['TSR_REVISION_ROOT']).resolve()==revision,'Wrong revision root'
 assert Path(os.environ['TSR_SOURCE_ROOT']).resolve()==project/'Scientific Reports Manuscript','Wrong source root'
 assert Path(os.environ['TSR_OUTPUT_ROOT']).resolve()==revision/'results','Wrong output root'
 lib=Path(os.environ['TSR_R_LIB']).resolve()
 assert os.path.commonpath([str(revision/'results/dependencies'),str(lib)])==str(revision/'results/dependencies'),'Library must be a separate project library'
 assert 'R-4.5.2-gcc-14.1.0-linux' in str(lib),'Wrong R/compiler library namespace'
 manifest_path=revision/'transfer/payload_manifest.json';manifest=validate_files(project,manifest_path,revision)
 if a.audit_only:
  print('LOCAL_AUDIT_OK',len(manifest['files'])+1,'files; no remote acceptance receipt written');return
 assert platform.system()=='Linux','Remote acceptance must run on Linux; use --audit-only locally'
 assert subprocess.check_output(['id','-un'],text=True).strip()=='hbulut','Acceptance must run as hbulut'
 receipt_path=revision/'results/environment/transfer_acceptance.json'
 fingerprint=dict(manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
  config_sha256=hashlib.sha256((revision/'config/truba.env').read_bytes()).hexdigest(),
  runtime_pins_sha256=hashlib.sha256((revision/'config/runtime_md5.csv').read_bytes()).hexdigest(),settings=REQUIRED,revision_root=str(revision))
 if a.verify_only:
  receipt=json.loads(receipt_path.read_text())
  assert all(receipt.get(k)==v for k,v in fingerprint.items()),'Stale transfer acceptance; rerun check_transfer.sh'
  print('TRANSFER_RECEIPT_OK');return
 for folder in ['results/environment','results/pre_truba','figures','tables']:(revision/folder).mkdir(parents=True,exist_ok=True)
 fingerprint.update(checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),user='hbulut',platform='Linux')
 receipt_path.write_text(json.dumps(fingerprint,indent=2)+'\n')
 print('TRANSFER_AND_CONFIGURATION_OK',len(manifest['files'])+1,'files; receipt',receipt_path)
 print('No R process, package installation or SLURM submission was performed.')
if __name__=='__main__':
 try:main()
 except (AssertionError,OSError,KeyError,ValueError) as e:print('TRANSFER_CHECK_FAILED:',e,file=sys.stderr);sys.exit(2)

"""Prepare a whitelist and SHA256 manifest. No network, transfer or job submission."""
from pathlib import Path
import csv,hashlib,json
revision=Path(__file__).resolve().parents[1];project=revision.parent
selected=set()
def add(path):
 p=project/path
 if not p.is_file() or p.is_symlink():raise RuntimeError('Missing/non-regular payload file: '+str(path))
 selected.add(Path(path).as_posix())
def files(folder,pattern):
 for p in (project/folder).glob(pattern):
  if p.is_file():add(p.relative_to(project))
# Submitted package source and fixed-fixture tests; no nested outputs/reproducibility copy.
pkg=Path('Scientific Reports Manuscript/RPackageNew/TSRClust')
for name in ['DESCRIPTION','NAMESPACE','LICENSE','README.md','CITATION']:add(pkg/name)
for directory,pattern in [('R','*.R'),('man','*.Rd'),('inst','*'),('tests','testthat.R'),('tests/testthat','*.R')]:files(pkg/directory,pattern)
# All loader dependencies, exactly the versioned submitted implementation/config/metrics.
source=Path('Scientific Reports Manuscript/reproducibility')
for folder in ['core','competitors','config','evaluation']:files(source/folder,'*.R')
rev=Path('revision_scientific_reports')
add(rev/'README.md')
for folder,patterns in [('R',['*.R']),('scripts',['*.R','*.py','*.sh']),('tests',['*.R','*.py']),('slurm',['*.sh','*.sbatch']),('config',['*.R','*.csv','truba.env','truba.env.example']),('reports',['*.md'])]:
 for pattern in patterns:files(rev/folder,pattern)
for name in ['local_tests.csv','workflow_checks.csv','public_local_alignment.csv','github_provenance.json','source_preservation_pre_truba.csv','truba_transfer_static_checks.csv']:
 add(rev/'reports'/name)
for name in ['workload.csv','workload_by_cell.csv']:add(rev/'tables'/name)
deps=rev/'config/dependencies'
for name in ['README.md','source_lock.csv','macos_binary_lock.csv']:add(deps/name)
files(deps,'*.DESCRIPTION')
lock=list(csv.DictReader((project/deps/'source_lock.csv').open()))
for entry in lock:
 add(deps/entry['source_file'])
 p=project/deps/entry['source_file']
 if hashlib.sha256(p.read_bytes()).hexdigest()!=entry['sha256']:raise RuntimeError('Archive SHA256 mismatch: '+entry['package'])
# The whitelist includes itself and the integrity manifest. The manifest cannot hash itself.
transfer=revision/'transfer';transfer.mkdir(exist_ok=True)
list_relative=rev/'transfer/files_from_project.txt';manifest_relative=rev/'transfer/payload_manifest.json'
all_files=sorted(selected|{list_relative.as_posix(),manifest_relative.as_posix()})
(project/list_relative).write_text('\n'.join(all_files)+'\n')
records=[]
for relative in sorted(selected|{list_relative.as_posix()}):
 p=project/relative;records.append(dict(path=relative,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
manifest=dict(schema='tsr-transfer-v1',created_date='2026-11-19',files=records,
 source_lock_sha256=hashlib.sha256((project/deps/'source_lock.csv').read_bytes()).hexdigest(),
 excluded=['all results/ trees','local R libraries','macOS package binaries','snapshots','historical projects','Git metadata','cached PACKAGES index','figures','temporary files'])
(project/manifest_relative).write_text(json.dumps(manifest,indent=2)+'\n')
size=sum(x['bytes'] for x in records)+(project/manifest_relative).stat().st_size
print('TRANSFER_BUNDLE_READY',len(all_files),'files',size,'bytes',round(size/1048576,2),'MiB')
print('No transfer or sbatch call was made.')

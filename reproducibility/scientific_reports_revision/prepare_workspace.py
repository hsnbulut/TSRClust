"""Prepare a disposable, portable workspace for existing read-only derivations.
Does not fit models, generate observations, install packages, or submit jobs.
"""
from pathlib import Path
import sys,shutil,gzip,csv,hashlib
D=Path(__file__).resolve().parent;G=D.parents[1]
if len(sys.argv)!=2:raise SystemExit('Usage: python3 prepare_workspace.py EMPTY_OUTPUT_DIRECTORY')
W=Path(sys.argv[1]).resolve()
if W.exists() and any(W.iterdir()):raise SystemExit('Use an empty directory; existing work is never overwritten')
W.mkdir(parents=True,exist_ok=True)
def copy(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 if src.is_dir():shutil.copytree(src,dst,dirs_exist_ok=True)
 else:shutil.copy2(src,dst)
def link(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True);dst.symlink_to(src.resolve(),target_is_directory=src.is_dir())
copy(D/'analysis_sources',W/'Revision1')
copy(D/'runtime',W/'revision_scientific_reports')
for src in (D/'tables').rglob('*.csv.gz'):
 rel=src.relative_to(D/'tables');dst=W/'Revision1'/rel.parent/'tables'/rel.with_suffix('').name;dst.parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(src,'rb') as f:dst.write_bytes(f.read())
for c in (W/'Revision1').glob('reviewer*/comment*'):
 for name in ['tables','evidence','evidence/qa','figures']:(c/name).mkdir(parents=True,exist_ok=True)
for src in (D/'inputs').glob('reviewer*/comment*/*.rds'):copy(src,W/'Revision1'/src.relative_to(D/'inputs').parent/'evidence'/src.name)
F=W/'Revision1/reviewer1/comment3/finalization';(F/'evidence').mkdir(parents=True,exist_ok=True);(F/'tables').mkdir(exist_ok=True)
O=W/'Scientific Reports Manuscript';O.mkdir()
link(G/'reproducibility',O/'reproducibility');link(G/'results',O/'github_results');link(G,O/'RPackageNew/TSRClust')
copy(D/'inputs/original_analysis_R',O/'analysis/R')
for st in ['architecture','stress','negative_control']:copy(D/f'inputs/original_rds/{st}_replications.rds',O/f'analysis/v2/results/final_sim_combined/{st}_replications.rds')
link(G/'results/summary_results',O/'analysis/v2/results/final_sim_summary');link(G/'results/real_data_results',O/'analysis/v2/results/realdata')
link(G/'reproducibility/real_data/02_realdata_v2.R',O/'analysis/v2/scripts/02_realdata_v2.R')
copy(D/'inputs/legacy_diagnostics/automacro_replications.rds',W/'TSRClust_SciReports_TRUBA/results/combined/automacro_replications.rds')
# The private original guards also included manuscripts and correspondence. Public
# guards intentionally cover only available scientific source/data, not those files.
archive=D/'results/TSRClust_Tier1_Analysis.tar.gz'
pins=[]
for p in sorted(W.rglob('*')):
 if p.is_file() and (p.suffix in ['.R','.rds'] or (p.suffix=='.csv' and 'revision_scientific_reports' in p.parts)):
  pins.append(dict(role='public scientific source/data',path=str(p),project_relative=str(p.relative_to(W)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pins.append(dict(role='immutable Tier1 archive',path=str(archive),project_relative='',bytes=archive.stat().st_size,sha256=hashlib.sha256(archive.read_bytes()).hexdigest()))
for c in (W/'Revision1').glob('reviewer*/comment*'):
 for target in [c/'evidence/source_inventory.csv']+([F/'evidence/scientific_inputs.csv'] if c.name=='comment3' and c.parent.name=='reviewer1' else []):
  target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(pins[0]));w.writeheader();w.writerows(pins)
print(W)
print('Set TSR_TIER1_ARCHIVE='+str(archive))
print('Run only the documented read-only analysis commands; no simulation launch is implied.')

from pathlib import Path
import csv,hashlib,sys,os
r=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((r/'config/scientific_integrity.csv').open()))
for x in rows:assert hashlib.sha256((r/x['path']).read_bytes()).hexdigest()==x['sha256'],'Frozen science changed: '+x['path']
if '--local-originals' in sys.argv:
 for x in csv.DictReader((r/'config/source_sha256.csv').open()):
  if x['sha256']!='UNAVAILABLE_CLOUD_FILE':assert hashlib.sha256((r.parent/x['path']).read_bytes()).hexdigest()==x['sha256'],'Original source/result changed: '+x['path']
print('SCIENTIFIC_METHOD_ORACLE_METRICS_AND_FROZEN_DESIGN_UNCHANGED',len(rows),'hashes')

out=Path(os.environ.get('TSR_ACCEPTANCE_RUN_DIR',r/'results/final_acceptance/local_20261119'));out.mkdir(parents=True,exist_ok=True)
with (out/'scientific_preservation.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['check','status','hash_count']);w.writerow(['frozen method/design/submitted package','PASS',len(rows)])

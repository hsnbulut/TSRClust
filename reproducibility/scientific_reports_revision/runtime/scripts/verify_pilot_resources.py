from pathlib import Path
import csv,json,sys
out=Path(sys.argv[1]);tasks=list(csv.DictReader(open(sys.argv[2])));rows=[]
for t in tasks:
 path=out/f"task_{int(t['task_id']):06d}_resources.json"
 if not path.exists():raise SystemExit('Missing resource file: '+str(path))
 z=json.loads(path.read_text());assert z['exit_code']==0 and z['process_peak_rss_bytes']>0 and z['cpu_seconds']>=0
 rows.append({k:z[k] for k in ['task_id','wall_seconds','cpu_seconds','process_peak_rss_bytes','exit_code']})
with (out/'task_resources.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('PILOT_PROCESS_RESOURCES_OK',len(rows),'process measurements; SLURM MaxRSS must additionally be collected with sacct')

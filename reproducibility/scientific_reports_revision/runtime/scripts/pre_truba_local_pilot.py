"""Only the existing 20-cell/two-replication pilot; never a Tier 1 launcher."""
from pathlib import Path
import subprocess,os,time,csv,json,sys,resource
root=Path(__file__).resolve().parents[1];os.chdir(root)
out=Path(os.environ.get('TSR_LOCAL_PILOT_OUTPUT','results/pilot_final_acceptance_local_20261119'));(out/'logs').mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env['TSR_ARTIFACT_DIR']=str((out/'artifacts').resolve());env['R_LIBS_USER']=str((root/'results/dependencies/R-4.5.2-local/library').resolve());env['TSR_R_LIB']=env['R_LIBS_USER']
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS']:env[key]='1'
tasks=list(csv.DictReader(open('config/tasks_pilot.csv')))
assert len(tasks)==20 and all(int(t['rep_end'])-int(t['rep_start'])==1 for t in tasks)
resources=[]
for task in tasks:
 tid=task['task_id'];start=time.monotonic()
 with (out/'logs'/f'task_{int(tid):06d}.log').open('w') as log:
  p=subprocess.Popen(['python3','scripts/run_task_measured.py','config/tasks_pilot.csv',tid,str(out)],env=env,stdout=log,stderr=subprocess.STDOUT)
  _,status,usage=os.wait4(p.pid,0);code=os.waitstatus_to_exitcode(status);p.returncode=code
 rss=usage.ru_maxrss if sys.platform=='darwin' else usage.ru_maxrss*1024
 resources.append(dict(task_id=int(tid),wall_seconds=time.monotonic()-start,child_cpu_seconds=usage.ru_utime+usage.ru_stime,peak_rss_bytes=rss,exit_code=code))
 with (out/'task_resources.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=resources[0].keys());w.writeheader();w.writerows(resources)
 print('PILOT_TASK',tid,'process_exit',code,'wall_seconds',round(resources[-1]['wall_seconds'],2),flush=True)
 if code:raise SystemExit(code)
for script in ['check_complete.R','aggregate.R']:
 subprocess.run(['bash','scripts/tsr_rscript.sh',f'scripts/{script}','config/tasks_pilot.csv',str(out)],env=env,check=True)

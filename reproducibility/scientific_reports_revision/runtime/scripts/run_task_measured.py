"""Run one frozen manifest task through the common R launcher and save atomic resources."""
from pathlib import Path
import argparse,datetime,json,os,subprocess,sys,time
r=Path(__file__).resolve().parents[1];os.chdir(r)
p=argparse.ArgumentParser();p.add_argument('manifest');p.add_argument('task_id',type=int);p.add_argument('output');a=p.parse_args()
out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
checkpoint=out/f'task_{a.task_id:06d}.rds';previous_mtime=checkpoint.stat().st_mtime_ns if checkpoint.exists() else None
start=time.monotonic();proc=subprocess.Popen([str(r/'scripts/tsr_rscript.sh'),'scripts/run_task.R',a.manifest,str(a.task_id),str(out)])
_,status,usage=os.wait4(proc.pid,0);code=os.waitstatus_to_exitcode(status);proc.returncode=code
z={'task_id':a.task_id,'wall_seconds':time.monotonic()-start,'cpu_seconds':usage.ru_utime+usage.ru_stime,
 'process_peak_rss_bytes':usage.ru_maxrss if sys.platform=='darwin' else usage.ru_maxrss*1024,'exit_code':code,
 'slurm':{k:os.environ.get(k,'') for k in ['SLURM_JOB_ID','SLURM_ARRAY_JOB_ID','SLURM_ARRAY_TASK_ID','SLURM_STEP_ID']},
 'successful_checkpoint_reused':previous_mtime is not None and checkpoint.exists() and checkpoint.stat().st_mtime_ns==previous_mtime,
 'measured_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
# Keep earlier process measurements, including a skipped-success retry.
path=out/f'task_{a.task_id:06d}_resources.json'
if path.exists():path.rename(out/(path.stem+'.previous_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%f')+'.json'))
tmp=path.with_suffix('.'+str(os.getpid())+'.tmp');tmp.write_text(json.dumps(z,indent=2)+'\n');tmp.replace(path)
raise SystemExit(code)

"""Measure one available-package pilot task without sysctl or elevated access."""
import csv,os,resource,subprocess,time
os.makedirs('results/resource_probe',exist_ok=True)
start=time.monotonic()
p=subprocess.run(['bash','scripts/tsr_rscript.sh','scripts/run_task.R','config/tasks_pilot.csv','9','results/resource_probe'],capture_output=True,text=True)
r=resource.getrusage(resource.RUSAGE_CHILDREN)
# Darwin ru_maxrss is bytes, Linux is KiB.
import sys
rss_bytes=r.ru_maxrss if sys.platform=='darwin' else r.ru_maxrss*1024
with open('reports/resource_probe.csv','w') as f:
 w=csv.writer(f);w.writerow(['task_id','wall_seconds','peak_child_rss_bytes','return_code']);w.writerow([9,time.monotonic()-start,rss_bytes,p.returncode])
print(p.stdout,p.stderr,sep='');print('peak child RSS bytes',rss_bytes)
raise SystemExit(p.returncode)

"""Local static/contract tests. Never calls a real scheduler or remote host."""
from pathlib import Path
import csv,os,subprocess,tempfile
r=Path(__file__).resolve().parents[1]
rows=[]
def passed(name):rows.append(dict(check=name,status='PASS'));print('PASS',name)
for path in list((r/'slurm').glob('*.sh'))+list((r/'slurm').glob('*.sbatch'))+list((r/'scripts').glob('*.sh')):
 subprocess.run(['bash','-n',str(path)],capture_output=True,check=True)
passed('all shell files parse')
for name in ['array.sbatch','install_and_verify.sbatch']:
 s=(r/'slurm'/name).read_text()
 for v in ['#SBATCH --account=hbulut','#SBATCH --qos=normal','#SBATCH --partition=debug','runtime_environment.sh']:assert v in s
passed('active jobs explicitly use hbulut normal debug and shared runtime')
s=(r/'slurm/runtime_environment.sh').read_text()
assert s.index('SLURM_JOB_ID')<s.index('module load')<s.index('tsr_rscript.sh --probe')
for field in ['TSR_GCC_MODULE','TSR_R_MODULE']:assert 'module load "$'+field+'"' in s
assert s.index('scontrol show hostnames')<s.index('module load')
installer=(r/'scripts/install_dependencies.R').read_text()
assert installer.index('Current host is not an allocated SLURM compute node')<installer.index('install.packages(')
passed('compute allocation and host membership are required before modules, R or Linux installation')
for name in ['submit_common.sh','submit_install.sh']:
 s=(r/'slurm'/name).read_text();assert '--account="$TSR_ACCOUNT"' in s and '--qos="$TSR_QOS"' in s
 assert s.index('check_transfer.py --verify-only')<s.rindex('sbatch ')
passed('submission functions check transfer receipt and pass account and QoS')
s=(r/'slurm/install_and_verify.sbatch').read_text();assert '--install scripts/install_dependencies.R' in s and 'submit_tier' not in s and '--reuse-existing' not in s
passed('acceptance job performs source installation only on compute and never submits studies')
# Replace python3 (receipt checker) and sbatch with test stubs. The rejected receipt
# must prevent any scheduler invocation. Real sbatch/SSH is never reached.
with tempfile.TemporaryDirectory(prefix='tsr-transfer-test-') as td:
 p=Path(td);flag=p/'unexpected_sbatch'
 for name,body in [('python3','#!/bin/sh\nexit 2\n'),('sbatch','#!/bin/sh\ntouch "'+str(flag)+'"\nexit 99\n')]:
  f=p/name;f.write_text(body);f.chmod(0o700)
 env=os.environ.copy();env['PATH']=str(p)+os.pathsep+env['PATH']
 proc=subprocess.run(['bash',str(r/'slurm/submit_install.sh')],env=env,capture_output=True,text=True)
 assert proc.returncode==2 and not flag.exists()
passed('no sbatch call occurs when transfer acceptance fails (stubbed test)')
with (r/'reports/truba_transfer_static_checks.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['check','status']);w.writeheader();w.writerows(rows)

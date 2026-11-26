"""Stubbed acceptance control-flow test. Never invokes real R/SLURM/install."""
from pathlib import Path
import os,subprocess,tempfile,csv
r=Path(__file__).resolve().parents[1];checks=[]
for scenario in ['valid','invalid_existing']:
 with tempfile.TemporaryDirectory(prefix='tsr-acceptance-control-') as td:
  root=Path(td);(root/'slurm').mkdir();(root/'scripts').mkdir();(root/'existing-library').mkdir();(root/'bin').mkdir()
  (root/'slurm/configuration.sh').write_text('export TSR_REVISION_ROOT="'+str(root)+'"\nexport TSR_OUTPUT_ROOT="'+str(root/'results')+'"\nexport TSR_R_LIB="'+str(root/'existing-library')+'"\n')
  (root/'slurm/runtime_environment.sh').write_text('true\n')
  log=root/'calls.txt'
  body='#!/bin/sh\nprintf "%s\\n" "$*" >> "'+str(log)+'"\n'
  if scenario=='invalid_existing':body+='case "$1" in scripts/verify_dependencies.R) exit 3;; esac\n'
  f=root/'scripts/tsr_rscript.sh';f.write_text(body);f.chmod(0o700)
  for name,body in [('python3','#!/bin/sh\nexit 0\n'),('uname','#!/bin/sh\nprintf "%s\\n" Linux\n')]:
   f=root/'bin'/name;f.write_text(body);f.chmod(0o700)
  env=os.environ.copy();env.update(TSR_REVISION_ROOT=str(root),TSR_ACCEPTANCE_VERIFY_ONLY='1',SLURM_JOB_ID='STUB_ONLY',PATH=str(root/'bin')+os.pathsep+env['PATH'])
  p=subprocess.run(['bash',str(r/'slurm/install_and_verify.sbatch')],env=env,capture_output=True,text=True)
  calls=log.read_text().splitlines();assert not any('install_dependencies.R' in x for x in calls)
  assert not any((root/'existing-library').iterdir())
  if scenario=='valid':assert p.returncode==0 and 'DEPENDENCIES_ALREADY_VALID' in p.stdout and any('minimal_competitors.R' in x for x in calls)
  else:assert p.returncode!=0 and len(calls)==1 and 'DEPENDENCY_PRECHECK_FAILED' in p.stderr
  checks.append([scenario,'PASS'])
out=Path(os.environ.get('TSR_ACCEPTANCE_RUN_DIR',r/'results/final_acceptance/local_20261119'));out.mkdir(parents=True,exist_ok=True)
with (out/'acceptance_reuse_contract.csv').open('w') as f:w=csv.writer(f);w.writerow(['check','status']);w.writerows(checks)
print('ACCEPTANCE_PRECHECK_AND_NO_REINSTALL_OK 2 (stubbed)')

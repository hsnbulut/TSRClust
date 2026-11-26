"""Static contract for central launcher routing, exports and frozen pilot coverage."""
from pathlib import Path
import ast,csv,subprocess,os
r=Path(__file__).resolve().parents[1]
for f in list((r/'scripts').glob('*.py'))+list((r/'tests').glob('*.py')):ast.parse(f.read_text())
launcher=(r/'scripts/tsr_rscript.sh').read_text();assert 'Rscript --vanilla "$@"' in launcher and 'assert_revision_environment()' in launcher
array=(r/'slurm/array.sbatch').read_text();assert '--export=ALL' in array and 'scripts/run_task_measured.py' in array
measured=(r/'scripts/run_task_measured.py').read_text();assert 'scripts/tsr_rscript.sh' in measured
for name in ['local_pilot.sh','resource_probe.py','pre_truba_local_pilot.py']:
 s=(r/'scripts'/name).read_text();assert 'tsr_rscript.sh' in s or 'run_task_measured.py' in s
for f in ['tests/workflow_checks.R','tests/library_child_process.R']:assert "normalizePath('scripts/tsr_rscript.sh')" in (r/f).read_text()
runtime=(r/'slurm/runtime_environment.sh').read_text();assert 'SLURM_EXPORT_ENV=ALL' in runtime
cells=list(csv.DictReader((r/'config/design.csv').open()));tasks=list(csv.DictReader((r/'config/tasks_pilot.csv').open()))
assert len(tasks)==20 and sum(int(t['rep_end'])-int(t['rep_start'])+1 for t in tasks)==40
picked=[next(c for c in cells if c['cell_id']==t['cell_id']) for t in tasks]
assert {'legacy_rule','pvalue_only','recalibrated_rule'}<=set(c['rule'] for c in picked)
assert {'OTRIMLE_fixedK','tEIGEN_fixedK','tEIGEN_BIC','TCLUST','TSRClust','DBSCAN_full','HDBSCAN_full','Trimmed_kmeans','Oracle_TCLUST'}<=set(c['method'] for c in picked)
print('CENTRAL_LAUNCHER_ARRAY_EXPORTS_AND_FROZEN_40_FIT_COVERAGE_OK')

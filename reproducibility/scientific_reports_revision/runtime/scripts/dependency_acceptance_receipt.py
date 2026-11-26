"""Hashed acceptance status and success receipt; no R, installation or submission."""
from pathlib import Path
import argparse, csv, datetime, hashlib, json, os, sys

r = Path(__file__).resolve().parents[1]


def fp(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def identity():
    return {'runtime_sha256': fp(r/'config/runtime_md5.csv'),
            'dependency_lock_sha256': fp(r/'config/dependencies/source_lock.csv'),
            'library': str(Path(os.environ['TSR_R_LIB']).resolve())}


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def status(out, state, **details):
    assert os.environ.get('SLURM_JOB_ID'), 'Acceptance status requires a compute acceptance job'
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    z = dict(schema='tsr-dependency-acceptance-status-v1', status=state,
             job_id=os.environ['SLURM_JOB_ID'], output=str(out),
             checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), **details)
    atomic_json(out/'acceptance_status.json', z)
    atomic_json(r/'results/environment/dependency_acceptance_latest.json', z)
    return z


def validate_reports(out):
    required = ['dependency_verification.csv', 'library_child_regression.csv',
                'startup_isolation_regression.csv', 'minimal_competitors.csv',
                'targeted_checks.csv', 'local_tests.csv', 'workflow_checks.csv',
                'scientific_preservation.csv', 'startup_acceptance_contract.csv']
    assert all((out/x).is_file() for x in required), 'Incomplete acceptance test reports'
    for name in required:
        rows = list(csv.DictReader((out/name).open())); assert rows, 'Empty report: ' + name
        if name == 'dependency_verification.csv':
            assert len(rows) == 17 and all(x['loaded'] == 'TRUE' and x['from_isolated_library'] == 'TRUE' and x['version_matches'] == 'TRUE' for x in rows), 'Dependency verification failed'
        else:
            assert all(x['status'] == 'PASS' for x in rows), 'Failed acceptance report: ' + name
    return required


def main():
    p = argparse.ArgumentParser()
    action = p.add_mutually_exclusive_group()
    action.add_argument('--write'); action.add_argument('--begin'); action.add_argument('--record-failure')
    p.add_argument('--stage', default='unknown'); p.add_argument('--exit-code', type=int, default=1)
    p.add_argument('--command', default=''); a = p.parse_args()
    receipt = r/'results/environment/dependency_acceptance.json'
    latest = r/'results/environment/dependency_acceptance_latest.json'
    if a.begin:
        status(a.begin, 'RUNNING', stage='initialization')
        print('DEPENDENCY_ACCEPTANCE_STARTED'); return 0
    if a.record_failure:
        out = Path(a.record_failure).resolve()
        diagnostics = []
        for file in sorted(out.glob('startup_*.json')):
            try:
                z = json.loads(file.read_text())
                if 'returncode' in z:
                    diagnostics.append(dict(path=str(file), **z))
            except (OSError, ValueError) as error:
                diagnostics.append(dict(path=str(file), diagnostic_read_error=str(error)))
        z = status(out, 'FAILED', stage=a.stage, exit_code=a.exit_code, command=a.command,
                   child_diagnostics=diagnostics,
                   dependency_precheck_log=str(out/'dependency_precheck.log'),
                   slurm_stdout=str(r/f'results/logs/install_accept_{os.environ["SLURM_JOB_ID"]}.out'),
                   slurm_stderr=str(r/f'results/logs/install_accept_{os.environ["SLURM_JOB_ID"]}.err'))
        print('DEPENDENCY_FAILURE_RECEIPT_WRITTEN ' + str(out/'acceptance_status.json')); return 0
    if a.write:
        assert os.environ.get('SLURM_JOB_ID'), 'Acceptance receipt requires the real compute acceptance job'
        out = Path(a.write).resolve(); required = validate_reports(out)
        z = identity(); z.update(job_id=os.environ['SLURM_JOB_ID'], output=str(out),
                                reports={x: fp(out/x) for x in required},
                                checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        receipt.parent.mkdir(parents=True, exist_ok=True)
        if receipt.exists():
            receipt.rename(receipt.with_name('dependency_acceptance.before_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%f') + '.json'))
        atomic_json(receipt, z); status(out, 'PASSED', stage='complete', exit_code=0)
        print('DEPENDENCY_RECEIPT_WRITTEN'); return 0
    if latest.exists():
        z = json.loads(latest.read_text())
        if z['status'] != 'PASSED':
            print('DEPENDENCY_ACCEPTANCE_' + z['status'] + ' job=' + z['job_id'] + ' stage=' + z['stage'], file=sys.stderr)
            print(json.dumps(z, indent=2), file=sys.stderr)
            return 2
    if not receipt.is_file():
        print('DEPENDENCY_ACCEPTANCE_NOT_READY: no successful receipt. Inspect results/logs/install_accept_JOB_ID.out/.err and results/acceptance/job_JOB_ID/acceptance_status.json (available for corrected reruns).', file=sys.stderr)
        return 2
    z = json.loads(receipt.read_text())
    assert all(z.get(k) == v for k, v in identity().items()), 'Dependency acceptance is stale for this code/library'
    assert all(fp(Path(z['output'])/x) == h for x, h in z['reports'].items()), 'Acceptance report changed or missing'
    print('DEPENDENCY_ACCEPTANCE_RECEIPT_OK job=' + z['job_id']); return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (AssertionError, OSError, KeyError, ValueError) as error:
        print('DEPENDENCY_ACCEPTANCE_NOT_READY: ' + str(error), file=sys.stderr)
        sys.exit(2)

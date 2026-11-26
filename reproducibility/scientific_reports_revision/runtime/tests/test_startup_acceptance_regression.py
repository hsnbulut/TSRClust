"""Regression for frontend selection, captured 127 diagnostics and failure receipts.
All 127/SLURM cases are synthetic: no scheduler, real R or installed package writes.
"""
from pathlib import Path
import contextlib, csv, importlib.util, io, json, os, shlex, shutil, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('startup_under_test', ROOT/'tests/test_r_startup_isolation.py')
startup = importlib.util.module_from_spec(spec); spec.loader.exec_module(startup)
rows = []


def passed(name):
    rows.append([name, 'PASS']); print('PASS ' + name)


with tempfile.TemporaryDirectory(prefix='tsr-acceptance-regression-') as td:
    t = Path(td).resolve(); home = t/'R home with spaces'; (home/'bin/exec').mkdir(parents=True)
    probe = t/'probe with spaces.R'; probe.write_text('# disposable command-construction fixture\n')
    raw = home/'bin/exec/R'
    raw.write_text('#!' + sys.executable + '\nimport os,sys\n'
                   'if not os.environ.get("SYNTHETIC_LOADER_READY"):\n'
                   ' print("synthetic child stdout",flush=True)\n'
                   ' print("synthetic loader error: test runtime unavailable",file=sys.stderr,flush=True)\n'
                   ' sys.exit(127)\n'
                   'print("FRONTEND_PREPARATION_OK")\n')
    raw.chmod(0o700)
    frontend = home/'bin/R'
    frontend.write_text('#!/bin/sh\nexport SYNTHETIC_LOADER_READY=1\nexec ' + shlex.quote(str(raw)) + ' "$@"\n')
    frontend.chmod(0o700)
    env = os.environ.copy(); env.pop('SYNTHETIC_LOADER_READY', None)
    old = startup.startup_command(home, probe, 'raw', t/'output with spaces', system='Darwin')
    stdout = io.StringIO(); stderr = io.StringIO()
    try:
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            startup.run_child(old, env, t/'diagnostics', 'synthetic_127', cwd=t)
        raise AssertionError('127 child was silently accepted')
    except RuntimeError:
        pass
    z = json.loads((t/'diagnostics/startup_synthetic_127.json').read_text())
    assert z['returncode'] == 127 and z['cwd'] == str(t) and z['probe_exists'] and z['shell'] is False
    assert z['stdout'] in stdout.getvalue() and z['stderr'] in stderr.getvalue()
    assert (t/'diagnostics/startup_synthetic_127.stderr.log').read_text() == z['stderr']
    passed('exit 127 stdout/stderr and command/cwd/probe diagnostics saved and printed before exception (synthetic)')
    command = startup.startup_command(home, probe, 'raw', t/'output with spaces', system='Linux')
    assert command[0] == str(frontend) and '--vanilla' not in command and '--no-environ' not in command
    assert '--file=' + str(probe) in command
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        p = startup.run_child(command, env, t/'diagnostics', 'frontend', cwd=t)
    assert p.returncode == 0 and 'FRONTEND_PREPARATION_OK' in p.stdout
    passed('Linux uses official bin/R frontend, preserves fixture startup and handles spaces (synthetic)')
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            startup.run_child([str(t/'missing executable')], env, t/'diagnostics', 'missing', cwd=t)
            raise AssertionError('Missing executable was accepted')
        except RuntimeError:
            pass
    assert json.loads((t/'diagnostics/startup_missing.json').read_text())['os_error'] == 'FileNotFoundError'
    passed('missing executable is diagnosed with saved stderr before exception')

    # Exercise the actual sbatch EXIT trap in a disposable shell harness.
    root = t/'acceptance job'; (root/'slurm').mkdir(parents=True); (root/'scripts').mkdir()
    (root/'config/dependencies').mkdir(parents=True); (root/'bin').mkdir(); (root/'existing-library').mkdir()
    sentinel = root/'existing-library/installed-package-sentinel'; sentinel.write_bytes(b'preserved existing packages')
    (root/'config/runtime_md5.csv').write_text('path,md5\n')
    (root/'config/dependencies/source_lock.csv').write_text('package,version\n')
    (root/'slurm/configuration.sh').write_text('export TSR_REVISION_ROOT=' + shlex.quote(str(root)) + '\nexport TSR_OUTPUT_ROOT="$TSR_REVISION_ROOT/results"\nexport TSR_R_LIB="$TSR_REVISION_ROOT/existing-library"\ncd "$TSR_REVISION_ROOT"\n')
    (root/'slurm/runtime_environment.sh').write_text('true\n')
    shutil.copy2(ROOT/'scripts/dependency_acceptance_receipt.py', root/'scripts/dependency_acceptance_receipt.py')
    launcher = root/'scripts/tsr_rscript.sh'
    launcher.write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$TSR_REVISION_ROOT/calls.txt"\nexit 0\n'); launcher.chmod(0o700)
    fake_python = root/'bin/python3'
    fake_python.write_text('#!' + sys.executable + '\nimport json,os,sys\nfrom pathlib import Path\n'
                          'if sys.argv[1]=="tests/test_startup_acceptance_regression.py": sys.exit(0)\n'
                          'if sys.argv[1]=="tests/test_r_startup_isolation.py":\n'
                          ' out=Path(os.environ["TSR_ACCEPTANCE_RUN_DIR"])\n'
                          ' (out/"startup_fixture_override.json").write_text(json.dumps(dict(returncode=127,stdout="synthetic stdout",stderr="original synthetic failure")))\n'
                          ' print("original synthetic failure",file=sys.stderr)\n'
                          ' sys.exit(127)\n'
                          'os.execv(' + repr(sys.executable) + ',[' + repr(sys.executable) + ']+sys.argv[1:])\n')
    fake_python.chmod(0o700)
    uname = root/'bin/uname'; uname.write_text('#!/bin/sh\nprintf "Linux\\n"\n'); uname.chmod(0o700)
    env.update(TSR_REVISION_ROOT=str(root), TSR_R_LIB=str(root/'existing-library'),
               SLURM_JOB_ID='SYNTHETIC_ONLY', TSR_ACCEPTANCE_VERIFY_ONLY='1',
               PATH=str(root/'bin') + os.pathsep + env['PATH'])
    p = subprocess.run(['bash', str(ROOT/'slurm/install_and_verify.sbatch')], cwd=root, env=env, text=True, capture_output=True)
    assert p.returncode == 127 and 'DEPENDENCY_ACCEPTANCE_FAILED' in p.stderr
    latest = root/'results/environment/dependency_acceptance_latest.json'
    z = json.loads(latest.read_text())
    assert z['status'] == 'FAILED' and z['stage'] == 'startup_isolation' and z['exit_code'] == 127
    assert z['child_diagnostics'][0]['stderr'] == 'original synthetic failure'
    assert 'test_r_startup_isolation.py' in z['command']
    assert sentinel.read_bytes() == b'preserved existing packages'
    assert 'install_dependencies.R' not in (root/'calls.txt').read_text()
    passed('actual acceptance EXIT trap preserves original stage/code/child error and existing library (synthetic)')
    checker = [sys.executable, str(root/'scripts/dependency_acceptance_receipt.py')]
    # Even an old successful receipt must not mask a newer failed job.
    (root/'results/environment/dependency_acceptance.json').write_text('{}\n')
    p = subprocess.run(checker, env=env, text=True, capture_output=True)
    assert p.returncode == 2 and 'startup_isolation' in p.stderr and 'original synthetic failure' in p.stderr and 'Traceback' not in p.stderr
    passed('receipt checker reports original failure and rejects stale earlier success')
    latest.unlink(); (root/'results/environment/dependency_acceptance.json').unlink()
    p = subprocess.run(checker, env=env, text=True, capture_output=True)
    assert p.returncode == 2 and 'no successful receipt' in p.stderr and 'Traceback' not in p.stderr
    passed('missing receipt produces actionable NOT_READY message without FileNotFoundError traceback')

    accepted = root/'results/acceptance/synthetic_success'; accepted.mkdir(parents=True)
    names = ['library_child_regression.csv', 'startup_isolation_regression.csv', 'minimal_competitors.csv',
             'targeted_checks.csv', 'local_tests.csv', 'workflow_checks.csv', 'scientific_preservation.csv',
             'startup_acceptance_contract.csv']
    for name in names:
        (accepted/name).write_text('check,status\nsynthetic,PASS\n')
    (accepted/'dependency_verification.csv').write_text('loaded,from_isolated_library,version_matches\n' + 'TRUE,TRUE,TRUE\n'*17)
    p = subprocess.run(checker + ['--write', str(accepted)], env=env, text=True, capture_output=True)
    assert p.returncode == 0, p.stderr
    p = subprocess.run(checker, env=env, text=True, capture_output=True)
    assert p.returncode == 0 and 'DEPENDENCY_ACCEPTANCE_RECEIPT_OK' in p.stdout
    passed('all nine passed report groups create and validate a successful receipt (synthetic)')
    p = subprocess.run(checker + ['--begin', str(root/'results/acceptance/new_pending')], env=env, text=True, capture_output=True)
    assert p.returncode == 0
    p = subprocess.run(checker, env=env, text=True, capture_output=True)
    assert p.returncode == 2 and 'DEPENDENCY_ACCEPTANCE_RUNNING' in p.stderr
    passed('a pending new rerun blocks use of an older successful receipt (synthetic)')

out = Path(os.environ.get('TSR_ACCEPTANCE_RUN_DIR', ROOT/'results/acceptance_fix_7016756/local'))
out.mkdir(parents=True, exist_ok=True)
with (out/'startup_acceptance_contract.csv').open('w') as f:
    w = csv.writer(f); w.writerow(['check', 'status']); w.writerows(rows)
print('STARTUP_ACCEPTANCE_REGRESSION_OK ' + str(len(rows)))

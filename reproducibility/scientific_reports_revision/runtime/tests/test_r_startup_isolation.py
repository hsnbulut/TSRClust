"""Fresh-R startup regression; disposable fixtures, no installation/profile edits."""
from pathlib import Path
import csv, hashlib, json, os, platform, shlex, shutil, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
ENV_KEYS = ('R_HOME', 'R_ARCH', 'TSR_R_LIB', 'R_LIBS', 'R_LIBS_USER',
            'R_LIBS_SITE', 'R_ENVIRON', 'R_ENVIRON_USER', 'R_PROFILE',
            'R_PROFILE_USER', 'PATH', 'LD_LIBRARY_PATH', 'R_LD_LIBRARY_PATH',
            'DYLD_FALLBACK_LIBRARY_PATH')


def startup_command(home, probe, mode, output, system=None):
    # Linux must use R's shell frontend (etc/ldpaths), not the internal ELF binary.
    # The Mac Framework frontend points at Current; retain the tested exact binary.
    executable = Path(home) / ('bin/exec/R' if (system or platform.system()) == 'Darwin' else 'bin/R')
    return [str(executable), '--no-save', '--no-restore', '--slave',
            '--file=' + str(Path(probe).resolve()), '--args', mode, str(output)]


def run_child(command, env, output_dir, label, cwd=ROOT, required=True):
    """Persist and print both streams BEFORE reporting an unsuccessful child."""
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    cwd = Path(cwd).resolve()
    probe = next((Path(x.split('=', 1)[1]) for x in command if x.startswith('--file=')), None)
    executable = shutil.which(command[0], path=env.get('PATH')) or command[0]
    diagnostic = dict(label=label, command=command, shell=False, cwd=str(cwd),
                      resolved_executable=executable,
                      executable_exists=Path(executable).is_file(),
                      executable_access=os.access(executable, os.X_OK),
                      probe=str(probe) if probe else None,
                      probe_exists=probe.is_file() if probe else None,
                      probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest() if probe and probe.is_file() else None,
                      environment={k: env.get(k) for k in ENV_KEYS})
    try:
        result = subprocess.run(command, cwd=cwd, env=env, text=True,
                                capture_output=True, check=False)
        diagnostic.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    except OSError as error:
        result = None
        diagnostic.update(returncode=None, stdout='', stderr=str(error), os_error=type(error).__name__)
    prefix = output_dir / ('startup_' + label)
    Path(str(prefix) + '.stdout.log').write_text(diagnostic['stdout'])
    Path(str(prefix) + '.stderr.log').write_text(diagnostic['stderr'])
    Path(str(prefix) + '.json').write_text(json.dumps(diagnostic, indent=2) + '\n')
    print('STARTUP_CHILD_COMMAND ' + ' '.join(shlex.quote(x) for x in command), flush=True)
    print('STARTUP_CHILD_CWD ' + str(cwd), flush=True)
    print('STARTUP_CHILD_RETURN_CODE ' + str(diagnostic['returncode']), flush=True)
    print('STARTUP_CHILD_STDOUT_BEGIN\n' + diagnostic['stdout'] + '\nSTARTUP_CHILD_STDOUT_END', flush=True)
    print('STARTUP_CHILD_STDERR_BEGIN\n' + diagnostic['stderr'] + '\nSTARTUP_CHILD_STDERR_END', file=sys.stderr, flush=True)
    if required and (result is None or result.returncode != 0):
        raise RuntimeError('Startup child failed: ' + label + '; see ' + str(prefix) + '.json')
    return result


def main():
    os.chdir(ROOT)
    lib = Path(os.environ['TSR_R_LIB']).resolve(); assert lib.is_dir()
    out = Path(os.environ.get('TSR_ACCEPTANCE_RUN_DIR', 'results/startup_isolation/regression_7016756')).resolve()
    out.mkdir(parents=True, exist_ok=True)
    home = os.environ.get('TSR_TEST_R_HOME')
    if not home:
        found = run_child(['bash', str(ROOT/'scripts/tsr_rscript.sh'), '-e', 'cat(R.home())'],
                          os.environ.copy(), out, 'discover_home')
        home = found.stdout.strip()
    home = Path(home).resolve()
    probe = (ROOT/'tests/startup_file_probe.R').resolve()
    assert probe.is_file(), 'Missing probe: ' + str(probe)
    original_files = [p for p in lib.rglob('*') if p.is_file()]
    before = {str(p): hashlib.md5(p.read_bytes()).hexdigest() for p in original_files}
    with tempfile.TemporaryDirectory(prefix='tsr-startup-fixture-') as td:
        t = Path(td); legacy = t/'legacy_Rlibs'; legacy.mkdir()
        fixture = t/'Renviron.user'; fixture.write_text('R_LIBS_USER="' + str(legacy) + '"\n')
        fixture_hash = hashlib.sha256(fixture.read_bytes()).hexdigest()
        env = os.environ.copy()
        env.update(R_HOME=str(home), TSR_R_LIB=str(lib), R_LIBS_USER=str(lib), R_LIBS='', R_LIBS_SITE='NULL',
                   R_ENVIRON='/dev/null', R_ENVIRON_USER=str(fixture), R_PROFILE='/dev/null', R_PROFILE_USER='/dev/null')
        if platform.system() == 'Linux':
            # Read-only forensic replay of the former command; diagnostic, NOT a pass criterion.
            old = [str(home/'bin/exec/R'), '--no-save', '--no-restore', '--slave',
                   '--file=' + str(probe), '--args', 'raw', str(out/'legacy_direct_exec.rds')]
            print('LEGACY_INTERNAL_BINARY_DIAGNOSTIC_ONLY; official frontend determines test acceptance', flush=True)
            old_result = run_child(old, env, out, 'legacy_direct_exec', required=False)
            if (old_result is None or old_result.returncode != 0) and shutil.which('ldd'):
                run_child([shutil.which('ldd'), str(home/'bin/exec/R')], env, out, 'legacy_ldd', required=False)
        baseline = out/'fixture_override.rds'
        # Deliberately omit --vanilla/--no-environ: the fixture must really be read.
        run_child(startup_command(home, probe, 'raw', baseline), env, out, 'fixture_override')
        state = list(csv.DictReader(Path(str(baseline)+'.paths.csv').open()))[0]
        assert state['R_version'] == '4.5.2' and Path(state['R_LIBS_USER']).resolve() == legacy.resolve() and Path(state['first_lib']).resolve() == legacy.resolve()
        repaired = out/'fixture_blocked.rds'
        command = ['bash', '-c', 'source "$1"; shift; exec "$@"', '_',
                   str(ROOT/'slurm/r_startup_environment.sh')] + startup_command(home, probe, 'normal', repaired)
        result = run_child(command, env, out, 'fixture_blocked')
        state = list(csv.DictReader(Path(str(repaired)+'.paths.csv').open()))[0]
        assert state['R_version'] == '4.5.2' and Path(state['R_LIBS_USER']).resolve() == lib and Path(state['first_lib']).resolve() == lib
        assert all(state[k] == '/dev/null' for k in ['R_ENVIRON', 'R_ENVIRON_USER', 'R_PROFILE', 'R_PROFILE_USER'])
        packages = list(csv.DictReader(Path(str(repaired)+'.packages.csv').open()))
        assert len(packages) == 17 and all(x['loaded'] == 'TRUE' and x['isolated'] == 'TRUE' and x['version_matches'] == 'TRUE' for x in packages)
        assert fixture_hash == hashlib.sha256(fixture.read_bytes()).hexdigest()
        assert before == {str(p): hashlib.md5(p.read_bytes()).hexdigest() for p in original_files}
        rows = [('real Renviron fixture overrides exported R_LIBS_USER in a fresh R4.5.2 child', 'PASS'),
                ('project policy preserves initial .libPaths in a fresh child without --vanilla', 'PASS'),
                ('all 17 exact pinned packages load from isolated library', 'PASS'),
                ('startup fixture and installed package files unchanged', 'PASS')]
        with (out/'startup_isolation_regression.csv').open('w') as f:
            w = csv.writer(f); w.writerow(['check', 'status']); w.writerows(rows)
        print('R_STARTUP_ISOLATION_OK 4; 17 PACKAGES; NO INSTALLATION')


if __name__ == '__main__':
    main()

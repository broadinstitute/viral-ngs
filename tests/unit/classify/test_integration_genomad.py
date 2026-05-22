# Integration / smoke tests for geNomad
#
# These are intentionally lightweight: a full geNomad run requires the
# multi-GB geNomad database, which is not bundled with the repo. The smoke
# test here verifies the binary is importable and runnable end-to-end
# without exercising classification, which is enough to catch dependency
# regressions (e.g., the conda-forge tensorflow 2.19.1 libtensorflow_cc.so.2
# executable-stack import failure that motivated the tensorflow<2.19 pin).

import shutil
import subprocess

import pytest


pytestmark = pytest.mark.skipif(
    not shutil.which('genomad'),
    reason='geNomad executable is not installed',
)


def test_genomad_version_runs_cleanly():
    '''genomad --version should exit 0 with a non-empty version string.

    Catches the class of regressions where the genomad package installs
    but its python/native deps (tensorflow, libtensorflow_cc, etc.) fail
    to import at runtime.
    '''
    result = subprocess.run(
        ['genomad', '--version'],
        capture_output=True,
        text=True,
        check=True,
    )
    output = (result.stdout or result.stderr).strip()
    assert output, 'genomad --version produced no output'
    assert 'genomad' in output.lower(), (
        f'unexpected genomad --version output: {output!r}'
    )


def test_genomad_end_to_end_help_runs_cleanly():
    '''genomad end-to-end --help should exit 0.

    A successful --help run exercises the same module-import path as a
    real classification (tensorflow, numpy, pyrodigal-gv, etc. all get
    loaded) without needing the multi-GB geNomad database. This is the
    most cost-effective signal that the installed environment is healthy.
    '''
    subprocess.run(
        ['genomad', 'end-to-end', '--help'],
        capture_output=True,
        text=True,
        check=True,
    )
# Testing Patterns

**Analysis Date:** 2026-02-11

## Test Framework

**Runner:**
- pytest 7.0+ (specified in `pyproject.toml`)
- Config: `pyproject.toml` with `[tool.pytest.ini_options]`
- pytest plugins: pytest-xdist (parallel execution), pytest-cov (coverage)

**Assertion Library:**
- unittest.TestCase methods (inherited pattern)
- Standard assertions: `self.assertEqual()`, `self.assertTrue()`, `self.assertFalse()`, `self.assertGreater()`
- Custom assertion helpers in `tests/__init__.py`: `assertEqualContents()`, `assertEqualFasta()`, `assertEqualSamHeaders()`

**Run Commands:**
```bash
# Run all tests (via Docker as per CLAUDE.md)
docker run --rm \
  -v $(pwd):/opt/viral-ngs/source \
  quay.io/broadinstitute/viral-ngs:main-core \
  pytest -rsxX -n auto /opt/viral-ngs/source/tests/unit

# With specific markers
pytest -m "not slow"  # Skip slow tests
pytest -m "slow" --runslow  # Run slow tests with flag

# With coverage
pytest --cov=viral_ngs tests/unit
```

**Pytest Configuration (pyproject.toml):**
```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
addopts = "-v --tb=short"
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
]
```

## Test File Organization

**Location:**
- Tests co-located by module: `tests/unit/core/test_tools_samtools.py` mirrors `src/viral_ngs/core/samtools.py`
- Pattern: `tests/unit/<package>/<test_module>.py` for `src/viral_ngs/<package>/<module>.py`
- Base utilities in `tests/__init__.py` and `tests/conftest.py`

**Naming:**
- Test files: `test_<module>.py`
- Test classes: `Test<FeatureName>` (e.g., `TestToolSamtools`, `TestSamtoolsImport`)
- Test methods: `test_<feature>` (e.g., `test_count_bam`, `test_fasta_index`, `test_import_paired_fastq_basic`)

**Structure:**
```
tests/
├── __init__.py              # Base test class, fixtures, helpers
├── conftest.py              # pytest configuration, global fixtures
└── unit/
    ├── core/
    │   ├── test_tools_samtools.py
    │   ├── test_tools_picard.py
    │   ├── test_util_misc.py
    │   └── ...
    └── ...
```

## Test Structure

**Suite Organization:**
- Tests inherit from `TestCaseWithTmp` (defined in `tests/__init__.py`)
- `TestCaseWithTmp` extends `unittest.TestCase` with temp directory handling
- Decorator: `@pytest.mark.usefixtures('tmpdir_class')`

**Example from test_tools_samtools.py:**
```python
class TestToolSamtools(TestCaseWithTmp):

    def test_count_bam(self):
        sam = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'simple.sam')
        n = viral_ngs.core.samtools.SamtoolsTool().count(sam, ['-S'])
        self.assertEqual(n, 2)

    def test_fasta_index(self):
        orig_ref = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'in.fasta')
        expected_fai = os.path.join(viral_ngs.core.file.get_test_input_path(self), 'in.fasta.fai')
        samtools = viral_ngs.core.samtools.SamtoolsTool()
        for ext in ('.fasta', '.fa'):
            inRef = viral_ngs.core.file.mkstempfname(ext)
            shutil.copyfile(orig_ref, inRef)
            outFai = inRef + '.fai'
            samtools.faidx(inRef)
            self.assertEqualContents(outFai, expected_fai)
```

**Patterns:**

### Setup/Teardown
- No explicit `setUp()`/`tearDown()` (uses class-scope temp directory fixture)
- Fixtures handle cleanup via context managers
- Fixture: `tmpdir_function` (autouse=True) sets tempfile and TMPDIR for each test

### Assertion Patterns
- File comparison: `self.assertEqualContents(file1, file2)` (helper method)
- FASTA comparison: `self.assertEqualFasta(fasta1, fasta2)` (compares sequence IDs and bases)
- BAM/SAM comparison: `assert_equal_bam_reads()` (converts to SAM and compares)
- Header validation: `self.assertEqualSamHeaders(bam1, bam2, other_allowed_values={})`

### Input File Patterns
- Test input files stored by test class: `get_test_input_path(self)` returns class-specific directory
- Usage: `os.path.join(viral_ngs.core.file.get_test_input_path(self), 'simple.sam')`
- Fallback paths without class: `viral_ngs.core.file.get_test_input_path()` (no self parameter)

## Mocking

**Framework:** `monkeypatch` fixture (pytest built-in), custom `monkeypatch_function_result` fixture

**Patterns:**
Custom context manager fixture in `conftest.py`:
```python
@pytest.fixture
def monkeypatch_function_result(monkeypatch):
    """Patches result of a function for specified args"""
    @contextlib.contextmanager
    def _set_function_result(f, *patch_args, **patch_kwargs):
        # Patches f for exactly one argument combination
        # Usage: with monkeypatch_function_result(func, arg1, kwarg1=val1,
        #                                         patch_result=result):
        patch_result = patch_kwargs.pop('patch_result', None)
        patch_exception = patch_kwargs.pop('patch_exception', None)
        patch_module = patch_kwargs.pop('patch_module', inspect.getmodule(f))
        # ... implementation patches function in patch_module
        yield
    return _set_function_result
```

**Usage Example from test_util_misc.py:**
```python
def test_available_cpu_count(monkeypatch_function_result):
    with monkeypatch_function_result(os.path.exists, "/sys/fs/cgroup/cgroup.controllers",
                                     patch_result=True, patch_module=os.path), \
         monkeypatch_function_result(viral_ngs.core.file.slurp_file,
                                     '/sys/fs/cgroup/cpu.max',
                                     patch_result="100000 100000"):
        # Test code using mocked functions
```

**What to Mock:**
- System files and environment: `/sys/fs/cgroup/`, `/proc/self/status`
- External commands when testing logic (not execution)
- File I/O when testing file-processing logic separately

**What NOT to Mock:**
- Tool installation/execution (use real tools in Docker container)
- File system operations (use real temp directories)
- Data processing functions (test with real data)

## Fixtures and Factories

**Test Data:**
- Input files in `tests/*/TestToolName/` directories (organized by test class)
- Files accessed via `get_test_input_path(self)` which maps class name to directory
- Test data preserved if `VIRAL_NGS_TMP_DIRKEEP` env var is set (for debugging)

**Locations:**
- `tests/` - global test utilities and conftest
- `tests/unit/` - unit test files

**Key Fixtures from conftest.py:**

### Temp Directory Fixtures (multi-scope)
```python
@pytest.fixture(scope='session')
def tmpdir_session(request, tmpdir_factory):
    """Create a session-scope temporary directory."""
    with _tmpdir_aux(str(tmpdir_factory.getbasetemp()),
                     'session', id(request.session)) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='module')
def tmpdir_module(request, tmpdir_session):
    """Create a module-scope temporary directory."""
    with _tmpdir_aux(tmpdir_session, 'module', request.module.__name__) as tmpdir:
        yield tmpdir

@pytest.fixture(scope='class')
def tmpdir_class(request, tmpdir_module):
    """Create a class-scope temporary directory."""
    with _tmpdir_aux(tmpdir_module, 'class',
                     request.cls.__name__ if request.cls else '__noclass__') as tmpdir:
        yield tmpdir

@pytest.fixture(autouse=True)
def tmpdir_function(request, tmpdir_class, monkeypatch):
    """Create a temporary directory and set it as TMPDIR for tempfile module."""
    with _tmpdir_aux(tmpdir_class, 'node', request.node.name) as tmpdir:
        monkeypatch.setattr(tempfile, 'tempdir', tmpdir)
        monkeypatch.setenv('TMPDIR', tmpdir)
        yield tmpdir
```

## Coverage

**Requirements:** None enforced (no minimum threshold in config)

**View Coverage:**
```bash
pytest --cov=viral_ngs --cov-report=html tests/unit
# Then open htmlcov/index.html
```

**.coveragerc Configuration:**
```ini
[run]
branch = True
omit = tools/conda-tools/*
disable_warnings = module-not-imported
relative_files = True

[report]
exclude_lines =
    pragma: no cover
    raise NotImplementedError
    def __repr__
    if self\.debug
    if 0:
    if __name__ == .__main__.:
```

## Test Types

**Unit Tests:**
- Scope: Individual functions and Tool methods
- Location: `tests/unit/core/test_*.py`
- Approach: Mock external dependencies, test with real test data files
- Example: `test_count_bam()` tests counting logic with real BAM file

**Integration Tests:**
- Not explicitly separated (would be in `tests/integration/`)
- Tool tests like `test_import_fastq()` verify end-to-end FASTQ->BAM conversion
- Tests call real tools (in Docker container) with test data

**E2E Tests:**
- Not found in codebase
- Would likely be shell scripts in `tests/` directory

## Common Patterns

**Async Testing:**
- Not used (synchronous subprocess execution)

**Error Testing:**
```python
def testBasicRunFailAndCatch(self):
    self.assertRaises(subprocess.CalledProcessError,
        viral_ngs.core.misc.run_and_print, ['cat', '/notdev/notnull'],
        silent=False, buffered=False, check=True)
```

**Parametrized Tests:**
From `test_tools.py`:
```python
@pytest.fixture(params=viral_ngs.core.all_tool_classes())
def tool_class(request):
    return request.param

def test_tool_install(tool_class):
    if IS_ARM and tool_class.__name__ in X86_ONLY_TOOLS:
        pytest.skip(f"{tool_class.__name__} requires x86-only package")
    t = tool_class()
    t.install()
    assert t.is_installed()
```

**Skipping Tests:**
- Platform-specific: `pytest.skip()` for ARM-only restrictions
- Slow tests: `@pytest.mark.slow` with `--runslow` flag to include

**Test Helpers in tests/__init__.py:**
- `assert_equal_contents(testCase, f1, f2)` - File comparison
- `assert_equal_bam_reads(testCase, bam1, bam2)` - BAM comparison (converts to SAM)
- `assert_md5_equal_to_line_in_file(testCase, file, checksum_file)` - Checksum validation
- `assert_valid_feature_table(testCase, tbl_file, fasta_file, temp_dir)` - Feature table validation (placeholder)
- `assert_equal_bam_reads()` - Helper with debug output for BAM comparison failures

## Slow Test Marker

**Usage:**
```python
@pytest.mark.slow
def test_slow_operation(self):
    # Long-running test
    pass
```

**Skip by default:**
- Tests marked `@pytest.mark.slow` are skipped unless `--runslow` flag is passed
- Logic in `conftest.py` `pytest_collection_modifyitems()` hook

## Test Discovery

**Default:**
- pytest finds `test_*.py` files in `tests/` directory (testpaths in pyproject.toml)
- Excludes files with executable bit set (asserted in `tests/__init__.py`)

---

*Testing analysis: 2026-02-11*

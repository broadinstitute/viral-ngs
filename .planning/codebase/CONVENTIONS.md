# Coding Conventions

**Analysis Date:** 2026-02-11

## Naming Patterns

**Files:**
- Module files use snake_case (e.g., `samtools.py`, `picard.py`, `misc.py`)
- Test files follow pattern `test_<module>.py` (e.g., `test_tools_samtools.py`, `test_util_misc.py`)
- Location: `src/viral_ngs/` for source modules, `tests/unit/core/` for tests

**Functions:**
- Use snake_case for all functions (e.g., `unambig_count`, `run_and_print`, `available_cpu_count`)
- No underscores for "private" functions (convention is snake_case, no leading underscore distinction)
- Method names in classes follow same snake_case pattern (e.g., `install_and_get_path`, `count_bam`)

**Variables:**
- Use snake_case (e.g., `inFastq1`, `outBam`, `tool_cmd`)
- Mix of snake_case and camelCase in codebase - follow existing patterns in module
- Parameters use snake_case: `sample_name`, `library_name`, `platform`, `sequencing_center`

**Types/Classes:**
- Use PascalCase for classes (e.g., `SamtoolsTool`, `Tool`, `FeatureSorter`, `QCError`)
- Exception classes inherit from standard exceptions (e.g., `QCError(RuntimeError)`, `InvalidBamHeaderError(ValueError)`)

**Module-level:**
- `__author__` string with email (e.g., `"dpark@broadinstitute.org"`)
- `log = logging.getLogger(__name__)` is standard module-level logger
- Constants use UPPERCASE (e.g., `MAX_INT32`, `TOOL_NAME = 'samtools'`)

## Code Style

**Formatting:**
- No explicit formatter configured (.ruff.toml, .black.toml, .flake8 not found)
- Code appears to follow PEP 8 conventions loosely
- Imports organized in groups with comments (built-ins, third-party, intra-project)

**Linting:**
- Pylint pragma comments used locally: `# pylint: disable=W0221` (e.g., in `samtools.py` method override)
- No project-wide linting enforcement detected

**Line Length:**
- Multi-line function parameters aligned vertically (e.g., `collapse_dup_strs_to_str_or_md5` uses aligned parameter formatting)
- Some lines exceed 100 characters, no strict enforcement visible

## Import Organization

**Order:**
1. Built-in modules (sys, os, re, logging, subprocess, tempfile, etc.)
2. Third-party packages (Bio, pysam, numpy, yaml, etc.)
3. Intra-project imports (from . import, from viral_ngs import)

**Path Aliases:**
- Import aliases used for clarity: `from . import file as util_file, misc as util_misc`
- Pattern: `from . import <module>` for same-package imports
- Full module paths for cross-package: `import viral_ngs.core.samtools`
- Package import in `__init__.py`: `from . import core` (makes `viral_ngs.core.*` accessible)

**Example from `samtools.py`:**
```python
import logging
import shutil
import os
import re
import os.path
import subprocess
import tempfile
import contextlib
from collections import OrderedDict
from decimal import *

import pysam

from . import file as util_file, misc as util_misc
from . import Tool, PrexistingUnixCommand
```

## Error Handling

**Patterns:**
- Raise exceptions with descriptive messages (e.g., `raise exc(message)`)
- Custom exceptions defined in `src/viral_ngs/core/errors.py`
- Use `chk()` helper function for assertions: `chk(condition, message='...', exc=RuntimeError)`
- Subprocess errors caught and re-raised as `subprocess.CalledProcessError`

**Custom Exceptions:**
- `QCError(RuntimeError)` - for QC step failures
- `InvalidBamHeaderError(ValueError)` - for malformed BAM headers
- `StringNotFoundException(Exception)` - when substring not found in files

**Error Propagation:**
- Methods call `subprocess.check_call()` or `subprocess.check_output()` to automatically raise on non-zero exit
- File operations use exception context managers: `with open(...) as f:` with try/finally for cleanup

## Logging

**Framework:** `logging` module (standard library)

**Patterns:**
- Module-level logger: `log = logging.getLogger(__name__)`
- Log levels used: debug, info, warning, error
- Usage: `log.debug()`, `log.warning()`, `log.error()`
- Debug logging for command execution (e.g., `log.debug(' '.join(tool_cmd))`)
- Warning for parameter renames: `log.warning('Config param {} has been renamed to {}')`
- Multiline context logging with perf metrics: `log.debug(f"PERF: filter_time={elapsed:.2f}s ...")`

**Configuration:**
- Logging setup in `src/viral_ngs/core/cmd.py` via `setup_logger(log_level)`
- Formatter: `"%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s"`

## Comments

**When to Comment:**
- Inline comments explain non-obvious logic (e.g., "cgroup v2 detected", "Header lines (start with @) are passed through unchanged")
- Docstrings for all public functions and classes
- TODOs mark incomplete work (e.g., "TO DO: much of this stuff can be all eliminated by using pysam instead")

**JSDoc/TSDoc:**
- Uses triple-quoted docstrings (not per-parameter style)
- Docstring format: narrative description of behavior
- Example from `collapse_dup_strs_to_str_or_md5`:
  ```python
  """
  Collapse multiple string values into one string value

  Given a list of values (ex. from a column in a duplicated group of >=2 rows):
    1) If all values are empty (""), return "".
    2) Otherwise, if there's exactly 1 unique value (non-empty or empty), return that value + (suffix, optionally)
    3) Otherwise (>=2 distinct values), join with delimiter...
  """
  ```

**Comments Style:**
- Use single quotes for inline docstrings: `'''Brief description.'''`
- Use triple double quotes for full docstrings: `"""Full description."""`
- Comments start with descriptive verb: "Return", "Check", "Create", "Validate"

## Function Design

**Size:**
- Functions vary from 1-line to 50+ lines
- Large functions are complex operations (e.g., `collapse_dup_strs_to_str_or_md5` ~50 lines for intricate logic)
- Tool methods typically 5-20 lines (delegates to subprocess)

**Parameters:**
- Use keyword arguments for tool parameters: `samtools.import_fastq(inFastq1, inFastq2, outBam, sample_name='TestSample', library_name='Alexandria', ...)`
- Optional parameters with defaults: `def run_and_print(..., silent=False, buffered=False, check=False, loglevel=None)`
- Lists passed as args (e.g., `args` parameter for tool commands)

**Return Values:**
- Most tool methods return None (side effects on files)
- Utility functions return computed values: `def unambig_count(seq): ... return sum(...)`
- Context managers use `yield` (e.g., `@contextlib.contextmanager`)
- Methods returning file handle: `def bam2fq_pipe(self, inBam, threads=None): ... return p`

**Type Hints:**
- No type hints used in codebase (Python 3.10+, but conventions don't include them)

## Module Design

**Exports:**
- Modules export all public functions/classes (no `__all__` found)
- Tool classes inherit from `Tool` or `PrexistingUnixCommand` base classes
- Import pattern: `from viral_ngs.core.samtools import SamtoolsTool` or `import viral_ngs.core.samtools`

**Barrel Files:**
- `__init__.py` in `src/viral_ngs/` imports core: `from . import core`
- Allows `import viral_ngs.core` access to submodules
- No other barrel patterns (index exports) detected

**Module Structure:**
- Single responsibility: `samtools.py` wraps samtools, `picard.py` wraps picard
- Tool wrappers inherit from `Tool` class (defined in `src/viral_ngs/core/__init__.py`)
- Utility modules are flat (functions only): `misc.py`, `file.py`

## Docstring Examples

From `misc.py`:
```python
def unique(items):
    ''' Return unique items in the same order as seen in the input. '''
    seen = set()
    for i in items:
        if i not in seen:
            seen.add(i)
            yield i
```

From `file.py`:
```python
def get_project_path():
    '''Return the absolute path of the top-level project, assumed to be the
       parent of the directory containing this script.

       When package is pip-installed (e.g., in Docker), the relative path from
       __file__ won't work. In that case, check for VIRAL_NGS_SOURCE_DIR env var
       or fall back to cwd if it contains a tests/ directory.
    '''
```

---

*Convention analysis: 2026-02-11*

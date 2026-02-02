#!/usr/bin/env python3
"""
Rewrite imports from flat viral-core style to viral_ngs package style.

Usage:
    python scripts/rewrite_imports.py src/viral_ngs/
    python scripts/rewrite_imports.py tests/
"""

import re
import sys
from pathlib import Path

# Module names that need viral_ngs prefix
TOP_LEVEL_MODULES = {
    'broad_utils', 'errors', 'file_utils', 'illumina',
    'priorities', 'read_utils', 'reports'
}

SUBPACKAGES = {'util', 'tools'}


def rewrite_imports(content: str) -> str:
    """Rewrite imports to use viral_ngs package prefix."""
    lines = content.split('\n')
    result = []

    for line in lines:
        # Skip comments and empty lines
        stripped = line.lstrip()
        if stripped.startswith('#') or not stripped:
            result.append(line)
            continue

        original_line = line

        # Pattern: import util.X -> from viral_ngs.util import X
        match = re.match(r'^(\s*)import (util|tools)\.(\w+)\s*$', line)
        if match:
            indent, pkg, module = match.groups()
            result.append(f'{indent}from viral_ngs.{pkg} import {module}')
            continue

        # Pattern: import util -> from viral_ngs import util
        match = re.match(r'^(\s*)import (util|tools)\s*$', line)
        if match:
            indent, pkg = match.groups()
            result.append(f'{indent}from viral_ngs import {pkg}')
            continue

        # Pattern: from util.X import Y -> from viral_ngs.util.X import Y
        match = re.match(r'^(\s*)from (util|tools)\.(\S+) import (.+)$', line)
        if match:
            indent, pkg, module, imports = match.groups()
            result.append(f'{indent}from viral_ngs.{pkg}.{module} import {imports}')
            continue

        # Pattern: from util import X -> from viral_ngs.util import X
        match = re.match(r'^(\s*)from (util|tools) import (.+)$', line)
        if match:
            indent, pkg, imports = match.groups()
            result.append(f'{indent}from viral_ngs.{pkg} import {imports}')
            continue

        # Pattern: import read_utils -> from viral_ngs import read_utils
        match = re.match(r'^(\s*)import (' + '|'.join(TOP_LEVEL_MODULES) + r')\s*$', line)
        if match:
            indent, module = match.groups()
            result.append(f'{indent}from viral_ngs import {module}')
            continue

        # Pattern: from errors import * -> from viral_ngs.errors import *
        match = re.match(r'^(\s*)from (' + '|'.join(TOP_LEVEL_MODULES) + r') import (.+)$', line)
        if match:
            indent, module, imports = match.groups()
            result.append(f'{indent}from viral_ngs.{module} import {imports}')
            continue

        # No match, keep original
        result.append(line)

    return '\n'.join(result)


def process_file(filepath: Path) -> bool:
    """Process a single file. Returns True if changes were made."""
    content = filepath.read_text()
    new_content = rewrite_imports(content)

    if content != new_content:
        filepath.write_text(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python rewrite_imports.py <directory>...")
        sys.exit(1)

    changed_count = 0
    for dir_path in sys.argv[1:]:
        root = Path(dir_path)
        for pyfile in root.rglob('*.py'):
            if process_file(pyfile):
                print(f"Rewrote: {pyfile}")
                changed_count += 1

    print(f"\nTotal files modified: {changed_count}")


if __name__ == '__main__':
    main()

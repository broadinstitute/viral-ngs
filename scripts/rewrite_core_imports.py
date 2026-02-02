#!/usr/bin/env python3
"""Rewrite imports in core/ to use relative imports within the package."""

import os
import re
import sys

def rewrite_imports(content, filename):
    """Rewrite import statements to use relative imports within core/."""
    lines = content.split('\n')
    new_lines = []

    for line in lines:
        original = line

        # from viral_ngs.tools import X -> from . import X
        match = re.match(r'^(\s*)from viral_ngs\.tools import (.+)$', line)
        if match:
            indent, imports = match.groups()
            line = f'{indent}from . import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.util import X -> from . import X
        match = re.match(r'^(\s*)from viral_ngs\.util import (.+)$', line)
        if match:
            indent, imports = match.groups()
            line = f'{indent}from . import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.util.X import Y -> from .X import Y
        match = re.match(r'^(\s*)from viral_ngs\.util\.(\w+) import (.+)$', line)
        if match:
            indent, module, imports = match.groups()
            line = f'{indent}from .{module} import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.tools.X import Y -> from .X import Y
        match = re.match(r'^(\s*)from viral_ngs\.tools\.(\w+) import (.+)$', line)
        if match:
            indent, module, imports = match.groups()
            line = f'{indent}from .{module} import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # import viral_ngs.tools.X -> from . import X
        match = re.match(r'^(\s*)import viral_ngs\.tools\.(\w+)$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}from . import {module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # import viral_ngs.util.X -> from . import X
        match = re.match(r'^(\s*)import viral_ngs\.util\.(\w+)$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}from . import {module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs import util -> (keep as comment, mark for manual fix)
        # This one is tricky - we'll replace it with imports of specific submodules
        match = re.match(r'^(\s*)from viral_ngs import util\s*$', line)
        if match:
            indent = match.groups()[0]
            # Replace with comment - the actual util usage needs analysis
            line = f'{indent}from . import file as util_file, misc as util_misc  # was: from viral_ngs import util'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs import tools -> replace with specific imports
        match = re.match(r'^(\s*)from viral_ngs import tools\s*$', line)
        if match:
            indent = match.groups()[0]
            line = f'{indent}from . import samtools, picard  # was: from viral_ngs import tools'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs import X (where X is a module in core)
        # Handle illumina, read_utils, etc.
        match = re.match(r'^(\s*)from viral_ngs import (read_utils|illumina|reports|file_utils|broad_utils|priorities|errors)\s*$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}from . import {module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        new_lines.append(line)

    return '\n'.join(new_lines)


def process_file(filepath):
    """Process a single Python file."""
    with open(filepath, 'r') as f:
        content = f.read()

    new_content = rewrite_imports(content, os.path.basename(filepath))

    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python rewrite_core_imports.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]
    modified = 0

    for filename in os.listdir(directory):
        if filename.endswith('.py') and filename != '__init__.py':
            filepath = os.path.join(directory, filename)
            if process_file(filepath):
                modified += 1

    print(f"\nModified {modified} files")


if __name__ == '__main__':
    main()

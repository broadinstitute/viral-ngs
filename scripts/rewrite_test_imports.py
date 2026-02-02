#!/usr/bin/env python3
"""Rewrite imports in tests/ to use viral_ngs.core.* paths."""

import os
import re
import sys


def rewrite_imports(content, filename):
    """Rewrite import statements to use viral_ngs.core.*."""
    lines = content.split('\n')
    new_lines = []

    for line in lines:
        original = line

        # from viral_ngs import tools -> import viral_ngs.core
        match = re.match(r'^(\s*)from viral_ngs import tools\s*$', line)
        if match:
            indent = match.groups()[0]
            line = f'{indent}import viral_ngs.core'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs import util -> import viral_ngs.core
        match = re.match(r'^(\s*)from viral_ngs import util\s*$', line)
        if match:
            indent = match.groups()[0]
            line = f'{indent}import viral_ngs.core'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.tools import X -> import viral_ngs.core.X
        match = re.match(r'^(\s*)from viral_ngs\.tools import (\w+)$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}import viral_ngs.core.{module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.tools import X, Y, Z -> multiple import lines
        match = re.match(r'^(\s*)from viral_ngs\.tools import (.+)$', line)
        if match and ',' in match.groups()[1]:
            indent, modules = match.groups()
            module_list = [m.strip() for m in modules.split(',')]
            line = '\n'.join([f'{indent}import viral_ngs.core.{m}' for m in module_list])
            print(f"  {filename}: '{original.strip()}' -> multiple imports")

        # from viral_ngs.util import X -> import viral_ngs.core.X
        match = re.match(r'^(\s*)from viral_ngs\.util import (\w+)$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}import viral_ngs.core.{module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.util.file import X -> from viral_ngs.core.file import X
        match = re.match(r'^(\s*)from viral_ngs\.util\.(\w+) import (.+)$', line)
        if match:
            indent, module, imports = match.groups()
            line = f'{indent}from viral_ngs.core.{module} import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.tools.X import Y -> from viral_ngs.core.X import Y
        match = re.match(r'^(\s*)from viral_ngs\.tools\.(\w+) import (.+)$', line)
        if match:
            indent, module, imports = match.groups()
            line = f'{indent}from viral_ngs.core.{module} import {imports}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs import X (where X is a module in core)
        match = re.match(r'^(\s*)from viral_ngs import (read_utils|illumina|reports|file_utils|broad_utils|priorities|errors)\s*$', line)
        if match:
            indent, module = match.groups()
            line = f'{indent}import viral_ngs.core.{module}'
            print(f"  {filename}: '{original.strip()}' -> '{line.strip()}'")

        # from viral_ngs.util.misc import X -> from viral_ngs.core.misc import X
        match = re.match(r'^(\s*)from viral_ngs\.util\.misc import (.+)$', line)
        if match:
            indent, imports = match.groups()
            line = f'{indent}from viral_ngs.core.misc import {imports}'
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
        print("Usage: python rewrite_test_imports.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]
    modified = 0

    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('.py'):
                filepath = os.path.join(root, filename)
                if process_file(filepath):
                    modified += 1

    print(f"\nModified {modified} files")


if __name__ == '__main__':
    main()

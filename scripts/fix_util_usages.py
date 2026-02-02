#!/usr/bin/env python3
"""Fix usages of util.X to util_X after import changes."""

import os
import re
import sys


def fix_usages(content, filename):
    """Fix usages of util.X to util_X."""
    original = content

    # Replace util.file. with util_file.
    content = re.sub(r'\butil\.file\.', 'util_file.', content)
    content = re.sub(r'\butil\.file\b(?!\.)', 'util_file', content)

    # Replace util.misc. with util_misc.
    content = re.sub(r'\butil\.misc\.', 'util_misc.', content)
    content = re.sub(r'\butil\.misc\b(?!\.)', 'util_misc', content)

    # Replace util.cmd. with util_cmd.
    content = re.sub(r'\butil\.cmd\.', 'util_cmd.', content)
    content = re.sub(r'\butil\.cmd\b(?!\.)', 'util_cmd', content)

    # Replace util.version. with util_version.
    content = re.sub(r'\butil\.version\.', 'util_version.', content)
    content = re.sub(r'\butil\.version\b(?!\.)', 'util_version', content)

    # Replace util.stats. with util_stats.
    content = re.sub(r'\butil\.stats\.', 'util_stats.', content)
    content = re.sub(r'\butil\.stats\b(?!\.)', 'util_stats', content)

    # Replace tools.X. with just X. for known tools
    for tool in ['samtools', 'picard', 'bwa', 'bbmap', 'minimap2', 'novoalign',
                 'gatk', 'cdhit', 'fastqc', 'prinseq', 'sambamba', 'splitcode',
                 'trimmomatic', 'mvicuna']:
        content = re.sub(rf'\btools\.{tool}\.', f'{tool}.', content)
        content = re.sub(rf'\btools\.{tool}\b(?!\.)', tool, content)

    # Replace tools.Tool, tools.InstallMethod, tools.PrexistingUnixCommand
    # with just Tool, InstallMethod, PrexistingUnixCommand
    for cls in ['Tool', 'InstallMethod', 'PrexistingUnixCommand',
                'installed_tools', 'all_tool_classes', 'get_tool_by_name',
                'skip_install_test', 'is_osx', 'iter_leaf_subclasses']:
        content = re.sub(rf'\btools\.{cls}\b', cls, content)

    if content != original:
        print(f"  Fixed usages in {filename}")

    return content


def process_file(filepath):
    """Process a single Python file."""
    with open(filepath, 'r') as f:
        content = f.read()

    new_content = fix_usages(content, os.path.basename(filepath))

    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_util_usages.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]
    modified = 0

    for filename in os.listdir(directory):
        if filename.endswith('.py'):
            filepath = os.path.join(directory, filename)
            if process_file(filepath):
                modified += 1

    print(f"\nModified {modified} files")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Fix usages in tests after import rewrites."""

import os
import re
import sys


def fix_usages(content, filename):
    """Fix usages after import changes."""
    original = content

    # Replace util.X with viral_ngs.core.X for known modules
    for module in ['file', 'misc', 'cmd', 'version', 'stats']:
        content = re.sub(rf'\butil\.{module}\.', f'viral_ngs.core.{module}.', content)
        content = re.sub(rf'\butil\.{module}\b(?!\.)', f'viral_ngs.core.{module}', content)

    # Replace tools.X with viral_ngs.core.X for known tool wrappers
    for tool in ['samtools', 'picard', 'bwa', 'bbmap', 'minimap2', 'novoalign',
                 'gatk', 'cdhit', 'fastqc', 'prinseq', 'sambamba', 'splitcode',
                 'trimmomatic', 'mvicuna']:
        content = re.sub(rf'\btools\.{tool}\.', f'viral_ngs.core.{tool}.', content)
        content = re.sub(rf'\btools\.{tool}\b(?!\.)', f'viral_ngs.core.{tool}', content)

    # Replace tools.Tool, tools.InstallMethod, etc.
    for cls in ['Tool', 'InstallMethod', 'PrexistingUnixCommand',
                'installed_tools', 'all_tool_classes', 'get_tool_by_name',
                'skip_install_test', 'is_osx', 'iter_leaf_subclasses']:
        content = re.sub(rf'\btools\.{cls}\b', f'viral_ngs.core.{cls}', content)

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
        print("Usage: python fix_test_usages.py <directory>")
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

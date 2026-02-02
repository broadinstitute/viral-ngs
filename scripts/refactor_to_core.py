#!/usr/bin/env python3
"""
Refactor imports from viral_ngs.tools.* and viral_ngs.util.* to viral_ngs.core.*
Uses 'import x.y.z' style as preferred.
"""

import re
import sys
from pathlib import Path

# Mapping of old module locations to new
# All tools and util modules are now directly under core
TOOLS_MODULES = [
    'bbmap', 'bwa', 'cdhit', 'fastqc', 'gatk', 'minimap2', 'mvicuna',
    'novoalign', 'picard', 'prinseq', 'sambamba', 'samtools', 'splitcode', 'trimmomatic'
]
UTIL_MODULES = ['cmd', 'file', 'misc', 'stats', 'version', 'illumina_indices']
TOP_MODULES = ['broad_utils', 'errors', 'file_utils', 'illumina', 'priorities', 'read_utils', 'reports']

# The tools __init__.py defines Tool, InstallMethod, PrexistingUnixCommand classes
# These need special handling - they'll be in core/__init__.py


def refactor_imports(content: str, filepath: Path) -> str:
    """Refactor imports to use viral_ngs.core.*"""

    # Pattern: from viral_ngs.tools import X -> import viral_ngs.core.X
    for mod in TOOLS_MODULES:
        # from viral_ngs.tools import X
        pattern = rf'^(\s*)from viral_ngs\.tools import {mod}\s*$'
        replacement = rf'\1import viral_ngs.core.{mod}'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

        # from viral_ngs.tools.X import Y -> from viral_ngs.core.X import Y
        pattern = rf'^(\s*)from viral_ngs\.tools\.{mod} import (.+)$'
        replacement = rf'\1from viral_ngs.core.{mod} import \2'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

    # Pattern: from viral_ngs.util import X -> import viral_ngs.core.X
    for mod in UTIL_MODULES:
        # from viral_ngs.util import X
        pattern = rf'^(\s*)from viral_ngs\.util import {mod}\s*$'
        replacement = rf'\1import viral_ngs.core.{mod}'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

        # from viral_ngs.util.X import Y -> from viral_ngs.core.X import Y
        pattern = rf'^(\s*)from viral_ngs\.util\.{mod} import (.+)$'
        replacement = rf'\1from viral_ngs.core.{mod} import \2'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

    # Pattern: from viral_ngs import X (top-level modules) -> import viral_ngs.core.X
    for mod in TOP_MODULES:
        pattern = rf'^(\s*)from viral_ngs import {mod}\s*$'
        replacement = rf'\1import viral_ngs.core.{mod}'
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

    # Pattern: from viral_ngs import tools -> import viral_ngs.core
    content = re.sub(
        r'^(\s*)from viral_ngs import tools\s*$',
        r'\1import viral_ngs.core',
        content, flags=re.MULTILINE
    )

    # Pattern: from viral_ngs import util -> import viral_ngs.core
    content = re.sub(
        r'^(\s*)from viral_ngs import util\s*$',
        r'\1import viral_ngs.core',
        content, flags=re.MULTILINE
    )

    # Now update usages to match new imports
    # tools.X -> viral_ngs.core.X (when tools was imported as package)
    content = re.sub(r'\btools\.Tool\b', 'viral_ngs.core.Tool', content)
    content = re.sub(r'\btools\.PrexistingUnixCommand\b', 'viral_ngs.core.PrexistingUnixCommand', content)
    content = re.sub(r'\btools\.CondaPackage\b', 'viral_ngs.core.CondaPackage', content)
    content = re.sub(r'\btools\.DownloadPackage\b', 'viral_ngs.core.DownloadPackage', content)

    # util.X -> viral_ngs.core.X (when util was imported as package)
    for mod in UTIL_MODULES:
        content = re.sub(rf'\butil\.{mod}\b', f'viral_ngs.core.{mod}', content)

    # For the specific pattern: picard.X -> viral_ngs.core.picard.X
    # (when module was imported with 'from viral_ngs.tools import picard')
    # This is complex because we changed 'from ... import picard' to 'import viral_ngs.core.picard'
    # So we need to update usages from 'picard.X' to 'viral_ngs.core.picard.X'
    # But this is tricky... let's be more surgical

    return content


def process_file(filepath: Path) -> bool:
    """Process a single file. Returns True if changes were made."""
    content = filepath.read_text()
    new_content = refactor_imports(content, filepath)

    if content != new_content:
        filepath.write_text(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python refactor_to_core.py <directory>...")
        sys.exit(1)

    changed_count = 0
    for dir_path in sys.argv[1:]:
        root = Path(dir_path)
        for pyfile in root.rglob('*.py'):
            if process_file(pyfile):
                print(f"Refactored: {pyfile}")
                changed_count += 1

    print(f"\nTotal files modified: {changed_count}")


if __name__ == '__main__':
    main()

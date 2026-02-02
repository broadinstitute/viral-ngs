#!/usr/bin/env python3
"""
Fix util imports for files that use util.X pattern.

Changes:
  from viral_ngs.util import cmd -> from viral_ngs import util
  from viral_ngs.util import file -> (removed, covered by above)
  from viral_ngs.util import misc -> (removed, covered by above)

Only applies to files that use util.X references.
"""

import re
import sys
from pathlib import Path


def fix_util_imports(content: str) -> str:
    """Fix util imports to use the package import pattern."""
    lines = content.split('\n')

    # Check if this file uses util.X pattern
    uses_util_pattern = bool(re.search(r'\butil\.(cmd|file|misc|stats|version)\b', content))

    if not uses_util_pattern:
        return content

    result = []
    util_import_added = False
    util_submodule_imports = []

    for line in lines:
        # Collect util submodule imports and replace with single util import
        match = re.match(r'^(\s*)from viral_ngs\.util import (\w+)\s*$', line)
        if match:
            indent = match.group(1)
            module = match.group(2)
            util_submodule_imports.append(module)
            if not util_import_added:
                result.append(f'{indent}from viral_ngs import util')
                util_import_added = True
            # Skip this line (don't add the from viral_ngs.util import X)
            continue

        result.append(line)

    return '\n'.join(result)


def process_file(filepath: Path) -> bool:
    """Process a single file. Returns True if changes were made."""
    content = filepath.read_text()
    new_content = fix_util_imports(content)

    if content != new_content:
        filepath.write_text(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_util_imports.py <directory>...")
        sys.exit(1)

    changed_count = 0
    for dir_path in sys.argv[1:]:
        root = Path(dir_path)
        for pyfile in root.rglob('*.py'):
            if process_file(pyfile):
                print(f"Fixed: {pyfile}")
                changed_count += 1

    print(f"\nTotal files modified: {changed_count}")


if __name__ == '__main__':
    main()

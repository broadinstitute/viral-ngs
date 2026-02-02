#!/usr/bin/env python3
"""
Fix module usage patterns to match imports.

Changes patterns like:
  tools.picard.X -> picard.X (when picard is imported)
  util.file.X -> file.X (when file is imported as util_file or file)
"""

import re
import sys
from pathlib import Path


def fix_usages(content: str) -> str:
    """Fix module usage patterns to match imports."""

    # For files that import 'from viral_ngs.tools import X':
    # Change tools.X.something to X.something

    # Find what's imported from tools
    tools_imports = re.findall(r'from viral_ngs\.tools import (\w+)', content)
    for module in tools_imports:
        # Replace tools.module.X with module.X
        pattern = rf'\btools\.{module}\.'
        replacement = f'{module}.'
        content = re.sub(pattern, replacement, content)

    # Also remove 'from viral_ngs import tools' lines that are no longer needed
    # if all tools.X usages are replaced
    if tools_imports and not re.search(r'\btools\.\w', content):
        content = re.sub(r'^from viral_ngs import tools\n', '', content, flags=re.MULTILINE)

    return content


def process_file(filepath: Path) -> bool:
    """Process a single file. Returns True if changes were made."""
    content = filepath.read_text()
    new_content = fix_usages(content)

    if content != new_content:
        filepath.write_text(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_module_usages.py <directory>...")
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

#!/usr/bin/env python3
"""Add imports for Tool and PrexistingUnixCommand where needed."""

import os
import re
import sys


def needs_import(content, name):
    """Check if the file uses a name but doesn't import it."""
    # Check if name is used (as class base, function call, etc.)
    pattern = rf'\b{name}\b'
    if not re.search(pattern, content):
        return False

    # Check if it's already imported
    import_pattern = rf'from \. import .*\b{name}\b'
    if re.search(import_pattern, content):
        return False

    return True


def add_imports(content, filename):
    """Add missing imports for base classes."""
    imports_needed = []

    for name in ['Tool', 'InstallMethod', 'PrexistingUnixCommand']:
        if needs_import(content, name):
            imports_needed.append(name)

    if not imports_needed:
        return content

    import_line = f"from . import {', '.join(imports_needed)}\n"

    # Find the last import line and add after it
    lines = content.split('\n')
    last_import_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('import ') or line.startswith('from '):
            last_import_idx = i

    # Insert after last import
    lines.insert(last_import_idx + 1, import_line.rstrip())

    print(f"  {filename}: Added import for {', '.join(imports_needed)}")
    return '\n'.join(lines)


def process_file(filepath):
    """Process a single Python file."""
    with open(filepath, 'r') as f:
        content = f.read()

    new_content = add_imports(content, os.path.basename(filepath))

    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python add_base_class_imports.py <directory>")
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

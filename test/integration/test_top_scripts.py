"""Test that the top-level scripts run"""

import os
import util.misc
import util.file

def test_top_scripts_help():
    """Test that calling each top-level script works, and prints a help message"""

    for top_script in ('assembly.py', 'broad_utils.py', 'illumina.py', 'interhost.py', 'intrahost.py',
                       'metagenomics.py', 'ncbi.py', 'read_utils.py', 'reports.py', 'taxon_filter.py'):
        result = util.misc.run_and_print([os.path.join(util.file.get_project_path(), top_script), '-h'], check=True)
        assert result.returncode == 0
        assert result.stdout

    

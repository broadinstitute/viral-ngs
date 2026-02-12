'''geNomad viral and plasmid identification tool'''
import logging
import os
import os.path
import shutil
import subprocess

from viral_ngs import core
from viral_ngs.core import file
from viral_ngs.core import misc
from builtins import super

log = logging.getLogger(__name__)

class Genomad(core.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(core.PrexistingUnixCommand(shutil.which('genomad'), require_executability=False))
        super(Genomad, self).__init__(install_methods=install_methods)

    def version(self):
        '''Parse version from genomad --version output'''
        try:
            result = subprocess.run(['genomad', '--version'], capture_output=True, text=True, check=True)
            version_line = result.stdout.strip() or result.stderr.strip()
            # geNomad outputs "geNomad, version X.Y.Z"
            if version_line:
                parts = version_line.split()
                if len(parts) >= 1:
                    # Return the last part, stripping any trailing comma
                    return parts[-1].rstrip(',')
            return 'unknown'
        except (subprocess.CalledProcessError, FileNotFoundError, IndexError):
            return 'unknown'

    def execute(self, command, args=None, options=None):
        '''Run a genomad command.

        Args:
          command: Command to run (e.g., 'genomad').
          args: List of positional args.
          options: Dict of keyword options.
        '''
        options = options or {}
        args = args or []

        cmd = [command]
        cmd.extend(args)

        for key, value in options.items():
            if value is None:
                cmd.append(key)
            else:
                cmd.extend([key, str(value)])

        log.debug('Calling %s: %s', command, ' '.join(cmd))

        subprocess.check_call(cmd)

    def _is_fasta_empty(self, fasta_file):
        '''Check if FASTA file is empty or does not exist.

        Args:
          fasta_file: Path to FASTA file.

        Returns:
          True if file is None, does not exist, or is empty; False otherwise.
        '''
        if fasta_file is None:
            return True
        if not os.path.exists(fasta_file):
            return True
        if os.path.getsize(fasta_file) == 0:
            return True

        # For small files, check content
        file_size = os.path.getsize(fasta_file)
        if file_size <= 1024:
            with open(fasta_file, 'r') as f:
                content = f.read()
                if content.strip() == '':
                    return True

        # For files >1024 bytes, assume non-empty for performance
        return False

    def end_to_end(self, in_fasta, db_path, out_dir, num_threads=None):
        '''Run geNomad end-to-end pipeline on input FASTA.

        Args:
          in_fasta: Input FASTA file with sequences to classify.
          db_path: Path to geNomad database directory.
          out_dir: Output directory for geNomad results.
          num_threads: Number of threads to use (optional).

        Raises:
          ValueError: If database path does not exist or is not a directory.
        '''
        # Validate database path first
        if not os.path.isdir(db_path):
            raise ValueError(f"Database path does not exist or is not a directory: {db_path}")

        # Create output directory
        file.mkdir_p(out_dir)

        # Check for empty input FASTA
        if self._is_fasta_empty(in_fasta):
            log.warning("Input FASTA is empty, skipping geNomad")
            return

        # Build options
        opts = {}
        if num_threads is not None:
            opts['--threads'] = misc.sanitize_thread_count(num_threads)

        # Execute geNomad with correct argument order: end-to-end INPUT OUTPUT DATABASE
        self.execute('genomad', ['end-to-end', in_fasta, out_dir, db_path], options=opts)

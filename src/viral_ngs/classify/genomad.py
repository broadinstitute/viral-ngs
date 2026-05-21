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
        '''Check if FASTA file is empty.

        Args:
          fasta_file: Path to FASTA file.

        Returns:
          True if file is empty; False otherwise.
        '''
        file_size = os.path.getsize(fasta_file)
        if file_size == 0:
            return True

        # For small files, check content
        if file_size <= 1024:
            with file.open_or_gzopen(fasta_file, 'rt') as f:
                content = f.read()
                if content.strip() == '':
                    return True

        # For files >1024 bytes, assume non-empty for performance
        return False

    def end_to_end(self, in_fasta, db_path, out_dir, num_threads=None,
                   cleanup=False, restart=False, filter_preset=None,
                   enable_score_calibration=False, composition=None,
                   min_score=None, max_fdr=None, min_number_genes=None,
                   max_uscg=None, splits=None):
        '''Run geNomad end-to-end pipeline on input FASTA.

        Args:
          in_fasta: Input FASTA file with sequences to classify.
          db_path: Path to geNomad database directory.
          out_dir: Output directory for geNomad results.
          num_threads: Number of threads to use (optional).
          cleanup: Delete intermediate files after execution.
          restart: Overwrite existing intermediate files.
          filter_preset: Summary filtering preset; one of conservative or relaxed.
          enable_score_calibration: Execute score calibration module.
          composition: Sample composition for score calibration.
          min_score: Minimum score to flag a sequence as virus or plasmid.
          max_fdr: Maximum accepted false discovery rate.
          min_number_genes: Minimum number of genes required for classification.
          max_uscg: Maximum allowed universal single copy genes.
          splits: Split data for MMseqs2 marker search to reduce memory usage.

        Raises:
          FileNotFoundError: If input FASTA does not exist.
          ValueError: If database path does not exist or is not a directory.
        '''
        # Validate inputs first
        if not in_fasta or not os.path.isfile(in_fasta):
            raise FileNotFoundError(in_fasta)
        if not os.path.isdir(db_path):
            raise ValueError(f"Database path does not exist or is not a directory: {db_path}")
        if filter_preset is not None:
            filter_preset = filter_preset.lower()
        if composition is not None:
            composition = composition.lower()
        self._validate_end_to_end_options(
            filter_preset=filter_preset,
            composition=composition,
            min_score=min_score,
            max_fdr=max_fdr,
            min_number_genes=min_number_genes,
            max_uscg=max_uscg,
            splits=splits,
        )

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
        if cleanup:
            opts['--cleanup'] = None
        if restart:
            opts['--restart'] = None
        if filter_preset == 'conservative':
            opts['--conservative'] = None
        elif filter_preset == 'relaxed':
            opts['--relaxed'] = None
        if enable_score_calibration:
            opts['--enable-score-calibration'] = None
        if composition is not None:
            opts['--composition'] = composition
        if min_score is not None:
            opts['--min-score'] = min_score
        if max_fdr is not None:
            opts['--max-fdr'] = max_fdr
        if min_number_genes is not None:
            opts['--min-number-genes'] = min_number_genes
        if max_uscg is not None:
            opts['--max-uscg'] = max_uscg
        if splits is not None:
            opts['--splits'] = splits

        # Execute geNomad with correct argument order: end-to-end INPUT OUTPUT DATABASE
        self.execute('genomad', ['end-to-end', in_fasta, out_dir, db_path], options=opts)

    def _validate_end_to_end_options(self, filter_preset=None, composition=None,
                                     min_score=None, max_fdr=None,
                                     min_number_genes=None, max_uscg=None,
                                     splits=None):
        valid_presets = (None, 'conservative', 'relaxed')
        if filter_preset not in valid_presets:
            raise ValueError("filter_preset must be one of: conservative, relaxed")

        explicit_filters = {
            'min_score': min_score,
            'max_fdr': max_fdr,
            'min_number_genes': min_number_genes,
            'max_uscg': max_uscg,
        }
        if filter_preset and any(value is not None for value in explicit_filters.values()):
            raise ValueError(
                "filter_preset cannot be combined with explicit filtering options"
            )

        valid_compositions = (None, 'auto', 'metagenome', 'virome')
        if composition not in valid_compositions:
            raise ValueError("composition must be one of: auto, metagenome, virome")
        if min_score is not None and not 0 <= min_score <= 1:
            raise ValueError("min_score must be between 0 and 1")
        if max_fdr is not None and not 0 <= max_fdr <= 1:
            raise ValueError("max_fdr must be between 0 and 1")
        if min_number_genes is not None and min_number_genes < 0:
            raise ValueError("min_number_genes must be non-negative")
        if max_uscg is not None and max_uscg < 0:
            raise ValueError("max_uscg must be non-negative")
        if splits is not None and splits < 0:
            raise ValueError("splits must be non-negative")

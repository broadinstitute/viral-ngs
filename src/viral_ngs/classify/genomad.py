'''geNomad viral and plasmid identification tool'''
import logging
import os
import os.path
import shutil
import subprocess
import tempfile

from viral_ngs import core
import viral_ngs.core.duckdb_utils as duckdb_utils
from viral_ngs.core import file
from viral_ngs.core import misc
from builtins import super
from . import contig_alignments

log = logging.getLogger(__name__)

VIRUS_SUMMARY_COLUMNS = [
    "seq_name", "length", "topology", "coordinates", "n_genes", "genetic_code",
    "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy",
]
READ_LABEL_COLUMNS = [
    "read_id", "read_length", "contig_id", "seq_name", "mapping_quality",
    "pct_identity", "pct_query_cov", "mapped_well", *VIRUS_SUMMARY_COLUMNS[1:],
]


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


def _load_virus_summary(con, virus_summary):
    con.execute(
        "CREATE TABLE virus_raw AS SELECT * FROM read_csv(?, delim='\t', header=true, all_varchar=true)",
        [virus_summary],
    )
    columns = {row[1] for row in con.execute("PRAGMA table_info('virus_raw')").fetchall()}
    missing = [name for name in VIRUS_SUMMARY_COLUMNS if name not in columns]
    if missing:
        raise ValueError("geNomad virus summary missing required columns: " + ", ".join(missing))
    duplicate = con.execute(
        "SELECT seq_name FROM virus_raw GROUP BY seq_name HAVING count(*) > 1 LIMIT 1"
    ).fetchone()
    if duplicate is not None:
        raise ValueError("geNomad virus summary has duplicate seq_name: {}".format(duplicate[0]))

    con.execute(
        r"""
        CREATE TABLE virus AS
        SELECT {},
            regexp_replace(seq_name, '\|provirus_[0-9]+_[0-9]+$', '') AS contig_id,
            CASE WHEN regexp_matches(coordinates, '^[0-9]+-[0-9]+$')
                 THEN CAST(split_part(coordinates, '-', 1) AS BIGINT) END AS region_start,
            CASE WHEN regexp_matches(coordinates, '^[0-9]+-[0-9]+$')
                 THEN CAST(split_part(coordinates, '-', 2) AS BIGINT) END AS region_end
        FROM virus_raw
        """.format(", ".join(VIRUS_SUMMARY_COLUMNS))
    )
    con.execute("DROP TABLE virus_raw")


def _write_empty_read_labels(output_tsv):
    with file.open_or_gzopen(output_tsv, "wt") as outf:
        outf.write("\t".join(READ_LABEL_COLUMNS) + "\n")
    log.info("%s reads labeled with a geNomad virus call", 0)


def label_reads_by_contig(
    aligned_bam, virus_summary, output_tsv, min_mapq=5, min_identity=90.0,
    min_query_cov=80.0, duckdb_memory_limit=None, work_dir=None,
):
    """Label each read with its best primary alignment overlapping a geNomad virus call."""
    for path in (aligned_bam, virus_summary):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    min_mapq = int(min_mapq)
    min_identity = float(min_identity)
    min_query_cov = float(min_query_cov)
    contig_alignments.validate_percent_threshold("min_identity", min_identity)
    contig_alignments.validate_percent_threshold("min_query_cov", min_query_cov)
    duckdb = duckdb_utils.import_duckdb("geNomad read-by-contig labeling")

    file.mkdir_p(os.path.dirname(os.path.abspath(output_tsv)))
    if work_dir:
        file.mkdir_p(work_dir)
    with tempfile.TemporaryDirectory(prefix="genomad_reads_", dir=work_dir) as tmp_dir:
        if os.path.getsize(virus_summary) == 0:
            log.info("geNomad virus summary %s is empty; no reads labeled", virus_summary)
            _write_empty_read_labels(output_tsv)
            return
        normalized, n_written, _ = contig_alignments.normalize_primary_alignments(aligned_bam, tmp_dir)
        if not n_written:
            log.info("%s primary alignments retained", 0)
            _write_empty_read_labels(output_tsv)
            return

        spill_dir = os.path.join(tmp_dir, "duckdb_spill")
        file.mkdir_p(spill_dir)
        result_tsv = os.path.join(tmp_dir, "reads_labeled.tsv")
        con = duckdb_utils.connect(spill_dir, memory_limit=duckdb_memory_limit, duckdb_module=duckdb)
        try:
            n_primary = contig_alignments.load_primary_alignments(
                con, normalized, min_mapq, min_identity, min_query_cov
            )
            log.info("%s primary alignments retained", n_primary)
            _load_virus_summary(con, virus_summary)
            # BAM spans are 0-based half-open; geNomad intervals are 1-based inclusive.
            con.execute(
                """
                CREATE TABLE result AS
                SELECT a.query_name AS read_id, a.query_length AS read_length, a.target_name AS contig_id,
                       a.mapping_quality, a.pct_identity, a.pct_query_cov, a.mapped_well,
                       v.* EXCLUDE (contig_id, region_start, region_end)
                FROM alignment_filtered a JOIN virus v ON a.target_name = v.contig_id
                    AND (v.region_start IS NULL OR (a.target_start < v.region_end AND a.target_end >= v.region_start))
                QUALIFY ROW_NUMBER() OVER (
                    PARTITION BY a.query_name ORDER BY a.mapping_quality DESC, a.pct_identity DESC,
                        TRY_CAST(v.virus_score AS DOUBLE) DESC NULLS LAST, a.source_order ASC, v.rowid ASC
                ) = 1
                """
            )
            n_reads = con.execute("SELECT count(*) FROM result").fetchone()[0]
            if not n_reads:
                _write_empty_read_labels(output_tsv)
                return
            con.execute(
                "COPY (SELECT {} FROM result) TO ? (FORMAT CSV, DELIMITER '\t', HEADER TRUE)".format(
                    ", ".join(READ_LABEL_COLUMNS)
                ),
                [result_tsv],
            )
        finally:
            con.close()

        with open(result_tsv, "rt") as inf, file.open_or_gzopen(output_tsv, "wt") as outf:
            shutil.copyfileobj(inf, outf)
        log.info("%s reads labeled with a geNomad virus call", n_reads)

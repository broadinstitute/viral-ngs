'''
Kallisto/kb-python classification tool.

This wrapper exposes viral-ngs kallisto commands while using kb-python's `kb`
CLI to orchestrate kallisto and bustools. Count outputs include a long-form
collapsed hit table (`counts.tsv`) and the kb-generated h5ad needed by
`kb extract`; extract outputs include the downstream read summary (`summary.tsv`).
'''
import csv
import gzip
import glob
import itertools
import json
import logging
import os
import os.path
import shutil
import subprocess

import anndata
import numpy as np

from viral_ngs import core
from viral_ngs.core import picard
from viral_ngs.core import samtools
from viral_ngs.core import file
from viral_ngs.core import misc
from builtins import super

log = logging.getLogger(__name__)
_FASTQ_SUFFIXES = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

class Kallisto(core.Tool):
    """Run kallisto workflows through kb-python and normalize public outputs."""

    SUBCOMMANDS = ['count', 'ref', 'extract']
    COUNTS_TSV = 'counts.tsv'
    H5AD = 'adata.h5ad'
    READ_HITS_TSV = 'read_hits.tsv'
    SUMMARY_TSV = 'summary.tsv'
    COUNTS_HEADER = ('sample_id', 'db_hit_id', 'count')
    READ_HITS_HEADER = ('read_id', 'db_hit_id')
    SUMMARY_HEADER = ('SAMPLE_ID', 'READ_ID', 'DB_ID', 'TAXONOMY_LINEAGE', 'TAXONOMY_NAME', 'SEQUENCE_LENGTH')

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(core.PrexistingUnixCommand(shutil.which('kb'), require_executability=False))
        super(Kallisto, self).__init__(install_methods=install_methods)

    def version(self):
        try:
            result = subprocess.run(
                ['micromamba', 'list', 'kb-python', '--json'],
                capture_output=True,
                text=True,
                check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return 'unknown'

        try:
            packages = json.loads(result.stdout or '[]')
        except json.JSONDecodeError:
            return 'unknown'

        for package in packages:
            if package.get('name') == 'kb-python' and package.get('version'):
                return package['version']

        return 'unknown'

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())


    def execute(self, command,output, args=None, options=None):
        '''Run a kallisto command via the kb executable.

        Args:
          command: Subcommand to run.
          output: Output file to send to command.
          args: List of positional args.
          options: List of keyword options. Values can be single items or lists for multi-value options.
        '''
        options = options or {}

        if output:
            options['-o'] = output
        args = args or []

        cmd = command.split()

        # Build options, handling both single values and lists
        for key, value in options.items():
            if value is None:
                # Empty flag like --aa
                cmd.append(key)
            elif isinstance(value, list):
                # Multi-value option like -ts target1 target2
                cmd.append(key)
                cmd.extend([str(v) for v in value])
            else:
                # Single value option
                cmd.extend([key, str(value)])
        
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        # Use Popen to capture both stdout and stderr for better error reporting
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        # Log output
        if stdout:
            log.info("kb output: %s", stdout)
        if stderr:
            log.info("kb stderr: %s", stderr)
            
        # Check for errors: the kb backend sometimes catches exceptions without proper exit codes
        # Look for Python tracebacks or subprocess errors in stderr
        has_error = stderr and ('Traceback (most recent call last)' in stderr or 'CalledProcessError' in stderr)
        
        # Raise exception if command failed or error detected
        if process.returncode != 0 or has_error:
            error_msg = f"Command failed with return code {process.returncode}: {' '.join(cmd)}"
            if stderr:
                error_msg += f"\nStderr: {stderr}"
            log.error(error_msg)
            raise subprocess.CalledProcessError(process.returncode if process.returncode != 0 else 1, cmd, output=stdout, stderr=stderr)
        
    def build(self, ref_fasta, index, workflow='standard', kmer_len=31,  protein=False, num_threads=None):
        '''Create a kallisto index.
        Args:
          ref_fasta: reference fasta file
          index: output index file
          workflow: one of 'standard', 'nac', 'kite', 'custom'
          kmer_len: kmer size (default 31)
          protein: ref_fasta file contains amino acid sequences
        '''
        # build db
        build_opts = {
            '-t': misc.sanitize_thread_count(num_threads)
        }
        if kmer_len:
            build_opts['-k'] = kmer_len
        if index:
            build_opts['-i'] = index
        if protein:
            build_opts['--aa'] = None
        if workflow:
            build_opts['--workflow'] = workflow
        self.execute('kb ref', None, args=[ref_fasta], options=build_opts)
        
    def classify(self, in_bam, index_file, out_dir, t2g_file, k=31, parity='single', technology='bulk', loom=False, protein=False, num_threads=None, sample_name=None):
        """Classify input reads (bam)

        Args:
          in_bam: unaligned reads (single or paired-end) in BAM or FASTQ format.
          index_file: kallisto index file
          out_dir: output directory
          t2g_file: transcript to gene mapping file
          k: kmer size (default 31)
          parity: one of 'single', 'paired'
          technology: one of '10xv2', '10xv3', '10xv3cr', '10xv3multiome', '10xv1', 'dropseq', 'indrop', 'sci-rna-seq', 'bulk'
          loom: output loom file
          protein: ref_fasta file contains amino acid sequences
          num_threads: number of threads to use
          sample_name: sample name to use in counts.tsv
        """
        opts = {
            '-i': index_file,
            '-g': t2g_file,
            '-t': misc.sanitize_thread_count(num_threads),
            '--parity': parity
        }
        if k:
            opts['-k'] = k
        if technology:
            opts['-x'] = technology
        opts['--h5ad'] = None
        if loom:
            opts['--loom'] = None
        if protein:
            opts['--aa'] = None

        if not os.path.exists(in_bam):
            raise FileNotFoundError(in_bam)

        sample_name = sample_name or self._sample_name_from_input(in_bam)
        if self._is_fastq(in_bam):
            self.execute('kb count', out_dir, args=[in_bam], options=opts)
        else:
            if samtools.SamtoolsTool().isEmpty(in_bam):
                self._write_empty_count_outputs(out_dir, sample_name)
                return

            tmp_fastq1 = file.mkstempfname('.1.fastq')
            tmp_fastq2 = file.mkstempfname('.2.fastq')
            tmp_fastq3 = file.mkstempfname('.s.fastq')
            tmp_interleaved = None
            try:
                # Do not convert this to samtools bam2fq unless we can figure out how to replicate
                # the clipping functionality of Picard SamToFastq
                picard_tool = picard.SamToFastqTool()
                picard_opts = {
                    'CLIPPING_ATTRIBUTE': picard.SamToFastqTool.illumina_clipping_attribute,
                    'CLIPPING_ACTION': 'X'
                }
                picard_tool.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                            picardOptions=picard.PicardTools.dict_to_picard_opts(picard_opts),
                            JVMmemory=picard_tool.jvmMemDefault)

                # Detect if input bam was paired by checking fastq 2
                if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
                    self.execute('kb count', out_dir, args=[tmp_fastq3], options=opts)
                else:
                    tmp_interleaved = file.mkstempfname('.interleaved.fastq')
                    with open(tmp_fastq1, 'rb') as fastq1, open(tmp_fastq2, 'rb') as fastq2, open(tmp_interleaved, 'wb') as interleaved:
                        while True:
                            read1 = [fastq1.readline() for _ in range(4)]
                            if not read1[0]:
                                break
                            if any(line == b'' for line in read1[1:]):
                                raise ValueError("Unexpected end of read 1 FASTQ while interleaving paired data")
                            read2 = [fastq2.readline() for _ in range(4)]
                            if any(line == b'' for line in read2):
                                raise ValueError("Unexpected end of read 2 FASTQ while interleaving paired data")
                            interleaved.writelines(read1)
                            interleaved.writelines(read2)
                        if fastq2.readline():
                            raise ValueError("Read 2 FASTQ contains extra data after interleaving paired data")

                    self.execute('kb count', out_dir, args=[tmp_interleaved], options=opts)
            finally:
                for path in (tmp_fastq1, tmp_fastq2, tmp_fastq3, tmp_interleaved):
                    if path and os.path.exists(path):
                        try:
                            os.unlink(path)
                        except OSError as e:
                            log.warning("Failed to delete temporary file %s: %s", path, e)

        self._finalize_count_outputs(out_dir, sample_name)
        
    def extract(self, in_bam, index_file, target_ids, out_dir, t2g_file, protein=False, num_threads=None, sample_name=None, id_to_tax_map=None, taxonomy_level='highest'):
        """Extracts reads mapping to target ids from input reads (bam)
        
        Args:
          in_bam: unaligned read to extract reads from (FASTQ or BAM)
          index_file: kallisto index file
          out_dir: output directory
          t2g_file: transcript to gene mapping file
          protein: ref_fasta file contains amino acid sequences
          target_ids: list of target ids to extract
          num_threads: number of threads to use
          sample_name: sample name to use in summary.tsv
          id_to_tax_map: optional mapping from target IDs to taxonomy lineage
          taxonomy_level: taxonomy name to report from the mapping, one of highest or deepest
        """
        target_ids = [target_id for target_id in (target_ids or []) if target_id]
        opts = {
            '-i': index_file,
            '-g': t2g_file,
            '-ts': target_ids,  # Pass as list for multi-value option
            '-t': misc.sanitize_thread_count(num_threads)
        }
        if protein:
            opts['--aa'] = None

        if not os.path.exists(in_bam):
            raise FileNotFoundError(in_bam)

        sample_name = sample_name or self._sample_name_from_input(in_bam)
        if not target_ids:
            self._write_empty_extract_outputs(out_dir)
            return

        if self._is_fastq(in_bam):
            self.execute('kb extract', out_dir, args=[in_bam], options=opts)
        else:
            if samtools.SamtoolsTool().isEmpty(in_bam):
                self._write_empty_extract_outputs(out_dir)
                return

            tmp_fastq1 = file.mkstempfname('.1.fastq')
            tmp_fastq2 = file.mkstempfname('.2.fastq')
            tmp_fastq3 = file.mkstempfname('.s.fastq')
            tmp_interleaved = None
            try:
                # Do not convert this to samtools bam2fq unless we can figure out how to replicate
                # the clipping functionality of Picard SamToFastq
                picard_tool = picard.SamToFastqTool()
                picard_opts = {
                    'CLIPPING_ATTRIBUTE': picard.SamToFastqTool.illumina_clipping_attribute,
                    'CLIPPING_ACTION': 'X'
                }
                picard_tool.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                            picardOptions=picard.PicardTools.dict_to_picard_opts(picard_opts),
                            JVMmemory=picard_tool.jvmMemDefault)

                # Detect if input bam was paired by checking fastq 2
                if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
                    self.execute('kb extract', out_dir, args=[tmp_fastq3], options=opts)
                else:
                    tmp_interleaved = file.mkstempfname('.interleaved.fastq')
                    with open(tmp_fastq1, 'rb') as fastq1, open(tmp_fastq2, 'rb') as fastq2, open(tmp_interleaved, 'wb') as interleaved:
                        while True:
                            read1 = [fastq1.readline() for _ in range(4)]
                            if not read1[0]:
                                break
                            if any(line == b'' for line in read1[1:]):
                                raise ValueError("Unexpected end of read 1 FASTQ while interleaving paired data")
                            read2 = [fastq2.readline() for _ in range(4)]
                            if any(line == b'' for line in read2):
                                raise ValueError("Unexpected end of read 2 FASTQ while interleaving paired data")
                            interleaved.writelines(read1)
                            interleaved.writelines(read2)
                        if fastq2.readline():
                            raise ValueError("Read 2 FASTQ contains extra data after interleaving paired data")

                    self.execute('kb extract', out_dir, args=[tmp_interleaved], options=opts)
            finally:
                for path in (tmp_fastq1, tmp_fastq2, tmp_fastq3, tmp_interleaved):
                    if path and os.path.exists(path):
                        try:
                            os.unlink(path)
                        except OSError as e:
                            log.warning("Failed to delete temporary file %s: %s", path, e)

        self._write_extract_tsvs(out_dir, target_ids, sample_name, id_to_tax_map=id_to_tax_map, taxonomy_level=taxonomy_level)

    @staticmethod
    def _is_fastq(input_path):
        return input_path.lower().endswith(_FASTQ_SUFFIXES)

    @staticmethod
    def _sample_name_from_input(input_path):
        sample_name = os.path.basename(input_path)
        for suffix in _FASTQ_SUFFIXES + ('.bam',):
            if sample_name.lower().endswith(suffix):
                return sample_name[:-len(suffix)]
        return os.path.splitext(sample_name)[0]

    def _finalize_count_outputs(self, out_dir, sample_name):
        h5ad_file = self._find_count_h5ad(out_dir)
        self._add_sample_metadata_to_h5ad(h5ad_file, sample_name=sample_name)

        exposed_h5ad = os.path.join(out_dir, self.H5AD)
        if os.path.abspath(h5ad_file) != os.path.abspath(exposed_h5ad):
            shutil.copyfile(h5ad_file, exposed_h5ad)

        self.write_counts_tsv_from_h5ad(exposed_h5ad, os.path.join(out_dir, self.COUNTS_TSV))

    def _find_count_h5ad(self, out_dir):
        h5ad_files = glob.glob(os.path.join(out_dir, "counts_unfiltered", "*.h5ad"))
        if len(h5ad_files) == 0:
            raise FileNotFoundError(f"No .h5ad file found in counts_unfiltered/ directory of {out_dir}")
        if len(h5ad_files) > 1:
            raise ValueError(f"Expected exactly one .h5ad file in counts_unfiltered/ directory of {out_dir}, found {len(h5ad_files)}")
        return h5ad_files[0]

    def _write_empty_count_outputs(self, out_dir, sample_name):
        file.mkdir_p(out_dir)

        adata = anndata.AnnData(np.zeros((1, 0), dtype=np.int64))
        adata.obs_names = [sample_name]
        adata.obs['sample'] = [sample_name]
        adata.obs['batch_name'] = [sample_name]

        exposed_h5ad = os.path.join(out_dir, self.H5AD)
        adata.write_h5ad(exposed_h5ad)
        self.write_counts_tsv_from_h5ad(exposed_h5ad, os.path.join(out_dir, self.COUNTS_TSV))

    def _write_empty_extract_outputs(self, out_dir):
        file.mkdir_p(out_dir)
        with file.open_or_gzopen(os.path.join(out_dir, self.READ_HITS_TSV), 'wt') as read_hits_outf:
            csv.writer(read_hits_outf, delimiter='\t', lineterminator='\n').writerow(self.READ_HITS_HEADER)
        with file.open_or_gzopen(os.path.join(out_dir, self.SUMMARY_TSV), 'wt') as summary_outf:
            csv.writer(summary_outf, delimiter='\t', lineterminator='\n').writerow(self.SUMMARY_HEADER)

    def write_counts_tsv_from_h5ad(self, h5ad_file, out_tsv):
        """Write long-form collapsed kallisto counts from a kb-generated h5ad file.

        The output has columns `sample_id`, `db_hit_id`, and `count`, with one
        row per nonzero collapsed hit.
        """
        adata = anndata.read_h5ad(h5ad_file)
        sample_ids = self._sample_ids_from_h5ad(adata)
        hit_ids = adata.var.index.tolist()

        with file.open_or_gzopen(out_tsv, 'wt') as outf:
            writer = csv.writer(outf, delimiter='\t', lineterminator='\n')
            writer.writerow(self.COUNTS_HEADER)
            if hasattr(adata.X, 'getrow'):
                for row_idx in range(adata.n_obs):
                    row = adata.X.getrow(row_idx).tocoo()
                    for col_idx, count in sorted(zip(row.col, row.data)):
                        if count != 0:
                            writer.writerow((sample_ids[row_idx], hit_ids[col_idx], self._format_count(count)))
            else:
                for row_idx, row in enumerate(adata.X):
                    for col_idx, count in enumerate(row):
                        if count != 0:
                            writer.writerow((sample_ids[row_idx], hit_ids[col_idx], self._format_count(count)))

    @staticmethod
    def _sample_ids_from_h5ad(adata):
        if 'sample' in adata.obs:
            return [str(sample) for sample in adata.obs['sample'].tolist()]
        return [str(obs_name) for obs_name in adata.obs_names.tolist()]

    @staticmethod
    def _format_count(count):
        as_float = float(count)
        if not as_float.is_integer():
            raise ValueError(f"Non-integer kallisto count: {count!r}")
        return str(int(as_float))

    def write_top_taxa_report_from_counts_tsv(self, counts_tsv, out_report, id_to_tax_map=None, target_taxon='Viruses'):
        """Write a ranked focal-taxon report from long-form kallisto counts.

        The output schema intentionally matches the historical kb/kallisto
        top-taxa report consumed by WDL scalar outputs.
        """
        counts_by_hit = self._read_counts_tsv(counts_tsv)
        taxonomy_map = self._load_id_to_tax_map(id_to_tax_map, taxonomy_level='deepest')
        rows = []

        for db_hit_id, hit_reads in sorted(counts_by_hit.items()):
            if hit_reads <= 0:
                continue
            if id_to_tax_map:
                taxonomy_lineage, taxonomy_name = self._taxonomy_for_hit(db_hit_id, taxonomy_map, taxonomy_map_provided=True)
                if not self._taxonomy_lineage_contains(taxonomy_lineage, target_taxon):
                    continue
                hit_name = taxonomy_name
            else:
                taxonomy_name = db_hit_id
                hit_name = db_hit_id

            rows.append({
                'palmdb_id': db_hit_id,
                'hit_id': hit_name,
                'hit_lowest_taxa_name': taxonomy_name,
                'hit_reads': int(hit_reads),
            })

        total_focal_count = sum(row['hit_reads'] for row in rows)
        if total_focal_count == 0:
            out_rows = [{
                'focal_taxon_name': target_taxon,
                'focal_taxon_count': 0,
                'palmdb_id': '',
                'hit_id': '',
                'hit_lowest_taxa_name': '',
                'hit_reads': 0,
                'pct_of_focal': 0.0,
            }]
        else:
            out_rows = []
            for row in sorted(rows, key=lambda item: (-item['hit_reads'], item['palmdb_id'])):
                out_row = {
                    'focal_taxon_name': target_taxon,
                    'focal_taxon_count': total_focal_count,
                    'pct_of_focal': 100.0 * row['hit_reads'] / total_focal_count,
                }
                out_row.update(row)
                out_rows.append(out_row)

        header = ('focal_taxon_name', 'focal_taxon_count', 'palmdb_id', 'hit_id', 'hit_lowest_taxa_name', 'hit_reads', 'pct_of_focal')
        with file.open_or_gzopen(out_report, 'wt') as outf:
            writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
            writer.writeheader()
            writer.writerows(out_rows)

    @classmethod
    def _read_counts_tsv(cls, counts_tsv):
        counts_by_hit = {}
        with file.open_or_gzopen(counts_tsv, 'rt') as inf:
            reader = csv.DictReader(inf, delimiter='\t')
            required_columns = {'sample_id', 'db_hit_id', 'count'}
            if reader.fieldnames is None:
                raise ValueError(f"Empty kallisto counts TSV: {counts_tsv}")
            missing_columns = required_columns.difference(reader.fieldnames)
            if missing_columns:
                raise ValueError(f"Kallisto counts TSV missing required columns: {', '.join(sorted(missing_columns))}")
            for row in reader:
                db_hit_id = row['db_hit_id'].strip()
                if not db_hit_id:
                    raise ValueError("Kallisto counts TSV contains an empty db_hit_id")
                count = cls._parse_counts_tsv_count(row['count'], db_hit_id)
                counts_by_hit[db_hit_id] = counts_by_hit.get(db_hit_id, 0) + count
        return counts_by_hit

    @staticmethod
    def _parse_counts_tsv_count(raw_count, db_hit_id):
        try:
            count = float(raw_count)
        except ValueError:
            raise ValueError(f"Invalid count for {db_hit_id!r}: {raw_count!r}")
        if count < 0:
            raise ValueError(f"Negative count for {db_hit_id!r}: {raw_count!r}")
        if not count.is_integer():
            raise ValueError(f"Non-integer count for {db_hit_id!r}: {raw_count!r}")
        return int(count)

    @staticmethod
    def _taxonomy_lineage_contains(taxonomy_lineage, target_taxon):
        target = target_taxon.strip().lower()
        return any(value.strip().lower() == target for value in taxonomy_lineage.split(';'))

    def _write_extract_tsvs(self, out_dir, target_ids, sample_name, id_to_tax_map=None, taxonomy_level='highest'):
        if not target_ids:
            raise ValueError("target_ids must be provided when writing read hit TSV")
        taxonomy_map = self._load_id_to_tax_map(id_to_tax_map, taxonomy_level)
        target_id_set = set(target_ids)
        out_tsv = os.path.join(out_dir, self.READ_HITS_TSV)
        summary_tsv = os.path.join(out_dir, self.SUMMARY_TSV)
        with file.open_or_gzopen(out_tsv, 'wt') as read_hits_outf, file.open_or_gzopen(summary_tsv, 'wt') as summary_outf:
            read_hits_writer = csv.writer(read_hits_outf, delimiter='\t', lineterminator='\n')
            summary_writer = csv.writer(summary_outf, delimiter='\t', lineterminator='\n')
            read_hits_writer.writerow(self.READ_HITS_HEADER)
            summary_writer.writerow(self.SUMMARY_HEADER)
            for fastq_path, db_hit_id in self._iter_extracted_fastqs(out_dir, target_id_set):
                taxonomy_lineage, taxonomy_name = self._taxonomy_for_hit(db_hit_id, taxonomy_map, id_to_tax_map is not None)
                for read_id, sequence_length in self._iter_fastq_records(fastq_path):
                    read_hits_writer.writerow((read_id, db_hit_id))
                    summary_writer.writerow((sample_name, read_id, db_hit_id, taxonomy_lineage, taxonomy_name, sequence_length))

    @classmethod
    def _load_id_to_tax_map(cls, id_to_tax_map, taxonomy_level='highest'):
        if taxonomy_level not in ('highest', 'deepest'):
            raise ValueError("taxonomy_level must be one of: highest, deepest")
        if not id_to_tax_map:
            return {}

        taxonomy_map = {}
        with file.open_or_gzopen(id_to_tax_map, 'rt') as inf:
            first_line = inf.readline()
            if not first_line:
                raise ValueError(f"Empty ID to taxonomy mapping file: {id_to_tax_map}")
            delimiter = cls._detect_mapping_delimiter(id_to_tax_map, first_line)
            first_row = next(csv.reader([first_line], delimiter=delimiter))
            has_header = cls._looks_like_taxonomy_header(first_row)
            if has_header:
                tax_start, tax_end = cls._taxonomy_column_bounds(first_row)
                row_iter = csv.reader(inf, delimiter=delimiter)
            else:
                tax_start, tax_end = cls._taxonomy_row_bounds(first_row)
                row_iter = itertools.chain([first_row], csv.reader(inf, delimiter=delimiter))

            for row in row_iter:
                if not row:
                    continue
                if len(row) <= tax_start:
                    raise ValueError(f"Malformed taxonomy mapping row for {row[0]!r}: expected taxonomy columns")
                db_hit_id = row[0]
                tax_values = cls._taxonomy_values_from_row(row, tax_start, tax_end)
                taxonomy_map[db_hit_id] = cls._format_taxonomy_values(tax_values, taxonomy_level)
        return taxonomy_map

    @staticmethod
    def _detect_mapping_delimiter(mapping_path, header_line):
        lower_path = mapping_path.lower()
        if lower_path.endswith(('.tsv', '.tsv.gz')):
            return '\t'
        if lower_path.endswith(('.csv', '.csv.gz')):
            return ','
        return '\t' if '\t' in header_line else ','

    @staticmethod
    def _looks_like_taxonomy_header(row):
        normalized = [column.strip().lower() for column in row]
        known_header_values = {
            'id',
            'palmdb_id',
            'db_hit_id',
            'gene_id',
            'transcript_id',
            'tax_level_1',
            'taxonomy',
            'strand',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species',
        }
        return any(column in known_header_values or column.startswith('tax_level_') for column in normalized)

    @staticmethod
    def _taxonomy_column_bounds(header):
        if len(header) < 2:
            raise ValueError("ID to taxonomy mapping must contain at least an ID column and one taxonomy column")
        normalized = [column.strip().lower() for column in header]
        tax_start = 2 if len(normalized) > 2 and normalized[0] == normalized[1] else 1
        tax_end = len(header) - 1 if Kallisto._is_strand_column_name(normalized[-1]) else len(header)
        if tax_start >= tax_end:
            raise ValueError("ID to taxonomy mapping does not contain taxonomy columns")
        return tax_start, tax_end

    @classmethod
    def _taxonomy_row_bounds(cls, row):
        if len(row) < 2:
            raise ValueError("ID to taxonomy mapping must contain at least an ID column and one taxonomy column")
        tax_start = 2 if len(row) > 2 and row[0].strip() == row[1].strip() else 1
        tax_end = len(row) - 1 if cls._is_strand_value(row[-1]) else len(row)
        if tax_start >= tax_end:
            raise ValueError("ID to taxonomy mapping does not contain taxonomy columns")
        return tax_start, tax_end

    @classmethod
    def _taxonomy_values_from_row(cls, row, tax_start, tax_end):
        tax_values = [value.strip() for value in row[tax_start:min(tax_end, len(row))]]
        if tax_values and cls._is_strand_value(tax_values[-1]):
            tax_values = tax_values[:-1]
        return [value for value in tax_values if value and value != '.']

    @staticmethod
    def _is_strand_column_name(column_name):
        return column_name in {'strand', 'sense', 'rna_type', 'genome_type'}

    @staticmethod
    def _is_strand_value(value):
        normalized = value.strip().lower()
        return normalized in {'+', '-', '+ssrna', '-ssrna', 'ssrna', 'dsrna', '+ssdna', '-ssdna', 'ssdna', 'dsdna'}

    @staticmethod
    def _format_taxonomy_values(tax_values, taxonomy_level):
        if not tax_values:
            return ('Unclassified RdRP', 'Unclassified RdRP')
        taxonomy_lineage = ';'.join(tax_values)
        taxonomy_name = tax_values[-1] if taxonomy_level == 'deepest' else tax_values[0]
        return taxonomy_lineage, taxonomy_name

    @staticmethod
    def _taxonomy_for_hit(db_hit_id, taxonomy_map, taxonomy_map_provided):
        if db_hit_id in taxonomy_map:
            return taxonomy_map[db_hit_id]
        if taxonomy_map_provided:
            return ('Unclassified RdRP', 'Unclassified RdRP')
        return ('', '')

    @staticmethod
    def _iter_extracted_fastqs(out_dir, target_ids):
        for db_hit_id in sorted(target_ids):
            target_dir = os.path.join(out_dir, db_hit_id)
            if not os.path.isdir(target_dir):
                continue
            for root, dirs, filenames in os.walk(target_dir):
                dirs[:] = sorted(dirs)
                for filename in sorted(filenames):
                    if filename.lower().endswith(_FASTQ_SUFFIXES):
                        yield os.path.join(root, filename), db_hit_id

    @staticmethod
    def _iter_fastq_records(fastq_path):
        open_fn = gzip.open if fastq_path.lower().endswith('.gz') else open
        with open_fn(fastq_path, 'rt') as inf:
            while True:
                header = inf.readline()
                if not header:
                    break
                sequence = inf.readline()
                plus = inf.readline()
                quality = inf.readline()
                if not sequence or not plus or not quality:
                    raise ValueError(f"Unexpected end of FASTQ record in {fastq_path}")
                if not header.startswith('@'):
                    raise ValueError(f"Expected FASTQ header line starting with @ in {fastq_path}: {header.rstrip()}")
                if not plus.startswith('+'):
                    raise ValueError(f"Expected FASTQ separator line starting with + in {fastq_path}: {plus.rstrip()}")
                yield header[1:].strip().split()[0], len(sequence.strip())

    def _add_sample_metadata_to_h5ad(self, h5ad_file, sample_name=None):
        """Add sample metadata to an h5ad file (in-place modification).
        
        Args:
          h5ad_file: path to h5ad file
          sample_name: optional sample name to use. If not provided, use the filename as fallback
        """
        if sample_name is None:
            sample_name = os.path.splitext(os.path.basename(h5ad_file))[0]
            log.warning(f"No sample name provided for {h5ad_file}, using filename as sample name: {sample_name}")

        adata = anndata.read_h5ad(h5ad_file)
        barcodes = None
        
        # Add sample metadata to observations
        adata.obs['sample'] = sample_name
        adata.obs['batch_name'] = sample_name
        
        # Add barcode info if available
        if barcodes is not None and len(barcodes) == adata.n_obs:
            adata.obs['batch_barcode'] = barcodes
        
        # Write back to the h5ad file
        adata.write_h5ad(h5ad_file)
        log.debug(f"Added sample metadata to {h5ad_file}: sample={sample_name}")

    def extract_hit_ids_from_h5ad(self, h5ad_file, threshold=1):
        """Parse h5ad file and extract all target IDs with 1 or more hits.
        
        Assumes h5ad contains a single sample (single row).

        Args:
          h5ad_file: path to h5ad file
          threshold: minimum count threshold to consider a target ID to return (default 1)

        Returns:
          List of target IDs (strings) with counts > 0
        """
        log.debug(f"Reading h5ad file: {h5ad_file}")
        adata = anndata.read_h5ad(h5ad_file)
        if adata.n_obs != 1:
            raise ValueError(
                f"Expected single-sample h5ad for kallisto extract target selection; found {adata.n_obs} observations"
            )

        gene_totals = self._h5ad_hit_totals(adata)
        gene_ids = adata.var.index.tolist()
        if threshold is None:
            return [gene_id for gene_id, count in zip(gene_ids, gene_totals) if count > 0]

        return [gene_id for gene_id, count in zip(gene_ids, gene_totals) if count >= threshold]

    @staticmethod
    def _h5ad_hit_totals(adata):
        gene_totals = adata.X.sum(axis=0)
        if hasattr(gene_totals, 'A1'):  # numpy matrix
            return gene_totals.A1
        if hasattr(gene_totals, 'toarray'):
            return gene_totals.toarray().ravel()
        return gene_totals

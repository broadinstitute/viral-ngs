'''
Centrifuger metagenomics classifier.
'''
import csv
import logging
import os
import shutil
import subprocess

import viral_ngs.core as core
import viral_ngs.core.file as util_file
import viral_ngs.core.misc as util_misc
import viral_ngs.core.picard as picard
import viral_ngs.core.samtools as samtools

log = logging.getLogger(__name__)

CLASSIFICATION_HEADER = [
    'readID',
    'seqID',
    'taxID',
    'score',
    '2ndBestScore',
    'hitLength',
    'queryLength',
    'numMatches',
]


class Centrifuger(core.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                core.PrexistingUnixCommand(
                    shutil.which('centrifuger'),
                    require_executability=False,
                )
            ]
        super().__init__(install_methods=install_methods)

    def version(self):
        # Centrifuger uses `-v` (not `--version`); passing `--version`
        # triggers a buffer overflow in the binary and aborts.
        try:
            result = subprocess.run(
                ['centrifuger', '-v'],
                capture_output=True,
                text=True,
                check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return 'unknown'
        return (result.stdout or result.stderr).strip() or 'unknown'

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())

    def execute(self, command, output=None, args=None, options=None):
        '''Run a Centrifuger command.'''
        args = args or []
        options = options or {}

        cmd = [command]
        for key, value in options.items():
            if value is None:
                cmd.append(key)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    cmd.extend([key, str(item)])
            else:
                cmd.extend([key, str(value)])
        cmd.extend(args)

        log.debug('Calling %s: %s', command, ' '.join(cmd))
        if output:
            with open(output, 'wt') as outf:
                subprocess.check_call(cmd, stdout=outf)
        else:
            subprocess.check_call(cmd)

    def build(self, db_prefix, taxonomy_tree, name_table, ref_fastas=None,
              ref_list=None, conversion_table=None, build_mem=None,
              num_threads=None):
        '''Build a Centrifuger index.'''
        if bool(ref_fastas) == bool(ref_list):
            raise ValueError('Specify exactly one of ref_fastas or ref_list')

        opts = {
            '--taxonomy-tree': taxonomy_tree,
            '--name-table': name_table,
            '-o': db_prefix,
        }
        if ref_fastas:
            if isinstance(ref_fastas, str):
                ref_fastas = [ref_fastas]
            opts['-r'] = ref_fastas
        else:
            opts['-l'] = ref_list
        if conversion_table:
            opts['--conversion-table'] = conversion_table
        if build_mem:
            opts['--build-mem'] = build_mem
        if num_threads:
            opts['-t'] = util_misc.sanitize_thread_count(num_threads)

        self.execute('centrifuger-build', options=opts)

    def classify(self, in_bam, db, output, k=1, unclassified_prefix=None,
                 classified_prefix=None, min_hitlen=None, hitk_factor=None,
                 merge_readpair=False, num_threads=None):
        '''Classify reads from an unaligned BAM file.'''
        if not os.path.isfile(in_bam):
            raise FileNotFoundError(in_bam)
        if samtools.SamtoolsTool().isEmpty(in_bam):
            log.warning('Input BAM is empty, skipping Centrifuger classification')
            with open(output, 'wt') as outf:
                outf.write('\t'.join(CLASSIFICATION_HEADER) + '\n')
            return

        tmp_fastq1 = util_file.mkstempfname('.1.fastq')
        tmp_fastq2 = util_file.mkstempfname('.2.fastq')
        tmp_fastq3 = util_file.mkstempfname('.s.fastq')

        try:
            picard_tool = picard.SamToFastqTool()
            picard_opts = {
                'CLIPPING_ATTRIBUTE': picard.SamToFastqTool.illumina_clipping_attribute,
                'CLIPPING_ACTION': 'X',
            }
            picard_tool.execute(
                in_bam,
                tmp_fastq1,
                tmp_fastq2,
                outFastq0=tmp_fastq3,
                picardOptions=picard.PicardTools.dict_to_picard_opts(picard_opts),
                JVMmemory=picard_tool.jvmMemDefault,
            )

            opts = {'-x': db}
            if num_threads:
                opts['-t'] = util_misc.sanitize_thread_count(num_threads)
            if k is not None:
                opts['-k'] = k
            if unclassified_prefix:
                opts['--un'] = unclassified_prefix
            if classified_prefix:
                opts['--cl'] = classified_prefix
            if min_hitlen is not None:
                opts['--min-hitlen'] = min_hitlen
            if hitk_factor is not None:
                opts['--hitk-factor'] = hitk_factor
            if merge_readpair:
                opts['--merge-readpair'] = None

            if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
                log.warning('running in single-end read mode')
                opts['-u'] = tmp_fastq3
            else:
                opts['-1'] = tmp_fastq1
                opts['-2'] = tmp_fastq2

            self.execute('centrifuger', output=output, options=opts)
        finally:
            for path in (tmp_fastq1, tmp_fastq2, tmp_fastq3):
                if os.path.exists(path):
                    os.unlink(path)

    def quant(self, db, classification, output, min_score=None,
              min_length=None, output_format=None):
        '''Run centrifuger-quant on a classification output.'''
        opts = {
            '-x': db,
            '-c': classification,
        }
        if min_score is not None:
            opts['--min-score'] = min_score
        if min_length is not None:
            opts['--min-length'] = min_length
        if output_format is not None:
            opts['--output-format'] = output_format

        self.execute('centrifuger-quant', output=output, options=opts)

    @staticmethod
    def classification_to_kraken2(classification, output):
        '''Convert Centrifuger classification TSV to Kraken2-style read TSV.

        This assumes the Centrifuger classification was produced with k=1.
        The output contains one classified Kraken2-style row per classified
        Centrifuger row and omits unclassified reads.
        '''
        previous_classified_read_id = None

        with util_file.open_or_gzopen(str(classification), 'rt') as inf, \
                util_file.open_or_gzopen(str(output), 'wt') as outf:
            header = inf.readline()
            if not header:
                return
            writer = csv.writer(outf, delimiter='\t', lineterminator='\n')

            for line_no, line in enumerate(inf, start=2):
                if not line.strip():
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 7:
                    raise ValueError(
                        'Malformed Centrifuger classification row '
                        f'at line {line_no}: expected at least 7 columns'
                    )

                read_id = cols[0]
                seq_id = cols[1]
                tax_id = cols[2]
                query_length = cols[6]

                if tax_id == '0' or seq_id == 'unclassified':
                    continue
                if read_id == previous_classified_read_id:
                    raise ValueError(
                        'Consecutive duplicate classified readID in Centrifuger output: '
                        f'{read_id}. Expected k=1 input.'
                    )

                previous_classified_read_id = read_id
                writer.writerow(['C', read_id, tax_id, query_length, 'N/A'])

    def kreport(self, db, classification, output, no_lca=False,
                show_zeros=False, is_count_table=False, min_score=None,
                min_length=None, report_score_data=False):
        '''Produce a Kraken-style report from a centrifuger classification.

        Wraps the centrifuger-kreport Perl script, which writes its output
        to stdout; the wrapper redirects stdout to `output`.

        Short-circuits on empty or header-only classification input
        (centrifuger-kreport exits 255 with `No sequence matches with
        given settings` in that case). This mirrors the empty-BAM
        short-circuit in classify() so callers can chain
        classify -> kreport unconditionally.
        '''
        if self._classification_is_empty(
                classification,
                has_header=not is_count_table):
            log.warning(
                'Classification file %s has no data rows, '
                'skipping centrifuger-kreport', classification,
            )
            with open(output, 'wt'):
                pass
            return

        opts = {'-x': db}
        if no_lca:
            opts['--no-lca'] = None
        if show_zeros:
            opts['--show-zeros'] = None
        if is_count_table:
            opts['--is-count-table'] = None
        if min_score is not None:
            opts['--min-score'] = min_score
        if min_length is not None:
            opts['--min-length'] = min_length
        if report_score_data:
            opts['--report-score-data'] = None

        self.execute(
            'centrifuger-kreport',
            output=output,
            args=[classification],
            options=opts,
        )

    @staticmethod
    def _classification_is_empty(path, has_header=True):
        '''Return True if `path` has no data rows.'''
        if os.path.getsize(path) == 0:
            return True
        with open(path, 'rt') as inf:
            if has_header:
                inf.readline()
            return not any(line.strip() for line in inf)

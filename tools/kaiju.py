'''
Kaiju - DNA-to-protein metagenomic classifier
'''
from builtins import super
import collections
import itertools
import logging
import os
import os.path
import shlex
import shutil
import subprocess
import tools

from Bio import SeqIO

import util.file

TOOL_VERSION = '1.6.3_yesimon'

log = logging.getLogger(__name__)


def read_a2t(fn, base_accession=True):
    if base_accession:
        accession_col = 0
    else:
        accession_col = 1
    d = {}
    with open(fn) as f:
        f.readline()  # Cannot use next(f) in python2
        for line in f.readlines():
            parts = line.split('\t')
            taxid = int(parts[2])
            accession = parts[accession_col]
            d[accession] = taxid
    return d


class Kaiju(tools.Tool):

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = [
                tools.CondaPackage("kaiju", version=TOOL_VERSION)
            ]
        super(Kaiju, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def build(self, db_prefix, protein_fastas, threads=None, options=None, option_string=None,
              tax_db=None, translate_accessions=False):
        '''Create a kaiju database.

        Args:
          db_prefix: Kaiju database prefix (file path not containing extension) to create.
          protein_fastas: List of input fasta files to process.
          tax_db: Contains accession2taxid used if translate_accessions is True.
          translate_accessions: If fasta IDs are accessions, translate to
            <accession>_<taxid> format kaiju expects.
        '''
        assert len(protein_fastas), ('Kaiju requires input files to create a database.')
        options = options or {}
        options['-a'] = 'protein'
        options['-n'] = util.misc.sanitize_thread_count(threads)
        options['-o'] = db_prefix

        temp_file = util.file.temp_catted_files(protein_fastas, prefix='kaiju_', suffix='.faa')

        db_fasta_fn = util.file.tempfname()
        with db_fasta_fn as db_fasta_fn:
            if translate_accessions:
                with temp_file as input_fasta:
                    with open(db_fasta_fn, 'w') as db_fasta:
                        prot2taxid_fn = os.path.join(tax_db, 'accession2taxid', 'prot.accession2taxid')
                        prot2taxid = read_a2t(prot2taxid_fn)
                        for record in SeqIO.parse(input_fasta, 'fasta'):
                            accession = record.id.split('.', 1)[0]
                            taxid = prot2taxid.get(accession)
                            if taxid is not None:
                                new_id = '_'.join([record.id, str(taxid)])
                                record.id = new_id
                                record.description = new_id
                            SeqIO.write(record, db_fasta, 'fasta')
                    option_string = db_fasta_fn
            else:
                option_string = temp_

            self.execute('mkbwt', options=options, option_string=option_string)
            self.execute('mkfmi', option_string=db_prefix)

    def classify(self, db, tax_db, in_bam, output_reads=None, output_report=None, rank='species', verbose=None, num_threads=None):
        assert output_reads or output_report
        output_ctx = util.file.tempfname()
        with util.file.tempfname() as temp_reads:
            # Fake loop to break out early if report is not wanted
            while True:
                tmp_fastq1 = util.file.mkstempfname('_1.fastq')
                tmp_fastq2 = util.file.mkstempfname('_2.fastq')
                picard = tools.picard.SamToFastqTool()
                picard_opts = {
                    'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
                    'CLIPPING_ACTION': 'X'
                }
                picard.execute(in_bam, tmp_fastq1, tmp_fastq2,
                            picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                            JVMmemory=picard.jvmMemDefault)

                nodes_dmp = os.path.join(tax_db, 'nodes.dmp')
                names_dmp = os.path.join(tax_db, 'names.dmp')

                opts = {
                    '-t': nodes_dmp,
                    '-f': db,
                    '-o': temp_reads,
                    '-z': util.misc.sanitize_thread_count(num_threads),
                    '-i': tmp_fastq1
                    }
                if verbose:
                    opts['-v'] = None

                if os.path.getsize(tmp_fastq2) >= 50:
                    opts['-j'] = tmp_fastq2

                self.execute('kaiju', options=opts)
                if not output_report:
                    break

                opts = {
                    '-r': rank,
                    '-p': True,
                    '-t': nodes_dmp,
                    '-n': names_dmp,
                    '-i': temp_reads,
                    '-o': output_report,
                    '-w': None,
                    '-x': None
                    }
                if verbose:
                    opts['-v'] = None
                self.execute('kaijuReport', options=opts)
                break
            if output_reads:
                subprocess.check_call(['pigz', '-9', temp_reads])
                shutil.move(temp_reads + '.gz', output_reads)

    def execute(self, command, options=None, option_string=None, return_stdout=False):
        '''Run a kaiju command

        Args:
          options: Dict of command line options to values. Set value to None
            for an option with no value.
          return_stdout: Whether to return stdout as well as in
            (exitcode, stdout).
        '''
        cmd = [command]
        if options:
            # We need some way to allow empty options args like --log, hence
            # we filter out on 'x is None'.
            cmd.extend([str(x) for x in itertools.chain(*options.items()) if x is not None])
        if option_string:
            cmd.extend(shlex.split(option_string))
        log.debug("Calling {}: {}".format(command, " ".join(cmd)))
        subprocess.check_call(cmd)

    def read_report(self, report_fn):
        report = collections.Counter()
        with open(report_fn) as f:
            f.readline()  # Cannot use next(f) in python2
            for line in f:
                if line.startswith('---'):
                    continue
                parts = line.strip().split('\t')
                percent = float(parts[0])
                reads = int(parts[1])
                fullname = parts[2]
                if 'cannot be assigned to a species' in fullname:
                    tax_id = 1
                elif 'unclassified' == fullname:
                    tax_id = 0
                else:
                    tax_id = parts[3]
                report[tax_id] = reads
        return report

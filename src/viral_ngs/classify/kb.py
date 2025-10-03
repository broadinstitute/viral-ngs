'''
kb_python classification tool
'''
import itertools
import logging
import os
import os.path
import shutil
import subprocess

import anndata

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc
from builtins import super

log = logging.getLogger(__name__)

class kb(tools.Tool):
    SUBCOMMANDS = ['count', 'ref', 'extract']

    def __init__(self, install_methods=None):
        if not install_methods:
            install_methods = []
            install_methods.append(tools.PrexistingUnixCommand(shutil.which('kb'), require_executability=False))
        super(kb, self).__init__(install_methods=install_methods)

    def version(self):
        return '1'

    @property
    def libexec(self):
        if not self.executable_path():
            self.install_and_get_path()
        return os.path.dirname(self.executable_path())


    def execute(self, command,output, args=None, options=None):
        '''Run a kb * command.

        Args:
          command: Subcommand to run.
          input: Input file kb_python/kallisto operates on.
          output: Output file to send to command.
          args: List of positional args.
          options: List of keyword options.
        '''
        options = options or {}

        if output:
            options['-o'] = output
        args = args or []

        cmd = [command]

        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
        cmd.extend(args)
        log.debug('Calling %s: %s', command, ' '.join(cmd))

        subprocess.check_call(cmd)
        
    def build(self, out_index, ref_fasta, aa=False, k=31, workflow_type='standard', num_threads=None):
        '''Create a kb_python index.
        Args:
          out_index: output index file
          ref_fasta: reference fasta file
          k: kmer size (default 31)
          aa: ref_fasta file contains amino acid sequences
          workflow_type: one of 'standard', 'nac', 'kite', 'custom'
        '''
        # build db
        build_opts = {
            '--threads': util.misc.sanitize_thread_count(num_threads)
        }
        if k:
            build_opts['-k'] = k
        if out_index:
            build_opts['-i'] = out_index
        if aa:
            build_opts['--aa'] = None
        if workflow_type:
            build_opts['--workflow-type'] = workflow_type
        self.execute('kb ref', None, None, args=[ref_fasta], options=build_opts)
        
    def classify(self, in_bam, index_file, out_dir, t2g_file, k=31, technology='bulk', h5ad=False, loom=False, num_threads=None):
        """Classify input reads (bam)

        Args:
          in_bam: unaligned reads
          index_file: kb_python index file
          out_dir: output directory
          t2g_file: transcript to gene mapping file
          k: kmer size (default 31)
          technology: one of '10xv2', '10xv3', '10xv3cr', '10xv3multiome', '10xv1', 'dropseq', 'indrop', 'sci-rna-seq', 'bulk'
          h5ad: output h5ad file
          loom: output loom file
          num_threads: number of threads to use
        """
        if tools.samtools.SamtoolsTool().isEmpty(in_bam):
            return

        opts = {
            '-i': index_file,
            '-g': t2g_file,
            '--threads': util.misc.sanitize_thread_count(num_threads)
        }
        if k:
            opts['-k'] = k
        if technology:
            opts['-x'] = technology
        if h5ad:
            opts['--h5ad'] = None
        if loom:
            opts['--loom'] = None   
            
            
        tmp_fastq1 = util.file.mkstempfname('.1.fastq')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq')
        tmp_fastq3 = util.file.mkstempfname('.s.fastq')
        # Do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
            log.warning("running in single-end read mode!")
            opts['--parity'] = "single"
            self.execute('kb count', out_dir, args=[tmp_fastq3], options=opts)
        else:
            opts['--parity'] = "paired"
            self.execute('kb count', out_dir, args=[tmp_fastq1, tmp_fastq2], options=opts)
            
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        os.unlink(tmp_fastq3)

    def extract(self, in_bam, index_file, target_ids, out_dir, t2g_file, aa=False, num_threads=None):
        """Extracts reads mapping to target ids from input reads (bam)
        
        Args:
          in_bam: unaligned read to extract reads from
          index_file: kb_python index file
          out_dir: output directory
          t2g_file: transcript to gene mapping file
          aa: ref_fasta file contains amino acid sequences
          target_ids: list of target ids to extract
          num_threads: number of threads to use
        """
        if tools.samtools.SamtoolsTool().isEmpty(in_bam):
            return

        opts = {
            '-i': index_file,
            '-g': t2g_file,
            '--kallisto': self.executable_path(),
            '--target-ids': ','.join(target_ids),
            '--threads': util.misc.sanitize_thread_count(num_threads)
        }
        if aa:
            opts['--aa'] = True
            
            
        tmp_fastq1 = util.file.mkstempfname('.1.fastq')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq')
        tmp_fastq3 = util.file.mkstempfname('.s.fastq')
        # Do not convert this to samtools bam2fq unless we can figure out how to replicate
        # the clipping functionality of Picard SamToFastq
        picard = tools.picard.SamToFastqTool()
        picard_opts = {
            'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
            'CLIPPING_ACTION': 'X'
        }
        picard.execute(in_bam, tmp_fastq1, tmp_fastq2, outFastq0=tmp_fastq3,
                       picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
                       JVMmemory=picard.jvmMemDefault)

        # Detect if input bam was paired by checking fastq 2
        if os.path.getsize(tmp_fastq2) < os.path.getsize(tmp_fastq3):
            opts['--parity'] = "single"
            self.execute('kb extract', out_dir, args=[tmp_fastq3], options=opts)
        else:
            ## TODO: Interleave the paired-end reads and pass through through 'kb extract'
            log.warning("kb extract does not support paired-end reads")
            self.execute('kb extract', out_dir, args=[tmp_fastq1], options=opts)
            
        os.unlink(tmp_fastq1)
        os.unlink(tmp_fastq2)
        os.unlink(tmp_fastq3)
        
    def merge_h5ads(self, in_h5ads, out_h5ad):
        """Merge multiple h5ad files into a single h5ad file
        
        Args:
          in_h5ads: list of input h5ad files
          out_h5ad: output h5ad file
        """
        assert len(in_h5ads) > 0, "no input h5ad files provided"
        if len(in_h5ads) == 1:
            shutil.copyfile(in_h5ads[0], out_h5ad)
            return
        
        adatas = [anndata.read_h5ad(f) for f in in_h5ads]
    
        # Check that all h5ad files have the same number of variables (genes)
        if len(adatas) > 1:
            n_vars_ref = adatas[0].n_vars
            var_names_ref = adatas[0].var_names
            for i, adata in enumerate(adatas[1:], 1):
                if adata.n_vars != n_vars_ref:
                    raise ValueError(f"Dimension mismatch: file {in_h5ads[0]} has {n_vars_ref} variables, "
                                f"but file {in_h5ads[i]} has {adata.n_vars} variables")
                if not adata.var_names.equals(var_names_ref):
                    raise ValueError(f"Variable names mismatch: file {in_h5ads[0]} and file {in_h5ads[i]} have different variable names")
                
        combined = anndata.concat(adatas, join='outer', axis=0, label='batch', fill_value=0)
        combined.write_h5ad(out_h5ad)
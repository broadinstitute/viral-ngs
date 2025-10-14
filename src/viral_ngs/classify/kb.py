'''
kb_python classification tool
'''
import itertools
import glob
import logging
import os
import os.path
import shutil
import subprocess

import anndata
import pandas as pd

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

        cmd = command.split()

        # We need some way to allow empty options args like --build, hence
        # we filter out on 'x is None'.
        cmd.extend([str(x) for x in itertools.chain(*options.items())
                    if x is not None])
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
            
        # Check for errors: kb_python sometimes catches exceptions without proper exit codes
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
        '''Create a kb_python index.
        Args:
          ref_fasta: reference fasta file
          index: output index file
          workflow: one of 'standard', 'nac', 'kite', 'custom'
          kmer_len: kmer size (default 31)
          protein: ref_fasta file contains amino acid sequences
        '''
        # build db
        build_opts = {
            '-t': util.misc.sanitize_thread_count(num_threads)
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
            '-t': util.misc.sanitize_thread_count(num_threads)
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

    def extract(self, in_bam, index_file, target_ids, out_dir, t2g_file, protein=False, num_threads=None):
        """Extracts reads mapping to target ids from input reads (bam)
        
        Args:
          in_bam: unaligned read to extract reads from (FASTQ or BAM)
          index_file: kb_python index file
          out_dir: output directory
          t2g_file: transcript to gene mapping file
          protein: ref_fasta file contains amino acid sequences
          target_ids: list of target ids to extract
          num_threads: number of threads to use
        """
        if tools.samtools.SamtoolsTool().isEmpty(in_bam):
            return

        opts = {
            '-i': index_file,
            '-g': t2g_file,
            '-ts': ','.join(target_ids),
            '-t': util.misc.sanitize_thread_count(num_threads)
        }
        if protein:
            opts['--aa'] = None
            
            
        tmp_fastq1 = util.file.mkstempfname('.1.fastq')
        tmp_fastq2 = util.file.mkstempfname('.2.fastq')
        tmp_fastq3 = util.file.mkstempfname('.s.fastq')
        tmp_interleaved = None
        try:
            if in_bam.lower().endswith('.fastq') or in_bam.lower().endswith('.fq') or in_bam.lower().endswith('.fastq.gz') or in_bam.lower().endswith('.fq.gz'):
                # Input is already FASTQ, use it directly (don't delete it later)
                self.execute('kb extract', out_dir, args=[in_bam], options=opts)
            else:
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
                    self.execute('kb extract', out_dir, args=[tmp_fastq3], options=opts)
                else:
                    tmp_interleaved = util.file.mkstempfname('.interleaved.fastq')
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
        except Exception as e:
            log.error("Error during kb extract: %s", e)
            raise
        finally:
            for path in (tmp_fastq1, tmp_fastq2, tmp_fastq3, tmp_interleaved):
                if path and os.path.exists(path):
                    try:
                        os.unlink(path)
                    except OSError as e:
                        log.warning("Failed to delete temporary file %s: %s", path, e)
        
        
    def merge_h5ads(self, in_count_tars, out_h5ad, tmp_dir_parent=None):
        """Merge multiple kb count output tarballs into a single h5ad file with sample metadata

        Args:
          in_count_tars: list of input kb count tarball files (tar.zst format)
          out_h5ad: output h5ad file path
          tmp_dir_parent: parent directory for temporary files (optional)
        """

        assert len(in_count_tars) > 0, "no input count tarballs provided"

        adatas = []

        for count_tar in in_count_tars:
            # Extract tarball to temporary directory
            with util.file.tmp_dir(dir=tmp_dir_parent) as tmp_dir:
                util.file.extract_tarball(count_tar, tmp_dir)

                # Find h5ad file in counts_unfiltered folder
                h5ad_files = glob.glob(os.path.join(tmp_dir, "counts_unfiltered", "*.h5ad"))
                assert len(h5ad_files) == 1, f"Expected exactly one .h5ad file in {count_tar}, found {len(h5ad_files)}"
                h5ad_file = h5ad_files[0]

                # Read sample name from matrix.cells file
                cells_file = os.path.join(tmp_dir, "matrix.cells")
                if os.path.exists(cells_file):
                    with open(cells_file, 'r') as f:
                        sample_name = f.read().strip()
                else:
                    # Fallback to tarball basename if matrix.cells doesn't exist
                    sample_name = os.path.splitext(os.path.splitext(os.path.basename(count_tar))[0])[0]
                    log.warning(f"matrix.cells not found in {count_tar}, using filename as sample name: {sample_name}")

                # Load h5ad file
                adata = anndata.read_h5ad(h5ad_file)

                # Add sample metadata to observations
                adata.obs['sample'] = sample_name
                adata.obs['batch_name'] = sample_name

                # Check if matrix.sample.barcodes exists and add barcode info
                barcodes_file = os.path.join(tmp_dir, "matrix.sample.barcodes")
                if os.path.exists(barcodes_file):
                    with open(barcodes_file, 'r') as f:
                        barcodes = [line.strip() for line in f]
                    if len(barcodes) == adata.n_obs:
                        adata.obs['batch_barcode'] = barcodes

                adatas.append(adata)

        # Handle single file case
        if len(adatas) == 1:
            log.warning("Only one count tarball provided - writing single file instead of merging")
            adatas[0].write_h5ad(out_h5ad)
            return

        # Check that all h5ad files have the same number of variables (genes)
        n_vars_ref = adatas[0].n_vars
        var_names_ref = adatas[0].var_names
        for i, adata in enumerate(adatas[1:], 1):
            if adata.n_vars != n_vars_ref:
                raise ValueError(f"Dimension mismatch: file {in_count_tars[0]} has {n_vars_ref} variables, "
                            f"but file {in_count_tars[i]} has {adata.n_vars} variables")
            if not adata.var_names.equals(var_names_ref):
                raise ValueError(f"Variable names mismatch: file {in_count_tars[0]} and file {in_count_tars[i]} have different variable names")

        # Merge all anndata objects
        combined = anndata.concat(adatas, join='outer', axis=0, label='batch', fill_value=0)
        combined.write_h5ad(out_h5ad)

    def parse_h5ad_counts(self, h5ad_file):
        """Parse h5ad file and return gene IDs with their total counts.

        Args:
          h5ad_file: path to h5ad file

        Returns:
          List of tuples (gene_id, count) sorted by gene ID
        """
        adata = anndata.read_h5ad(h5ad_file)

        counts_mtx = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

        gene_ids = adata.var.index.tolist()
        gene_totals = counts_mtx.sum(axis=0)

        if hasattr(gene_totals, 'A1'):  # numpy matrix
            gene_totals = gene_totals.A1

        return list(zip(gene_ids, gene_totals))

    def extract_hit_ids_from_h5ad(self, h5ad_file):
        """Parse h5ad file and extract all target IDs with 1 or more hits.
        
        Assumes h5ad contains a single sample (single row).

        Args:
          h5ad_file: path to h5ad file

        Returns:
          List of target IDs (strings) with counts > 0
        """
        adata = anndata.read_h5ad(h5ad_file)

        # Get count matrix and sum across all observations (rows)
        counts_mtx = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        gene_totals = counts_mtx.sum(axis=0)

        if hasattr(gene_totals, 'A1'):  # numpy matrix
            gene_totals = gene_totals.A1

        # Extract gene IDs and filter to those with counts > 0
        gene_ids = adata.var.index.tolist()
        hit_ids = [gene_id for gene_id, count in zip(gene_ids, gene_totals) if count > 0]

        return hit_ids

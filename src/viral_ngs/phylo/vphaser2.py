'''
    V-Phaser 2 variant caller
'''

import logging
import subprocess
import os
import tempfile
import shutil
import pysam
import tools
import util.file

log = logging.getLogger(__name__)

TOOL_NAME = "vphaser2"
TOOL_VERSION = "2.0"

class Vphaser2Tool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(Vphaser2Tool, self).__init__(install_methods=install_methods)

    def execute(self, inBam, outDir, numThreads=None):    # pylint: disable=W0221
        # index inBam
        samtoolsTool = tools.samtools.SamtoolsTool()
        samtoolsTool.index(inBam)

        log.debug("V-Phaser2 output directory: '%s'", outDir)

        cmd = [self.install_and_get_path(), '-i', inBam, '-o', outDir]
        cmd_str = ' '.join(cmd)
        envCopy = os.environ.copy()
        if numThreads is not None:
            envCopy['OMP_NUM_THREADS'] = str(numThreads)
            cmd_str = 'OMP_NUM_THREADS=%d ' % numThreads + cmd_str
        log.debug(cmd_str)

        # Use check_output instead of check_call so that we get error information
        #    if the executable can't run on GitHub Actions.
        # Also has the effect of suppressing informational messages from vphaser,
        #    which is probably a good thing.
        try:
            subprocess.check_output(cmd, env=envCopy, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as ex:
            log.error(ex.output if type(ex.output)==str else ex.output.decode('utf-8'))
            raise

    def iterate(self, inBam, numThreads=None):
        """
        Run V-Phaser 2 on inBam. Interate through lines in files
            CHROM.var.raw.txt in order of chroms in the inBam header.
        For each line yield:
         [CHROM, Ref_Pos, Var, Cons, Strd_bias_pval, Type, Var_perc,
          SNP_or_LP_Profile1, SNP_or_LP_Profile2, ...]
        """
        #outdir = tempfile.mkdtemp('vphaser2')
        with util.file.tmp_dir(prefix='vphaser2') as outdir:
            try:
                self.execute(inBam, outdir, numThreads)
            finally:
                # these V-Phaser droppings cause problems if they persist
                bti = inBam + '.bti'
                if os.path.isfile(bti):
                    os.unlink(bti)
            chromNames = pysam.Samfile(inBam).references
            for chromName in chromNames:
                outfile = os.path.join(outdir, chromName + '.var.raw.txt')
                if not os.path.exists(outfile):
                    continue
                with open(outfile, 'rt') as inf:
                    for line in inf:
                        if not line.startswith('#'):
                            yield [chromName] + line.strip().split()

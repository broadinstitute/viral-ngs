'''
    KMC - K-mer counter
'''

import logging
import os
import os.path
import subprocess
import shutil
import tempfile
import random

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = 'kmc'
TOOL_VERSION = '2.3.0'

log = logging.getLogger(__name__)

class KmcTool(tools.Tool):
    """Tool wrapper for KMC k-mer counter"""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION,
                                                  verifycmd='kmc 2>&1 | grep "K-Mer Counter" > /dev/null')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args, tool_sfx='', stdout=None):    # pylint: disable=W0221
        assert tool_sfx in ('', '_dump', '_tools')
        tool_cmd = [self.install_and_get_path()+tool_sfx] + args
        print('tool_cmd='+str(tool_cmd))
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout)
        if stdout:
            stdout.close()

    def _dbName(self, kmerDb):
        assert kmerDb.endswith('.kmc')
        return os.path.splitext(kmerDb)[0]

    def build_kmer_db(self, kmerDb, inFiles, kmerSize, kmerOpts='', maxMem=12, threads=1, kmc_opts=''):
        '''Build kmer database, from inFiles (which can be fasta or fastq).'''
        assert inFiles, "no input files!"
        if isinstance(inFiles,str): inFiles=[inFiles]
        assert len(set([ os.path.splitext(f)[1] for f in inFiles ]))==1, "mixing different input types not supported"
        with tempfile.TemporaryDirectory('_kmc_build') as kmc_tmp_dir:

            # convert any input bam files to fasta
            if inFiles[0].endswith('.bam'):
                inFilesNew=[]
                samtools=tools.samtools.SamtoolsTool()
                for f in inFiles:
                    reads1 = util.file.mkstempfname(suffix='.1.fa', directory=kmc_tmp_dir)
                    reads2 = util.file.mkstempfname(suffix='.2.fa', directory=kmc_tmp_dir)
                    reads0 = util.file.mkstempfname(suffix='.0.fa', directory=kmc_tmp_dir)
                    samtools.bam2fa(inBam=f, outFa1=reads1, outFa2=reads2, outFa0=reads0)
                    inFilesNew += [reads1, reads2, reads0]
                inFiles=inFilesNew

            log.debug('kmc_tmp_dir=' + kmc_tmp_dir)
            args=['-v', '-m'+str(maxMem), '-t'+str(threads), '-k'+str(kmerSize)]
            args += kmerOpts.split()
            args += kmc_opts.split()

            if not all([ f.endswith('.fq') or f.endswith('.fastq') or f.endswith('.fq.gz') or f.endswith('.fastq.gz')
                         for f in inFiles ]): 
                args += [ '-fm' ]
            with util.file.tempfname('.kmcInFiles.txt') as inFilesList:
                with tempfile.TemporaryDirectory('_kmcBuildDb') as kmcTmpDir:
                    if len(inFiles)>1:
                        with open(inFilesList, 'wt') as out:
                            out.write( '\n'.join(inFiles) + '\n' )
                        inFiles = '@'+inFilesList
                    else:
                        inFiles = inFiles[0]
                    args += [ inFiles, self._dbName(kmerDb), kmcTmpDir ]

                    self.execute( tool_sfx='', args = args )
                    util.file.touch_empty(kmerDb)

    def intersect(self, kmerDb1, kmerDb2, kmerDbOut, kmerDb1_opts='', kmerDb2_opts='', threads=1):
        '''Compute the intersection of two kmer sets'''
        args=['-v', '-t'+str(threads), 'intersect']
        args += [ self._dbName(kmerDb1) ] + kmerDb1_opts.split() + [ self._dbName(kmerDb2) ] + kmerDb2_opts.split() + [self._dbName(kmerDbOut)]
        self.execute( tool_sfx='_tools', args=args )

    def reduce(self, kmerDbIn, kmerDbOut, kmerDbIn_opts='', kmerDbOut_opts='', threads=1):
        '''Reduce a kmer set'''
        args=['-v', '-t'+str(threads), 'reduce']
        args += [ self._dbName(kmerDbIn) ] + kmerDbIn_opts.split() + [ self._dbName(kmerDbOut) ] + kmerDbOut_opts.split()
        self.execute( tool_sfx='_tools', args=args )


    def histogram(self, kmerDb, histogramFile, kmerDbOpts='-ci1'):
        '''Output kmer histogram'''
        args=['-v', 'histogram', self._dbName(kmerDb)]+kmerDbOpts.split()+[histogramFile]
        self.execute( tool_sfx='_tools', args=args )


    def dump(self, kmerDb, kmerListFile, kmerDbOpts='-ci1'):
        '''Dump kmers and counts from kmerDb to kmerListFile'''
        args = ['dump', self._dbName(kmerDb)] + kmerDbOpts.split() + [kmerListFile]
        self.execute( tool_sfx='_tools', args=args )
                
            






                      
                      

    


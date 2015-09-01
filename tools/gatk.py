'''
    GATK genotyping toolkit from the Broad Institute
    
    This software has different licenses depending on use cases.
    As such, we do not have an auto-downloader. The user must have GATK
    pre-installed on their own and available in $GATK_PATH.
'''

import tools, tools.picard, tools.samtools, util.file
import logging, os, os.path, subprocess, tempfile

log = logging.getLogger(__name__)

class GATKTool(tools.Tool) :
    jvmMemDefault = '2g'
    def __init__(self, path=None):
        self.tool_version = None
        install_methods = []
        for jarpath in [path, os.environ.get('GATK_PATH')]:
            if jarpath:
                if not jarpath.endswith('.jar'):
                    jarpath = os.path.join(jarpath, 'GenomeAnalysisTK.jar')
                install_methods.append(tools.PrexistingUnixCommand(
                    jarpath, verifycmd='java -jar %s --version' % jarpath,
                    verifycode=0, require_executability=False))
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def execute(self, command, gatkOptions=[], JVMmemory=None) :
        if JVMmemory==None:
            JVMmemory = self.jvmMemDefault
        toolCmd = ['java',
            '-Xmx' + JVMmemory,
            '-Djava.io.tmpdir=' + tempfile.tempdir,
            '-jar', self.install_and_get_path(),
            '-T', command] + list(map(str, gatkOptions))
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)
    
    def dict_to_gatk_opts(self, options) :
        return ["%s=%s" % (k,v) for k,v in options.items()]

    def version(self):
        if self.tool_version==None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        cmd = ['java', '-jar', self.install_and_get_path(), '--version']
        self.tool_version = subprocess.check_output(cmd).strip()
    
    def ug(self, inBam, refFasta, outVcf,
            options=["--min_base_quality_score", 15, "-ploidy", 4],
            JVMmemory=None, threads=1):

        if int(threads) < 1:
            threads = 1
        opts = ['-I', inBam, '-R', refFasta, '-o', outVcf,
            '-glm', 'BOTH',
            '--baq', 'OFF',
            '--useOriginalQualities',
            '-out_mode', 'EMIT_ALL_SITES',
            '-dt', 'NONE',
            '--num_threads', threads,
            '-stand_call_conf', 0,
            '-stand_emit_conf', 0,
            '-A', 'AlleleBalance',
        ]
        self.execute('UnifiedGenotyper', opts + options, JVMmemory=JVMmemory)

    def local_realign(self, inBam, refFasta, outBam, JVMmemory=None, threads=1):
        intervals = util.file.mkstempfname('.intervals')
        opts = ['-I', inBam, '-R', refFasta, '-o', intervals]
        self.execute('RealignerTargetCreator', opts, JVMmemory=JVMmemory)
        opts = ['-I', inBam, 
                '-R', refFasta, 
                '-targetIntervals', intervals, 
                '-o', outBam,
                '--num_threads', threads,
               ]
        self.execute('IndelRealigner', opts, JVMmemory=JVMmemory)
        os.unlink(intervals)

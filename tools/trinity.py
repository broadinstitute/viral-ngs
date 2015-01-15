'''
    The Trinity RNA-SEQ assembler
    
    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import logging, os, os.path, subprocess, tempfile, shutil
import tools, util.file

tool_version = "2011-11-26"
trinityVersion = "trinityrnaseq_r{}".format(tool_version)
url = "http://sourceforge.net/projects/trinityrnaseq/files/{}.tgz".format(
          trinityVersion)

log = logging.getLogger(__name__)


class TrinityTool(tools.Tool) :
    def __init__(self, install_methods = None) :
        if install_methods == None :
            install_methods = [
                DownloadAndBuildTrinity(url, trinityVersion + '/Trinity.pl')
                ]
        tools.Tool.__init__(self, install_methods = install_methods)
    
    def version(self) :
        return tool_version
    
    def execute(self, inFastq1, inFastq2, outFasta, min_contig_length=300):
        outdir = tempfile.mkdtemp(prefix='trinity-')
        cmd = [self.install_and_get_path(),
            '--CPU', '1',
            '--min_contig_length', str(min_contig_length),
            '--seqType', 'fq',
            '--left', inFastq1,
            '--right', inFastq2,
            '--output', outdir]
        log.debug(' '.join(cmd))
        subprocess.check_call(cmd)
        shutil.copyfile(os.path.join(outdir, 'Trinity.fasta'), outFasta)
        shutil.rmtree(outdir, ignore_errors=True)

class DownloadAndBuildTrinity(tools.DownloadPackage) :
    def post_download(self) :
        trinityDir = os.path.join(self.destination_dir, trinityVersion)
        if tool_version == "2011-11-26" :
            # Chrysalis doesn't compile. Need to add an include file.
            badFilePath = os.path.join(trinityDir, 'Chrysalis', 'analysis',
                                       'RunButterfly.cc')
            os.rename(badFilePath, badFilePath + '.orig')
            with open(badFilePath, 'wt') as outf:
                outf.write('#include <unistd.h>\n')
                with open(badFilePath+'.orig', 'rt') as inf:
                    for line in inf:
                        outf.write(line)
            
            # Trinity.pl insists on Java 1.6, but Java >= 1.6 is fine
            badFilePath = os.path.join(trinityDir, 'Trinity.pl')
            os.rename(badFilePath, badFilePath + '.orig')
            with open(badFilePath, 'wt') as outf:
                with open(badFilePath+'.orig', 'rt') as inf:
                    for line in inf:
                        if line == 'unless ($java_version =~ /java version \"1\.6\./) {\n':
                            outf.write('$java_version =~ /java version \"1\.(\d+)\./;\n')
                            outf.write('unless ($1 >= 6) {\n')
                        else:
                            outf.write(line)
            
        # Now we can make:
        os.system('cd "{}" && make -s'.format(trinityDir))
        shutil.rmtree(os.path.join(trinityDir, 'sample_data'), ignore_errors=True)
    
    def verify_install(self) :
        if not tools.DownloadPackage.verify_install(self) :
            return False
        # Verify that chrysalis and inchworm were built
        trinityDir = os.path.join(self.destination_dir, trinityVersion)
        chrysalisPath = os.path.join(trinityDir, 'Chrysalis', 'Chrysalis')
        inchwormPath  = os.path.join(trinityDir, 'Inchworm', 'src', 'inchworm')
        for path in [chrysalisPath, inchwormPath] :
            if not os.access(path, (os.X_OK | os.R_OK)) :
                log.debug('{} was not built.'.format(path))
                self.installed = False
        self.installed = True
        return self.installed

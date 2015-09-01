'''
    tbl2asn - from NCBI
    http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/
'''

__author__ = "dpark@broadinstitute.org"

import logging
import tools
import util.file
import os
import os.path
import subprocess
import gzip

url = 'ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/{os}.tbl2asn.gz'

log = logging.getLogger(__name__)


class Tbl2AsnTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [DownloadGzipBinary(url.format(os=get_bintype()), 'tbl2asn')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return None

    def execute(self, templateFile, inputDir, outputDir=None,
                source_quals=[], comment=None, verification='vb',
                file_type='s', structured_comment_file=None,
                per_genome_comment=False):

        toolCmd = [self.install_and_get_path(), '-t', templateFile]

        if inputDir:
            toolCmd += ['-p', inputDir]
        if outputDir:
            toolCmd += ['-r', outputDir]
        if source_quals:
            toolCmd.append('-j')
            toolCmd.append(' '.join("[{}={}]".format(k, v) for k, v in source_quals))
        if comment:
            toolCmd += ['-y', comment]
        if structured_comment_file:
            toolCmd += ['-w', structured_comment_file]
        if verification:
            toolCmd += ['-V', verification]
        if file_type:
            toolCmd += ['-a', file_type]
        if per_genome_comment:
            toolCmd += ['-X', 'C']

        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)


def get_bintype():
    uname = os.uname()
    if uname[0] == 'Darwin':
        return 'mac'
    elif uname[0] == 'Linux':
        if uname[4] == 'x86_64':
            return 'linux64'
        else:
            return 'linux'
    else:
        raise Exception('unsupported OS')


class DownloadGzipBinary(tools.DownloadPackage):

    def unpack(self, download_dir):
        util.file.mkdir_p(self.destination_dir)
        if (self.download_file.endswith('.gz') and
                not self.download_file.endswith('.tar.gz')):
            with gzip.open(os.path.join(download_dir, self.download_file)) as inf:
                with open(self.targetpath, 'wb') as outf:
                    outf.write(inf.read())
            os.chmod(self.targetpath, 0o755)

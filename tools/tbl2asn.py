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

TOOL_NAME = "tbl2asn"
TOOL_VERSION = "25.6" # quirk: versions error-out one year after their compilation date
TOOL_URL = 'ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/{os}.tbl2asn.gz'

log = logging.getLogger(__name__)


class Tbl2AsnTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
            install_methods.append(DownloadGzipBinary(TOOL_URL.format(os=get_bintype()), 'tbl2asn'))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return None

    # pylint: disable=W0221
    def execute(
        self,
        templateFile,
        inputDir,
        outputDir=None,
        source_quals=None,
        comment=None,
        verification='vb',
        file_type='s',
        structured_comment_file=None,
        per_genome_comment=False
    ):
        source_quals = source_quals or []

        tool_cmd = [self.install_and_get_path(), '-t', templateFile]

        if inputDir:
            tool_cmd += ['-p', inputDir]
        if outputDir:
            tool_cmd += ['-r', outputDir]
        if source_quals:
            tool_cmd.append('-j')
            tool_cmd.append(' '.join("[{}={}]".format(k, v) for k, v in source_quals))
        if comment:
            tool_cmd += ['-y', comment]
        if structured_comment_file:
            tool_cmd += ['-w', structured_comment_file]
        if verification:
            tool_cmd += ['-V', verification]
        if file_type:
            tool_cmd += ['-a', file_type]
        if per_genome_comment:
            tool_cmd += ['-X', 'C']

        log.debug(' '.join(tool_cmd))

        # tbl2asn has a fun quirk where if the build is more than a year old
        # it exits with a non-zero code and tells you to upgrade
        # See: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/demo/tbl2asn.c#L9674
        # We can try to work around this by examining the output for the upgrade message
        try:
            subprocess.check_output(tool_cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            old_version_expected_output = "This copy of tbl2asn is more than a year old.  Please download the current version."
            if old_version_expected_output in e.output.decode('UTF-8'):
                pass
            else:
                raise

    # pylint: enable=W0221


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
        if (self.download_file.endswith('.gz') and not self.download_file.endswith('.tar.gz')):
            with gzip.open(os.path.join(download_dir, self.download_file)) as inf:
                with open(self.targetpath, 'wb') as outf:
                    outf.write(inf.read())
            os.chmod(self.targetpath, 0o755)

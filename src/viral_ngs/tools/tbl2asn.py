'''
    tbl2asn - from NCBI
    http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/
'''

__author__ = "dpark@broadinstitute.org"

import logging, tools, util.file
import os, os.path, subprocess

url = 'ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/'

log = logging.getLogger(__name__)

class Tbl2AsnTool(tools.Tool):
    def __init__(self, install_methods = None):
        if install_methods == None:
            install_methods = []
        fname = "{}.tbl2asn.gz".format(get_bintype())
        raise NotImplementedError()
        tools.Tool.__init__(self, install_methods = install_methods)

    def version(self):
        return None

    def execute(self, args):
        pass

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





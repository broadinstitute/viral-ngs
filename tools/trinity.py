'''
    The Trinity RNA-SEQ assembler
    
    This uses an older version of Trinity that uses an older
    assembly algorithm that works better with highly diverse
    viral genomes.
'''

import logging, tools, util.file
import os, os.path, subprocess

tool_version = "2011-11-26"
url = "http://sourceforge.net/projects/trinityrnaseq/files/" \
    + "trinityrnaseq_r{}.tgz/download".format(tool_version)

log = logging.getLogger(__name__)




#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases.
'''

__author__ = "PLACEHOLDER"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)



if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

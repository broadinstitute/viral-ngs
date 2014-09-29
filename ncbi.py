#!/usr/bin/env python
'''This script contains a number of utilities for submitting our analyses
to NCBI's Genbank and SRA databases.

Requires python >= 2.7 and BioPython.  On the Broad cluster, it is known
to work with the Python-2.7 and Python-3.4 dotkits.'''

__author__ = "PLACEHOLDER"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging
import util.cmd, util.file, util.vcf, util.misc

log = logging.getLogger(__name__)
global_tool_paths = {}



if __name__ == '__main__':
	util.cmd.main_argparse(__commands__, __doc__)

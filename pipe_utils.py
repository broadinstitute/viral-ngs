#!/usr/bin/env python
"""
Utility functions to aid pipeline execution. 
"""

__author__ = "tomkinsc@broadinstitute.org"
__commands__ = []

import argparse
import logging
import os
import os.path
import re
import csv
import shutil
import subprocess
import tempfile
from collections import defaultdict, OrderedDict

import util.cmd
import util.file

log = logging.getLogger(__name__)







# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
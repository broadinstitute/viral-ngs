'''Utilities for working with the scripts directory'''

__author__ = "irwin@broadinstitute.org"

import os

def get_scripts_path() :
	'Return the path to the scripts directory.'
	return os.path.dirname(__file__)
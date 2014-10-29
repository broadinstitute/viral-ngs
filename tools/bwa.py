" tool for bwa "

__author__ = "hlevitin@broadinstitute.org"

import tools
import util.file
import os, logging

log = logging.getLogger(__name__)

# magic vars for now, later can set with config variables
# legacy version is version used in pipeline recipes (see old_scripts dir)
# current is lates version as of 8/27/2014
USE_CURRENT = True
DOWNLOAD_URL = {
	'legacy':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.2.tar.bz2',
	'current':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2'
	}

URL = DOWNLOAD_URL['current'] if USE_CURRENT else DOWNLOAD_URL['legacy']
BWA_DIR = '.'.join( [ x for x in URL.split("/")[-1].split('.') if
						x != "tar" and x != "bz2" and x != "gz"])

class Bwa(tools.Tool) :
	def __init__(self, install_methods = None) :
		log.debug("BWA_DIR: {}".format(BWA_DIR))
		if install_methods == None :
			install_methods = []
			install_methods.append( tools.DownloadPackage(
				URL, "{}/bwa".format(BWA_DIR),
				post_download_command="cd {}; make".format(BWA_DIR)))
			tools.Tool.__init__(self, install_methods = install_methods)

	def version(self) :
		return ''.join([c for c in BWA_DIR if c.isdigit() or c=='.'])

	def execute(self, subcommand, args=[], options={}, option_string="",
			post_cmd=""):
		"""
		args are required arguments for the specified bwa subcommand
			(order matters for bwa execution)
		options may be specified as key-value pairs of the form (flag: value)
			Leading dashes, ('-' or '--'), should be included in the key
			For flags without value arguments, value should equal the empty str
			(order does not matter for bwa execution)
		option_string spefifies options in a preformatted string.
			An alternative to options, but may be use in conjuction as well.
		post_cmd is appended to the end of the command.  It is intended to be
			used as a pipe ("| <other shell command>"), or to store output
			( "> output.sai")
		"""
		arg_str = " ".join(args)
		option_str = '{} {}'.format(' '.join([ "{} {}".format(k, v) for k, v in
												options.items() ]),
									option_string
									)
		cmd =  "{self.exec_path} {subcommand} {option_str} {arg_str} {post_cmd}" \
			.format(**locals())
		log.debug("Calling bwa with cmd: {}".format(cmd))
		return os.system(cmd)

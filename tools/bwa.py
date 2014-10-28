" tool for bwa "


__author__ = "hlevitin@broadinstitute.org"

import tools
import util.file
import os, logging

log = logging.getLogger(__name__)

# magic vars for now, later can set with config
# legacy version is version used in pipeline recipes (see old_scripts dir)
# current is lates version as of 8/27/2014
USE_CURRENT = True
BWA_URL = {
	'legacy':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.2.tar.bz2',
	'current':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2'
	}

URL = BWA_URL['current'] if USE_CURRENT else BWA_URL['legacy']
BWA_DIR = '.'.join( [ x for x in URL.split("/")[-1].split('.') if
						x != "tar" and x != "bz2" and x != "gz"])
log.debug("BWA_DIR: %s" % BWA_DIR)

class Bwa(tools.Tool) :
	def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods = []
			install_methods.append( tools.DownloadPackage(
				BWA_URL['current'], "%s/bwa" % BWA_DIR,
				post_download_command="cd %s; make" % (BWA_DIR, BWA_DIR)))
			tools.Tool.__init__(self, install_methods = install_methods)

	def version(self) :
		return ''.join([c for c in BWA_DIR if c.isdigit() or c=='.'])

	def execute(self, subcommand, args, options={}):
		arg_str = " ".join(args)
		option_str = " ".join([ "%s %s" % (k, v) for k, v in
								options.iteritems()])
		cmd = "./%s/bwa %s %s %s" % (BWA_DIR, subcommand, option_str, arg_str)
		os.system(cmd)

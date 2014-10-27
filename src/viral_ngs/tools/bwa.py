" tool for bwa "


__author__ = "hlevitin@broadinstitute.org"

import tools, util.file
import os



# magic vars for now, later can set with config
USE_CURRENT = True
BWA_URL = {
	'legacy':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.6.2.tar.bz2/download',
	'current':
	    'http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download'
	}

URL = BWA_URL['current'] if USE_CURRENT else BWA_URL['legacy']
BWA_DIR = '.'.join( [ x for x in URL.split("/")[-2].split('.') if
						x != "tar" and x != "bz2" and x != "gz"])


class Bwa(tools.Tool) :
    def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods.append( tools.DownloadPackage(
				BWA_URL['current'], "%s/bwa" % BWA_DIR,
				post_download_command = 'cd %s; make' % BWA_DIR))
			tools.Tool.__init__(self, install_methods = install_methods)

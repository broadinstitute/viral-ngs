"tools.Tool for trimmomatic."

import tools, util.file
import os, tempfile


trimmomaticBroadUnixPath = '/idi/sabeti-scratch/kandersen/bin/trimmomatic/trimmomatic-0.32.jar'

trimmomaticURL = 'http://www.usadellab.org/cms/uploads/supplementary/' \
				 'Trimmomatic/Trimmomatic-0.32.zip'

class TrimmomaticTool(tools.Tool) :
	def __init__(self, install_methods = None) :
		if install_methods == None :
			install_methods = []
			install_methods.append(tools.PrexistingUnixCommand(trimmomaticBroadUnixPath))
			install_methods.append(DownloadAndBuildTrimmomatic())
		tools.Tool.__init__(self, install_methods = install_methods)

class DownloadAndBuildTrimmomatic(tools.DownloadPackage) :
	def __init__(self) :
		buildDir = util.file.get_build_path()
		targetpath = os.path.join(buildDir, 'Trimmomatic-0.32', 'trimmomatic-0.32.jar')
		download_dir = tempfile.tempdir
		unpack_dir = buildDir
		tools.DownloadPackage.__init__(self, url = trimmomaticURL, targetpath = targetpath,
									   download_dir = download_dir, unpack_dir = unpack_dir,
									   requireExecutability = False)

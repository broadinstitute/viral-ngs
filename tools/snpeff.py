''' snpEff - a tool for annotating genetic consequences of variants in VCF format

	http://snpeff.sourceforge.net/
'''

import tools
import util.vcf, util.file
import os, tempfile

class SnpEff(tools.Tool):
	''' 
	'''
	def __init__(self, install_methods=None):
		if install_methods==None:
			install_methods = [BroadUnix(), DownloadAndConfigJar()]
		super(SnpEff, self).__init__(install_methods = install_methods)
	def version(self):
		return "4.0"
	def execute(self, inVcf, outVcf, commands, java_flags='-Xmx2g'):
		raise NotImplementedError  #### TO DO
		if outVcf.endswith('.gz'):
			tmpVcf = util.file.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
		else:
			tmpVcf = outVcf
		cmdline = "cat %s | java -Xmx2g -Djava.io.tmpdir=%s -jar %s/snpEff.jar eff -c %s/snpEff.config -o %s %s -treatAllAsProteinCoding false -noLog -ud 0 -noStats > %s" % (
			inVcf, tempfile.tempdir, snpEffPath, snpEffPath, (vcf_out and 'vcf' or 'txt'), genome, tmpVcf)
		if inVcf.endswith('.gz'):
			cmdline = 'z'+cmdline
		log.info("vcf_snpEff: %s" % cmdline)
		assert not os.system(cmdline)
		if outVcf.endswith('.gz'):
			if vcf_out:
				util.vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=tabixPath)
			else:
				cmdline = "%s/bgzip -c %s > %s" % (tabixPath, tmpVcf, outVcf)
				assert not os.system(cmdline)
				cmdline = "%s/tabix %s -f -p vcf" % (tabixPath, outVcf)
				assert not os.system(cmdline)
			os.unlink(tmpVcf)


class BroadUnix(tools.PrexistingUnixCommand):
	def __init__(self, path='/idi/sabeti-data/software/snpEff/snpEff_3.6-dev'):
		super(BroadUnix, self).__init__(path=path,
			verifycmd='java -Xmx50M -jar %s/snpEff.jar -h' % path)

class DownloadAndConfigJar(tools.DownloadPackage):
	def __init__(self, url='http://downloads.sourceforge.net/project/snpeff/snpEff_v4_0_core.zip',
		targetpath='build/snpEff', download_dir=tempfile.tempdir, unpack_dir='build'):
		super(DownloadAndConfigJar, self).__init__(url=url, targetpath=targetpath,
			download_dir=download_dir, unpack_dir=unpack_dir,
			verifycmd='java -Xmx50M -jar %s/snpEff.jar -h' % targetpath)
	def post_download(self):
		self.unpack()
		# other stuff here to set up config file and d/l some genomes
		#### TO DO



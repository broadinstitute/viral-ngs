''' snpEff - a tool for annotating genetic consequences of variants in VCF format

	http://snpeff.sourceforge.net/
'''

import tool, util.vcf, util.file

class SnpEff(tool.Tool):
	''' 
	'''
	def __init__(self, install_methods = [BroadUnix(), DownloadJar()]):
		super(SnpEff, self).__init__(install_methods = install_methods)
	def version(self):
		return "4.0"
	def execute(self, args):
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


class BroadUnix(tool.PrexistingUnixCommand):
	def __init__(self, path='/idi/sabeti-data/software/snpEff/snpEff_3.6-dev'):
		super(BroadUnix, self).__init__(path=path,
			verifycmd='java -Xmx50M -jar %s/snpEff.jar -h' % path)

class DownloadJar(tool.InstallMethod):
	''' #### TO DO
	'''
	def __init__(self):
		pass
		
	def is_attempted(self):
		pass
	def attempt_install(self):

		#	http://downloads.sourceforge.net/project/snpeff/snpEff_v4_0_core.zip

		pass
	def is_installed(self):
		pass
	def executable_path(self):
		pass





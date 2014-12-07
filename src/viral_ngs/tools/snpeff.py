'''
snpEff - a tool for annotating genetic consequences of variants in VCF format
http://snpeff.sourceforge.net/
'''

import pysam
import tools
import util.vcf, util.file
import os, tempfile, logging

log = logging.getLogger(__name__)

class SnpEff(tools.Tool):
    def __init__(self, install_methods=None, install_genomes=None):
        if install_methods==None:
            install_methods = [BroadUnix(), DownloadAndConfigJar()]
        if install_genomes==None:
            install_genomes = [SnpEffGenome('zebov.sl',
                'Zaire ebolavirus Sierra Leone G3686.1',
                ['KM034562.1'],
                ['http://www.ncbi.nlm.nih.gov/nuccore/661348725'])]
        self.install_genomes = install_genomes
        super(SnpEff, self).__init__(install_methods = install_methods)

    def version(self):
        return "4.0"

    def install(self):
        super(SnpEff, self).install()
        for g in self.install_genomes:
            g.install_genome(self)

    def has_genome(self, genome):
        pass ### TODO

    def execute(self, args, java_flags='-Xmx2g', pre_pipe='', post_pipe=''):
        cmdline = ' '.join([
            prepipe,
            'java', java_flags,
            '-Djava.io.tmpdir={}'.format(tempfile.tempdir),
            '-jar', '{}/snpEff.jar'.format(self.executable_path),
            args,
            post_pipe
            ])
        os.system(cmdline)

    def eff_vcf(self, inVcf, outVcf, genome, java_flags='-Xmx2g',
            in_format='vcf', out_format='vcf', eff_options=''):
        """
        TODO: docstring here
        """
        if outVcf.endswith('.vcf.gz'):
            tmpVcf = util.file.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
        else:
            tmpVcf = outVcf

        args = ' '.join([
                'eff',
                    '-c', '{}/snpEff.config'.format(self.executable_path()),
                    '-i', in_format,
                    '-o', out_format,
                    genome,
                    '-treatAllAsProteinCoding false',
                    '-noLog',
                    '-ud 0',
                    '-noStats',
                    eff_options
                ])

        if inVcf.endswith('.gz'):
            pre_pipe = "zcat {} | ".format(inVcf)
        else:
            pre_pipe = "cat {} | ".format(inVcf)
        post_pipe = " > {}".format(tmpVcf)
        self.execute(args, java_flags=java_flags, pre_pipe=pre_pipe,
                post_pipe=post_pipe)
        
        if outVcf.endswith('.vcf.gz'):
            pysam.tabix_compress(tmpVcf, outVcf, force=True)
            pysam.tabix_index(outVcf, force=True, preset='vcf')
            os.unlink(tmpVcf)


class BroadUnix(tools.PrexistingUnixCommand):
    """
    commented out in init, remove?:
    verifycmd=
        'java -Xmx50M -jar %s/snpEff.jar -h -noLog &> /dev/null' % path,
    verifycode=65280 xxx unfortunately snpEff is not returning a meaningful
        code)
    """
    def __init__(self, path='/idi/sabeti-data/software/snpEff/snpEff_3.6-dev'):
        super(BroadUnix, self).__init__(path=path)

    def verify_install(self):
        super(BroadUnix, self).verify_install()


class DownloadAndConfigJar(tools.DownloadPackage):
    """
    commented out in init, remove?:
    verifycmd=
        'java -Xmx50M -jar %s/snpEff.jar -h -noLog &> /dev/null' % path,
    verifycode=65280 xxx unfortunately snpEff is not returning a meaningful
        code)
    """

    DEFAULT_URL = \
        'http://downloads.sourceforge.net/project/snpeff/snpEff_v4_0_core.zip'

    def __init__(self, url=DEFAULT_URL, target_rel_path='snpEff',
            destination_dir=None):
        super(DownloadAndConfigJar, self).__init__(url=url,
                target_rel_path=target_rel_path,
                destination_dir=destination_dir)


class SnpEffGenome:
    def __init__(self, name, desc, chroms, chrom_urls, codonTableMap={}):
        pass
    def install_genome(self, tool):
        pass # check that it doesn't exist already
        #zebov.k.genome : Zaire ebolavirus Kissidougou
        #zebov.k.reference : http://www.ncbi.nlm.nih.gov/nuccore/KJ660346.1
        #vibrio.genome : Vibrio Cholerae
        #   vibrio.chromosomes : NC_002505.1, NC_002506.1
        #   vibrio.NC_002505.1.codonTable : Bacterial_and_Plant_Plastid
        #   vibrio.NC_002506.1.codonTable : Bacterial_and_Plant_Plastid

        # cd snpEffpath
        # java -jar snpEff.jar build -genbank -v genomeid

'''
snpEff - a tool for annotating genetic consequences of variants in VCF format
http://snpeff.sourceforge.net/
'''

import pysam
import tools
import util.vcf, util.file
import os, tempfile, logging, subprocess

log = logging.getLogger(__name__)

URL = 'http://downloads.sourceforge.net/project/snpeff/snpEff_v4_1_core.zip'

class SnpEff(tools.Tool):
    jvmMemDefault = '4g'

    def __init__(self, install_methods=None, install_genomes=None):
        if not install_methods:
            install_methods = [tools.DownloadPackage(
                URL, 'snpEff/snpEff.jar', require_executability=False)]
        if not install_genomes:
            install_genomes = [SnpEffGenome('zebov.sl',
                'Zaire ebolavirus Sierra Leone G3686.1',
                ['KM034562.1'],
                ['http://www.ncbi.nlm.nih.gov/nuccore/661348725'])]
        self.install_genomes = install_genomes
        self.known_dbs = set()
        self.installed_dbs = set()
        super(SnpEff, self).__init__(install_methods = install_methods)

    def version(self):
        return "4.1"

    def install(self):
        super(SnpEff, self).install()
        for g in self.install_genomes:
            g.install_genome(self)

    def execute(self, command, args, JVMmemory=None, stdin=None, stdout=None):
        if JVMmemory==None:
            JVMmemory = self.jvmMemDefault
        toolCmd = ['java',
            '-Xmx' + JVMmemory,
            '-Djava.io.tmpdir=' + tempfile.tempdir,
            '-jar', self.install_and_get_path(),
            command] + args
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd, stdin=stdin, stdout=stdout)

    def has_genome(self, genome):
        if not self.known_dbs:
            for row in self.available_databases():
                pass
        return genome in self.installed_dbs

    def download_db(self, dbname, verbose=False):
        opts = [dbname]
        if verbose:
            opts.append('-v')
        self.execute('download', opts)
        self.known_dbs.add(dbname)
        self.installed_dbs.add(dbname)
    
    def available_databases(self):
        toolCmd = ['java', '-jar', self.install_and_get_path(), 'databases']
        split_points = []
        keys = ['Genome', 'Organism', 'Status', 'Bundle', 'Database']
        self.installed_dbs = set()
        self.known_dbs = set()
        for line in subprocess.check_output(toolCmd, universal_newlines=True).split('\n'):
            line = line.strip()
            if not split_points:
                if not line.startswith('Genome'):
                    raise Exception()
                split_points = list(line.index(key) for key in keys)
            elif not line.startswith('----'):
                indexes = split_points + [len(line)]
                row = dict((keys[i], line[indexes[i]:indexes[i+1]].strip()) for i in range(len(split_points)))
                self.known_dbs.add(row['Genome'])
                if row.get('Status')=='OK':
                    self.installed_dbs.add(row['Genome'])
                yield row

    def annotate_vcf(self, inVcf, genome, outVcf, JVMmemory=None):
        """
        TODO: docstring here
        """
        if outVcf.endswith('.vcf.gz'):
            tmpVcf = util.file.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
        elif outVcf.endswith('.vcf'):
            tmpVcf = outVcf
        else:
            raise Exception("invalid input")

        args = [
            genome,
            '-treatAllAsProteinCoding', 'false',
            '-t',
            '-noLog',
            '-ud', '0',
            '-noStats'
            ]
        with open(tmpVcf, 'wt') as outf:
            self.execute('ann', args, JVMmemory=JVMmemory, stdout=outf)
        
        if outVcf.endswith('.vcf.gz'):
            pysam.tabix_compress(tmpVcf, outVcf, force=True)
            pysam.tabix_index(outVcf, force=True, preset='vcf')
            os.unlink(tmpVcf)



class SnpEffGenome:
    def __init__(self, id, desc, chroms, data_dir, build_opts='', codonTableMap={}):
        self.id=id
        self.desc=desc
        self.chroms=chroms
        self.data_dir=data_dir
        self.build_opts=build_opts
        self.codonTableMap=codonTableMap
    def has_genome(self):
        pass
    def install_genome(self, tool):
        if self.has_genome():
            return
        pass # check that it doesn't exist already
        #zebov.k.genome : Zaire ebolavirus Kissidougou
        #zebov.k.reference : http://www.ncbi.nlm.nih.gov/nuccore/KJ660346.1
        #vibrio.genome : Vibrio Cholerae
        #   vibrio.chromosomes : NC_002505.1, NC_002506.1
        #   vibrio.NC_002505.1.codonTable : Bacterial_and_Plant_Plastid
        #   vibrio.NC_002506.1.codonTable : Bacterial_and_Plant_Plastid
        
        # cd snpEffpath
        # java -jar snpEff.jar build -genbank -v genomeid

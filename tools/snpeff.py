'''
snpEff - a tool for annotating genetic consequences of variants in VCF format
http://snpeff.sourceforge.net/
'''

# built-ins
import hashlib
import tempfile
import os, tempfile, logging, subprocess

# third-party
import pysam

# module-specific
import tools, util.file, util.genbank

log = logging.getLogger(__name__)

URL = 'http://downloads.sourceforge.net/project/snpeff/snpEff_v4_1_core.zip'

class SnpEff(tools.Tool):
    jvmMemDefault = '4g'

    def __init__(self, install_methods=None, extra_genomes=['KJ660346.2']):
        if not install_methods:
            install_methods = [DownloadAndTweakSnpEff(URL, extra_genomes)]
        self.known_dbs = set()
        self.installed_dbs = set()
        super(SnpEff, self).__init__(install_methods = install_methods)

    def version(self):
        return "4.1"

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

    def create_db(self, accessions, emailAddress, JVMmemory):
        sortedAccessionString = ", ".join(sorted(accessions))
        databaseId = hashlib.sha256(sortedAccessionString).hexdigest()

        # if the database is not installed, we need to make it
        if not self.has_genome(databaseId):
            config_file = os.path.join(os.path.dirname(self.install_and_get_path()), 'snpEff.config')
            dataDir = get_data_dir(config_file)

            # if the data directory specified in the config is absolute, use it
            # otherwise get the data directory relative to the location of the config file
            if os.path.isabs(dataDir):
                outputDir = os.path.join(dataDir, databaseId)
            else:
                outputDir = os.path.realpath(os.path.join(os.path.dirname(config_file), dataDir, databaseId))

            #tempDir = tempfile.gettempdir()
            records = util.genbank.fetch_full_records_from_genbank(accessions, outputDir, emailAddress, forceOverwrite=True, combinedFilePrefix="genes", removeSeparateFiles=False)
            combinedGenbankFilepath = records[0]

            add_genomes_to_snpeff_config_file(config_file, [(databaseId, sortedAccessionString, sortedAccessionString)])
            self.known_dbs.add(databaseId)
            self.installed_dbs.add(databaseId)

            args = [
                '-genbank',
                '-v', databaseId
                ]
            self.execute('build', args, JVMmemory=JVMmemory)         

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

    def annotate_vcf(self, inVcf, genomes, outVcf, emailAddress, JVMmemory=None):
        """
        Annotate variants in VCF file with translation consequences using snpEff.
        """
        if outVcf.endswith('.vcf.gz'):
            tmpVcf = util.file.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
        elif outVcf.endswith('.vcf'):
            tmpVcf = outVcf
        else:
            raise Exception("invalid input")

        sortedAccessionString = ", ".join(sorted(genomes))
        databaseId = hashlib.sha256(sortedAccessionString).hexdigest()


        genomeToUse = ""
        # backward compatability for where a single genome name is provided
        if self.has_genome(genomes[0]):
            genomeToUse = genomes[0]

        # if the hash of the accessions passed in is not present in the genomes db
        if not self.has_genome(databaseId):
            self.create_db(genomes, emailAddress, JVMmemory)

        if not genomeToUse and self.has_genome(databaseId):
            genomeToUse = databaseId
        else:
            raise Exception()
        
        args = [
            '-treatAllAsProteinCoding', 'false',
            '-t',
            '-noLog',
            '-ud', '0',
            '-noStats',
            '-noShiftHgvs',
            genomeToUse,
            inVcf
            ]

        with open(tmpVcf, 'wt') as outf:
            self.execute('ann', args, JVMmemory=JVMmemory, stdout=outf)
        
        if outVcf.endswith('.vcf.gz'):
            pysam.tabix_compress(tmpVcf, outVcf, force=True)
            pysam.tabix_index(outVcf, force=True, preset='vcf')
            os.unlink(tmpVcf)


def get_data_dir(config_file):
    dataDir = ""
    with open(config_file, 'rt') as inf:
        for line in inf:
            if line.strip().startswith('data.dir'):
                dataDir = line[line.find("=")+1:].strip()
                break
    return dataDir

def add_genomes_to_snpeff_config_file(config_file, new_genomes):
    """
        new_genomes is a 3-tuple (g,d,c):
            where g is the genome name
                  d is the description
                  c is a comma-separated list of chromosomes
    """
    genomes = set()
    with open(config_file, 'rt') as inf:
        for line in inf:
            if not line.startswith('#') and line.strip():
                i = line.find('.genome : ')
                if i>=0:
                    genomes.add(line[:i])
    with open(config_file, 'at') as outf:
        for (g, d, c) in new_genomes:
            if g not in genomes:
                outf.write('{}.genome : {}\n'.format(g, d))
                if g != c:
                    outf.write("\t{}.chromosomes : {}\n".format(g, c))

class DownloadAndTweakSnpEff(tools.DownloadPackage):
    def __init__(self, url, extra_genomes=[]):
        self.extra_genomes = extra_genomes
        super(DownloadAndTweakSnpEff, self).__init__(
            url, 'snpEff/snpEff.jar', require_executability=False)
    def post_download(self):
        config_file = os.path.join(self.destination_dir, 'snpEff', 'snpEff.config')
        add_genomes_to_snpeff_config_file(config_file, zip(self.extra_genomes, self.extra_genomes, self.extra_genomes))

'''
snpEff - a tool for annotating genetic consequences of variants in VCF format
http://snpeff.sourceforge.net/
'''

# built-ins
import hashlib
import os
import tempfile
import logging
import subprocess
import shutil

# third-party
import pysam

# module-specific
import tools
import util.file
import util.misc
import phylo.genbank

_log = logging.getLogger(__name__)

TOOL_NAME = 'snpeff'
TOOL_VERSION = '4.3.1t'

URL = 'http://downloads.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip'


class SnpEff(tools.Tool):

    def __init__(self, install_methods=None, extra_genomes=None):
        self.jvmMemDefault = '4g'
        extra_genomes = extra_genomes or ['KJ660346.2']
        if not install_methods:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable="snpEff", version=TOOL_VERSION))
            install_methods.append(DownloadAndTweakSnpEff(URL, extra_genomes))
        self.known_dbs = set()
        self.installed_dbs = set()
        super(SnpEff, self).__init__(install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, JVMmemory=None, stdin=None, stdout=None, stderr=None):    # pylint: disable=W0221
        if not JVMmemory:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        if self.install_and_get_path().endswith(".jar"):
            tool_cmd = [
                'java', '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-jar', self.install_and_get_path(),
                command
            ] + args
        else:
            tool_cmd = [
                self.install_and_get_path(), '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), command
            ] + args

        _log.debug(' '.join(tool_cmd))
        return util.misc.run_and_print(tool_cmd, stdin=stdin, stderr=stderr, buffered=True, silent=command in ("databases","build"), check=True)

    def has_genome(self, genome):
        if not self.known_dbs:
            for _ in self.available_databases():
                pass

        return genome in self.installed_dbs

    def download_db(self, dbname, verbose=False):
        opts = [dbname]
        if verbose:
            opts.append('-v')
        self.execute('download', opts)
        self.known_dbs.add(dbname)
        self.installed_dbs.add(dbname)

    def create_db(self, accessions, emailAddress=None, JVMmemory=None):
        sortedAccessionString = ", ".join([phylo.genbank.parse_accession_str(acc) for acc in sorted(accessions)])
        databaseId = hashlib.sha256(sortedAccessionString.encode('utf-8')).hexdigest()[:55]

        # if the database is not installed, we need to make it
        if not self.has_genome(databaseId):
            config_file = os.path.join(os.path.dirname(os.path.realpath(self.install_and_get_path())), 'snpEff.config')
            data_dir = get_data_dir(config_file)

            # if the data directory specified in the config is absolute, use it
            # otherwise get the data directory relative to the location of the config file
            if os.path.isabs(data_dir):
                outputDir = os.path.join(data_dir, databaseId)
            else:
                outputDir = os.path.realpath(os.path.join(os.path.dirname(config_file), data_dir, databaseId))

            phylo.genbank.fetch_full_records_from_genbank(
                sorted(accessions), 
                outputDir,
                emailAddress,
                forceOverwrite=True,
                combinedFilePrefix="genes",
                removeSeparateFiles=False
            )

            # create a temp config file in the same location as the original
            # since the data dir for the built database has a relative
            # location by default
            tmp_config_file = config_file+".temp"

            # build a new snpEff database using the downloaded genbank file
            # Note that if the snpEff command fails for some reason,
            # we should take care to not include an entry for the database
            # in the config file, which is the record of installed databases.
            # In the event the build fails, the config file must reflect
            # databases present in the data_dir, otherwise calls to execute snpEff will
            # fail if the data_dir lacks databases for entries listed in the config file.

            # make a temp copy of the old config file
            shutil.copyfile(config_file, tmp_config_file)

            # add the new genome to the temp config file
            # since snpEff will not attempt to build unless the
            # database is in the config file
            add_genomes_to_snpeff_config_file(tmp_config_file, [(databaseId, sortedAccessionString, sortedAccessionString)])

            try:
                args = ['-genbank', '-v', databaseId, "-c", tmp_config_file]
                self.execute('build', args, JVMmemory=JVMmemory)
            except:
                # remove temp config file if the database build failed
                shutil.unlink(tmp_config_file)
                raise

            # copy the temp config including the built database 
            # if the execute('build') command did not raise an exception
            shutil.move(tmp_config_file, config_file)
            self.known_dbs.add(databaseId)
            self.installed_dbs.add(databaseId)

    def available_databases(self):
        # do not capture stderr, since snpEff writes 'Picked up _JAVA_OPTIONS'
        # which is not helpful for reading the stdout of the databases command
        with open(os.devnull, "wb") as devnull:
            command_ps = self.execute("databases", args=[], stderr=devnull)

            split_points = []
            keys = ['Genome', 'Organism', 'Status', 'Bundle', 'Database']
            self.installed_dbs = set()
            self.known_dbs = set()
            for line in command_ps.stdout.decode("utf-8").splitlines():
                line = line.strip()
                if not split_points:
                    if not line.startswith('Genome'):
                        raise Exception()
                    split_points = list(line.index(key) for key in keys)
                elif not line.startswith('----'):
                    indexes = split_points + [len(line)]
                    row = dict((keys[i], line[indexes[i]:indexes[i + 1]].strip()) for i in range(len(split_points)))
                    self.known_dbs.add(row['Genome'])
                    if row.get('Status') == 'OK':
                        self.installed_dbs.add(row['Genome'])
                    yield row

    def annotate_vcf(self, inVcf, genomes, outVcf, emailAddress=None, JVMmemory=None):
        """
        Annotate variants in VCF file with translation consequences using snpEff.
        """
        if outVcf.endswith('.vcf.gz'):
            tmpVcf = util.file.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
        elif outVcf.endswith('.vcf'):
            tmpVcf = outVcf
        else:
            raise Exception("invalid input")

        sortedAccessionString = ", ".join([phylo.genbank.parse_accession_str(acc) for acc in sorted(genomes)])
        databaseId = hashlib.sha256(sortedAccessionString.encode('utf-8')).hexdigest()[:55]

        genomeToUse = ""

        # if we don't have the genome, by name (snpEff official) or by hash (custom)
        if (not self.has_genome(databaseId)):
            if (not self.has_genome(genomes[0])):
                _log.info("Checking for snpEff database online...")
                # check to see if it is available for download, and if so install it
                for row in self.available_databases():
                    if (genomes[0].lower() in row['Genome'].lower()) or (
                        genomes[0].lower() in row['Bundle'].lower()
                    ) or (
                        genomes[0].lower() in row['Organism'].lower()
                    ):
                        self.download_db(row['Genome'])

        # backward compatability for where a single genome name is provided
        if self.has_genome(genomes[0]):
            genomeToUse = genomes[0]
        else:
            # if the hash of the accessions passed in is not present in the genomes db
            if not self.has_genome(databaseId):
                self.create_db(genomes, emailAddress, JVMmemory)

            if self.has_genome(databaseId):
                genomeToUse = databaseId

        if not genomeToUse:
            raise Exception()

        args = [
            '-treatAllAsProteinCoding', 'false', '-t', '-noLog', '-ud', '0', '-noStats', '-noShiftHgvs', genomeToUse,
            os.path.realpath(inVcf)
        ]

        command_ps = self.execute('ann', args, JVMmemory=JVMmemory)
        if command_ps.returncode == 0:
            with open(tmpVcf, 'wt') as outf:
               outf.write(command_ps.stdout.decode("utf-8"))

            if outVcf.endswith('.vcf.gz'):
                pysam.tabix_compress(tmpVcf, outVcf, force=True)
                pysam.tabix_index(outVcf, force=True, preset='vcf')
                os.unlink(tmpVcf)
        else:
            raise subprocess.CalledProcessError(cmd=command_ps.args, returncode=command_ps.returncode, output=command_ps.stdout)


def get_data_dir(config_file):
    data_dir = ""
    with open(config_file, 'rt') as inf:
        for line in inf:
            if line.strip().startswith('data.dir'):
                data_dir = line[line.find("=") + 1:].strip()
                break
    return data_dir


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
                if i >= 0:
                    genomes.add(line[:i])
    with open(config_file, 'at') as outf:
        for (g, d, c) in new_genomes:
            if g not in genomes:
                outf.write('{}.genome : {}\n'.format(g, d))
                if g != c:
                    outf.write("\t{}.chromosomes : {}\n".format(g, c))


class DownloadAndTweakSnpEff(tools.DownloadPackage):

    def __init__(self, url, extra_genomes=None):
        extra_genomes = extra_genomes or []

        self.extra_genomes = extra_genomes
        super(DownloadAndTweakSnpEff, self).__init__(url, 'snpEff/snpEff.jar', require_executability=False)

    def post_download(self):
        config_file = os.path.join(self.destination_dir, 'snpEff', 'snpEff.config')
        add_genomes_to_snpeff_config_file(config_file, zip(self.extra_genomes, self.extra_genomes, self.extra_genomes))

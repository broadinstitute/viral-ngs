'''
MEGAN - Metagenome Analyzer

MEGAN needs to be installed through the installer and a license also needs to
be obtained out of band.

The user needs to set environment variables:
MEGAN_PATH
MEGAN_LICENSE_PATH
MEGAN_DATA_PATH

The data directory can be populated with unzipped archives from
http://ab.inf.uni-tuebingen.de/data/software/megan5/download/welcome.html

If installed on MacOS, the likely MEGAN_PATH will be:
/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub
'''
from builtins import super
import os
import subprocess
import tempfile
import textwrap
import tools

SPECIES_PROJECTION = textwrap.dedent(
    '''
load taxGIFile='{gi_taxid}';
import blastFile='{blast_file}' meganFile='{megan_file}' maxMatches=100 minScore=50.0 maxExpected=0.01 topPercent=10.0 minSupportPercent=0.1 minSupport=1 minComplexity=0.3 useMinimalCoverageHeuristic=false useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=false textStoragePolicy=Embed blastFormat=BlastTAB mapping='Taxonomy:BUILT_IN=true,Taxonomy:GI_MAP=true';

collapse rank=Species;
compute profile=Projection rank=Species minPercent=1.0;
select rank=Species;
export what=CSV format=taxonname_count separator=tab counts=summarized file='{species_file}';
quit;
'''
)


@tools.skip_install_test()
class Megan(tools.Tool):
    """ Tool wrapper for the MEGAN metagenome analyzer
    """
    JVM_DEFAULT_MEM = '2G'

    def __init__(self, path=None):
        self.tool_version = None
        for megan_path in [path, os.environ.get('MEGAN_PATH')]:
            if not megan_path:
                continue
            if os.path.isdir(megan_path):
                megan_path = os.path.join(megan_path, 'MEGAN')
            install_methods = [
                tools.PrexistingUnixCommand(
                    megan_path,
                    verifycmd='{} -h > /dev/null'.format(megan_path),
                    verifycode=0,
                    require_executability=True
                )
            ]

        self.license_file = os.environ.get('MEGAN_LICENSE_PATH')
        self.data_dir = os.environ.get('MEGAN_DATA_PATH')

        super(Megan, self).__init__(install_methods=install_methods)

    @property
    def gi_taxid(self):
        return os.path.join(self.data_dir, 'gi_taxid-March2015X.bin')

    def execute(self, commands, memory=None, shell=False):
        with tempfile.NamedTemporaryFile(mode='w', prefix='megan_commands_', suffix='.txt') as command_file:
            command_file.write(commands)
            command_file.flush()

            memory = memory or Megan.JVM_DEFAULT_MEM
            megan = self.install_and_get_path()
            env = os.environ.copy()
            # Changing mem dynamically like this actually requires a slightly
            # modified version of MEGAN to work.
            env['INSTALL4J_ADD_VM_PARAMS'] = '-Xmx{}'.format(memory)
            megan_cmd = [
                megan, '--commandLineMode', '--licenseFile', self.license_file, '--commandFile', command_file.name
            ]
            if os.uname()[0] == 'Darwin':
                # OS X is an Aqua app which ignores $DISPLAY so we're going to
                # have to open the GUI.
                stdout = subprocess.check_output(megan_cmd, env=env)
                return 0, stdout
            elif 'linux' in os.uname()[0]:
                cmd = ['xvfb-run', '--auto-display']
                cmd.extend(megan_cmd)
                stdout = subprocess.check_output(cmd, env=env, shell=shell)
                return 0, stdout

    def species_projection(self, blast_file, output_file, memory=None):
        '''Calculate a species projection.

        Calculate a species projection of a blast m8 tab file with recommended
        filtering and support parameters from the manual. This requires the
        gi_taxid file downloaded to assign proper names.
        '''
        # We don't really care about the megan file since it's an on disk cache
        # of the imported blast file to make future imports faster. We just
        # want to do one shot analysis
        with tempfile.NamedTemporaryFile(prefix='megan_', suffix='.rma') as megan_file:
            commands = SPECIES_PROJECTION.format(
                gi_taxid=self.gi_taxid,
                blast_file=blast_file,
                megan_file=megan_file.name,
                species_file=output_file
            )
            self.execute(commands, memory=memory)

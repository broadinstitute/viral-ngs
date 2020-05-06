'''
    MAFFT - Multiple alignment program for amino acid or nucleotide sequences
    http://mafft.cbrc.jp/alignment/software/
'''

__author__ = "tomkinsc@broadinstitute.org"

from Bio import SeqIO
import logging
import tools
import util.file
import util.misc
import os
import os.path
import subprocess

TOOL_NAME = "mafft"
TOOL_VERSION = '7.464'

_log = logging.getLogger(__name__)


class MafftTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def __seqIdsAreAllUnique(self, filePath, inputFormat="fasta"):
        seqIds = []
        with open(filePath, "r") as inFile:
            fastaFile = SeqIO.parse(inFile, inputFormat)
            seqIds = [x.id for x in fastaFile]

        # collapse like IDs using set()
        if len(seqIds) > len(set(seqIds)):
            raise LookupError(
                "Not all sequence IDs in input are unique for file: {}".format(
                    os.path.basename(filePath)
                )
            )

    # pylint: disable=W0221
    def execute(
        self,
        inFastas,
        outFile,
        localpair,
        globalpair,
        preservecase,
        reorder,
        outputAsClustal,
        maxiters,
        gapOpeningPenalty=None,
        offset=None,
        threads=None,
        verbose=True,
        retree=None
    ):
        inputFileName = ""
        tempCombinedInputFile = ""

        # get the full paths of input and output files in case the user has specified relative paths
        inputFiles = []
        for f in inFastas:
            inputFiles.append(os.path.abspath(f))
        outFile = os.path.abspath(outFile)

        # ensure that all sequence IDs in each input file are unique
        # (otherwise the alignment result makes little sense)
        # we can check before combining to localize duplications to a specific file
        for filePath in inputFiles:
            self.__seqIdsAreAllUnique(filePath)

        # if multiple fasta files are specified for input
        if len(inputFiles) > 1:
            # combined specified input files into a single temp FASTA file so MAFFT can read them
            tempFileSuffix = ""
            for filePath in inputFiles:
                tempFileSuffix += "__" + os.path.basename(filePath)
            tempCombinedInputFile = util.file.mkstempfname('__combined.{}'.format(tempFileSuffix))
            with open(tempCombinedInputFile, "w") as outfile:
                for f in inputFiles:
                    with open(f, "r") as infile:
                        outfile.write(infile.read())
            inputFileName = tempCombinedInputFile
        # if there is only one file specified, just use it
        else:
            inputFileName = inputFiles[0]

        # check that all sequence IDs in a file are unique
        self.__seqIdsAreAllUnique(inputFileName)

        # change the pwd, since the shell script that comes with mafft depends on the pwd
        # being correct
        pwdBeforeMafft = os.getcwd()
        os.chdir(os.path.dirname(self.install_and_get_path()))

        # build the MAFFT command
        tool_cmd = [self.install_and_get_path()]

        if not (retree or localpair or globalpair):
            tool_cmd.append("--auto")
        tool_cmd.extend(["--thread", "{}".format(util.misc.sanitize_thread_count(threads))])

        if localpair and globalpair:
            raise Exception("Alignment type must be either local or global, not both.")

        if localpair:
            tool_cmd.append("--localpair")
            if not maxiters:
                maxiters = 1000
        if globalpair:
            tool_cmd.append("--globalpair")
            if not maxiters:
                maxiters = 1000
        if preservecase:
            tool_cmd.append("--preservecase")
        if reorder:
            tool_cmd.append("--reorder")
        if retree:
            tool_cmd.extend(["--retree", "{}".format(retree)])
        if gapOpeningPenalty:
            tool_cmd.extend(["--op", "{penalty}".format(penalty=gapOpeningPenalty)])
        if offset:
            tool_cmd.extend(["--ep", "{num}".format(num=offset)])
        if not verbose:
            tool_cmd.append("--quiet")
        if outputAsClustal:
            tool_cmd.append("--clustalout")
        if maxiters:
            tool_cmd.extend(["--maxiterate", "{iters}".format(iters=maxiters)])

        tool_cmd.append(inputFileName)

        _log.debug(' '.join(tool_cmd))

        # run the MAFFT alignment
        with open(outFile, 'w') as outf:
            util.misc.run_and_save(tool_cmd, outf=outf)

        if len(tempCombinedInputFile):
            # remove temp FASTA file
            os.unlink(tempCombinedInputFile)

        # restore pwd
        os.chdir(pwdBeforeMafft)

        return outFile
    # pylint: enable=W0221

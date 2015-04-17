'''
    MAFFT - Multiple alignment program for amino acid or nucleotide sequences
    http://mafft.cbrc.jp/alignment/software/
'''

__author__ = "tomkinsc@broadinstitute.org"


from Bio import SeqIO
import logging, tools, util.file
import os, os.path, subprocess

tool_version = '7.220'
url = 'http://mafft.cbrc.jp/alignment/software/mafft-{ver}-{os}.{ext}'

log = logging.getLogger(__name__)

class MafftTool(tools.Tool):
    def __init__(self, install_methods = None):
        if install_methods == None:
            install_methods = []
            mafft_os                = get_mafft_os()
            mafft_bitdepth          = get_mafft_bitdepth()
            mafft_archive_extension = get_mafft_archive_extension(mafft_os)
            binaryPath              = get_mafft_binary_path(mafft_os, mafft_bitdepth)
            binaryDir               = get_mafft_binary_path(mafft_os, mafft_bitdepth, full=False)
            
            target_rel_path = '{binPath}'.format(binPath=binaryPath)
            verify_command  = 'cd {dir}/mafft-{ver}/{binDir} && {dir}/mafft-{ver}/{binPath} --version > /dev/null 2>&1'.format(dir=util.file.get_build_path(), ver=tool_version, binPath=binaryPath, binDir=binaryDir) 
            destination_dir = '{dir}/mafft-{ver}'.format(dir=util.file.get_build_path(), ver=tool_version)

            install_methods.append(
                    tools.DownloadPackage( url.format(ver=tool_version, os=mafft_os, ext=mafft_archive_extension),
                    target_rel_path=target_rel_path,
                    destination_dir=destination_dir,
                    verifycmd=verify_command))

        tools.Tool.__init__(self, install_methods = install_methods)

    def version(self):
        return tool_version

    def __seqIdsAreAllUnique(self, filePath, inputFormat="fasta"):
        fastaFile = SeqIO.parse(open(filePath, "rU"), inputFormat)
        seqIds = [x.id for x in fastaFile]

        # collapse like IDs using set()
        if len(seqIds) > len(set(seqIds)):
            raise LookupError("Not all sequence IDs in input are unique for file: {}".format(os.path.basename(filePath)))

    def execute(self, inFastas, outFile, localpair, globalpair, preservecase, reorder, 
                outputAsClustal, maxiters, gapOpeningPenalty=None, offset=None, threads=-1, verbose=True):

        inputFileName         = ""
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
        if len(inputFiles)>1:
            # combined specified input files into a single temp FASTA file so MAFFT can read them
            tempFileSuffix = ""
            for filePath in inputFiles:
                tempFileSuffix += "__" + os.path.basename(filePath)
            tempCombinedInputFile = util.file.mkstempfname('__combined.{}'.format(tempFileSuffix))
            with open(tempCombinedInputFile, "w") as outfile:
                for f in inputFiles:
                    with open(f, "r") as infile:
                        outfile.write(infile.read())
                #outFile.close()
            inputFileName = tempCombinedInputFile
        # if there is only once file specified, just use it
        else:
            inputFileName = inputFiles[0]

        # check that all sequence IDs in a file are unique
        self.__seqIdsAreAllUnique(inputFileName)

        # change the pwd, since the shell script that comes with mafft depends on the pwd
        # being correct
        pwdBeforeMafft = os.getcwd()
        os.chdir(os.path.dirname(self.install_and_get_path()))

        # build the MAFFT command
        toolCmd = [self.install_and_get_path()]

        toolCmd.append("--auto")
        if threads >= 1 or threads == -1:
            toolCmd.append("--thread {}".format(threads))
        else:
            raise Exception("invalid value for threads: {}".format(threads))

        if localpair and globalpair:
            raise Exception("Alignment type must be either local or global, not both.")

        if localpair:
            toolCmd.append("--localpair")
            if not maxiters:
                maxiters = 1000
        if globalpair:
            toolCmd.append("--globalpair")
            if not maxiters:
                maxiters = 1000
        if preservecase:
            toolCmd.append("--preservecase")
        if reorder:
            toolCmd.append("--reorder")
        if gapOpeningPenalty:
            toolCmd.append("--op {penalty}".format(penalty=gapOpeningPenalty))
        if offset:
            toolCmd.append("--ep {num}".format(num=offset))
        if not verbose:
            toolCmd.append("--quiet")
        if outputAsClustal:
            toolCmd.append("--clustalout")
        if maxiters:
            toolCmd.append("--maxiterate {iters}".format(iters=maxiters))
        
        toolCmd.append(inputFileName)

        log.debug(' '.join(toolCmd))

        # run the MAFFT alignment
        with open(outFile, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)

        if len(tempCombinedInputFile):
            # remove temp FASTA file
            os.unlink(tempCombinedInputFile)

        # restore pwd
        os.chdir(pwdBeforeMafft)

def get_mafft_os():
    uname = os.uname()
    if uname[0] == "Darwin":
        return "mac"
    if uname[0] == "Linux":
        return "linux"

def get_mafft_archive_extension(mafft_os):
    if mafft_os == "mac":
        return "zip"
    elif mafft_os == "linux":
        return "tgz"

def get_mafft_bitdepth():
    uname = os.uname()
    if uname[4] == "x86_64":
        return "64"
    if uname[4] in ['i386','i686',"x86"]:
        return "32"

def get_mafft_binary_path(mafft_os, bitdepth, full=True):
    mafftPath = ""

    if mafft_os == "mac":
        mafftPath += "mafft-mac/"
    elif mafft_os == "linux":
        mafftPath += "mafft-linux{bit}".format(bit=bitdepth) + "/"

    if full:
        mafftPath += "mafft.bat"

    return mafftPath





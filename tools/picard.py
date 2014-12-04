"Tools in the 'Picard' suite."
import logging, os, os.path, subprocess, tempfile
import tools, util.file

tool_version = '1.126'
url = 'https://github.com/broadinstitute/picard/releases/download/' \
        + '{ver}/picard-tools-{ver}.zip'.format(ver=tool_version)
# Note: Version 1.126 is latest as of 2014-12-02
# Note: /seq/software/picard/{versionnumber}/ does not correspond with github release numbers!

log = logging.getLogger(__name__)

class PicardTools(tools.Tool) :
    """Base class for tools in the picard suite."""
    jvmMemDefault = '2g'
    def __init__(self, install_methods = None) :
        if install_methods == None :
            target_rel_path = 'picard-tools-{}/picard.jar'.format(tool_version)
            install_methods = [
                tools.DownloadPackage(url, target_rel_path, require_executability=False)]
        tools.Tool.__init__(self, install_methods = install_methods)
    def version(self) :
        return tool_version
    def execute(self, command, picardOptions=[], JVMmemory=None) :
        if JVMmemory==None:
            JVMmemory = self.jvmMemDefault
        toolCmd = ['java',
            '-Xmx' + JVMmemory,
            '-Djava.io.tmpdir=' + tempfile.tempdir,
            '-jar', self.install_and_get_path(),
            command] + picardOptions
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)
    def dict_to_picard_opts(self, options) :
        return ["%s=%s" % (k,v) for k,v in options.items()]

class RevertSamTool(PicardTools) :
    subtoolName = 'RevertSam'
    def execute(self, inBam, outBam,
                picardOptions=[], JVMmemory=None) :
        opts = ['INPUT='+inBam, 'OUTPUT='+outBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class MarkDuplicatesTool(PicardTools) :
    subtoolName = 'MarkDuplicates'
    def execute(self, inBams, outBam, outMetrics=None,
                picardOptions=[], JVMmemory=None) :
        if not outMetrics :
            outMetrics = util.file.mkstempfname('.metrics')
        opts = ['INPUT='+bam for bam in inBams] + [
            'OUTPUT='+outBam, 'METRICS_FILE='+outMetrics]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class SamToFastqTool(PicardTools) :
    subtoolName = 'SamToFastq'
    def execute(self, inBam, outFastq1, outFastq2,
                picardOptions=[], JVMmemory=None) :
        opts = ['INPUT='+inBam,
            'FASTQ='+outFastq1, 'SECOND_END_FASTQ='+outFastq2,
            'VALIDATION_STRINGENCY=SILENT']
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class FastqToSamTool(PicardTools) :
    subtoolName = 'FastqToSam'
    def execute(self, inFastq1, inFastq2, sampleName, outBam,
                picardOptions=[], JVMmemory=None) :
        opts = ['FASTQ='+inFastq1, 'FASTQ2='+inFastq2,
             'OUTPUT='+outBam, 'SAMPLE_NAME='+sampleName]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class SortSamTool(PicardTools) :
    subtoolName = 'SortSam'
    valid_sort_orders = ['unsorted', 'queryname', 'coordinate']
    default_sort_order = 'coordinate'
    def execute(self, inBam, outBam, sort_order = default_sort_order,
                picardOptions=[], JVMmemory=None) :
        if sort_order not in self.valid_sort_orders :
            raise Exception("invalid sort order")
        opts = ['INPUT='+inBam, 'OUTPUT='+outBam, 'SORT_ORDER='+sort_order]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class CreateSequenceDictionaryTool(PicardTools) :
    subtoolName = 'CreateSequenceDictionary'
    def execute(self, inFasta, outDict=None, overwrite=False,
                picardOptions=[], JVMmemory=None) :
        if not outDict:
            if inFasta.lower().endswith('.fa'):
                outDict = inFasta[:-3]
            elif inFasta.lower().endswith('.fasta'):
                outDict = inFasta[:-6]
            else:
                raise Exception("bad input")
        if os.path.isfile(outDict):
            if overwrite:
                os.unlink(outDict)
            else:
                return
        opts = ['REFERENCE='+inFasta, 'OUTPUT='+outDict]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class ExtractIlluminaBarcodesTool(PicardTools) :
    subtoolName = 'ExtractIlluminaBarcodes'
    def execute(self, basecalls_dir, lane, read_structure, barcode_file, metrics,
                picardOptions=[], JVMmemory=None) :
        opts = ['BASECALLS_DIR='+basecalls_dir, 'LANE='+str(lane),
            'READ_STRUCTURE='+read_structure,
            'BARCODE_FILE='+barcode_file,
            'METRICS_FILE='+metrics]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class IlluminaBasecallsToSamTool(PicardTools) :
    subtoolName = 'IlluminaBasecallsToSam'
    jvmMemDefault = '54g'
    def execute(self, basecalls_dir, lane, run_barcode, read_structure, library_params,
                picardOptions=[], JVMmemory=None) :
        opts = ['BASECALLS_DIR='+basecalls_dir, 'LANE='+str(lane),
            'READ_STRUCTURE='+read_structure,
            'RUN_BARCODE='+run_barcode, 'LIBRARY_PARAMS='+library_params]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)
        

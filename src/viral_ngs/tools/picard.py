'''
    Tools in the Picard suite
'''

import logging
import os
import os.path
import tempfile
import shutil
import pysam
import tools
import util.file
import util.misc

TOOL_NAME = "picard"
TOOL_VERSION = '1.126'
TOOL_URL = 'https://github.com/broadinstitute/picard/releases/download/' \
    + '{ver}/picard-tools-{ver}.zip'.format(ver=TOOL_VERSION)
# Note: Version 1.126 is latest as of 2014-12-02
# Note: /seq/software/picard/{versionnumber}/ does not correspond with github release numbers!

_log = logging.getLogger(__name__)


class PicardTools(tools.Tool):
    """Base class for tools in the picard suite."""
    jvmMemDefault = '2g'

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else None

        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, executable=self.subtool_name, version=TOOL_VERSION))

            target_rel_path = 'picard-tools-{}/picard.jar'.format(TOOL_VERSION)
            install_methods.append(tools.DownloadPackage(TOOL_URL, target_rel_path, require_executability=False))
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if JVMmemory is None:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        if self.install_and_get_path().endswith(".jar"):
            tool_cmd = [
                'java', '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.tempdir, '-jar', self.install_and_get_path(),
                command
            ] + picardOptions
        else:
            tool_cmd = [
                self.install_and_get_path(), '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.tempdir, command
            ] + picardOptions
        _log.debug(' '.join(tool_cmd))
        util.misc.run_and_print(tool_cmd, check=True)

    @staticmethod
    def dict_to_picard_opts(options):
        return ["%s=%s" % (k, v) for k, v in options.items()]


class RevertSamTool(PicardTools):
    subtoolName = 'RevertSam'

    def execute(self, inBam, outBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []
        opts = ['INPUT=' + inBam, 'OUTPUT=' + outBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class MarkDuplicatesTool(PicardTools):
    subtoolName = 'MarkDuplicates'

    def execute(self, inBams, outBam, outMetrics=None, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if not outMetrics:
            outMetrics = util.file.mkstempfname('.metrics')
        opts = ['INPUT=' + bam for bam in inBams] + ['OUTPUT=' + outBam, 'METRICS_FILE=' + outMetrics]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class SamToFastqTool(PicardTools):
    subtoolName = 'SamToFastq'
    illumina_clipping_attribute = 'XT'

    def execute(self, inBam, outFastq1, outFastq2, picardOptions=None, JVMmemory=None):  # pylint: disable=W0221
        picardOptions = picardOptions or []
        opts = ['FASTQ=' + outFastq1, 'SECOND_END_FASTQ=' + outFastq2,
                'INPUT=' + inBam, 'VALIDATION_STRINGENCY=SILENT']
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

    def per_read_group(self, inBam, outDir, picardOptions=None, JVMmemory=None):
        picardOptions = picardOptions or []

        opts = ['INPUT=' + inBam, 'OUTPUT_DIR=' + outDir, 'OUTPUT_PER_RG=true']
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class FastqToSamTool(PicardTools):
    subtoolName = 'FastqToSam'

    def execute(self,
                inFastq1,
                inFastq2,
                sampleName,
                outBam,
                picardOptions=None,
                JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if inFastq2:
            opts = ['FASTQ=' + inFastq1, 'FASTQ2=' + inFastq2, 'OUTPUT=' + outBam, 'SAMPLE_NAME=' + sampleName]
        else:
            opts = ['FASTQ=' + inFastq1, 'OUTPUT=' + outBam, 'SAMPLE_NAME=' + sampleName]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class SortSamTool(PicardTools):
    subtoolName = 'SortSam'
    valid_sort_orders = ['unsorted', 'queryname', 'coordinate']
    default_sort_order = 'coordinate'

    def execute(self,
                inBam,
                outBam,
                sort_order=default_sort_order,
                picardOptions=None,
                JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if sort_order not in self.valid_sort_orders:
            raise Exception("invalid sort order")
        opts = ['INPUT=' + inBam, 'OUTPUT=' + outBam, 'SORT_ORDER=' + sort_order]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class MergeSamFilesTool(PicardTools):
    subtoolName = 'MergeSamFiles'

    def execute(self, inBams, outBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        opts = ['INPUT=' + bam for bam in inBams] + ['OUTPUT=' + outBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class FilterSamReadsTool(PicardTools):
    ''' TO DO: it might be desirable to replace this tool with a
        non-Picard/non-Java approach that uses samtools/pysam, sqlite,
        and O(1) memory.
    '''
    subtoolName = 'FilterSamReads'
    jvmMemDefault = '4g'

    def execute(self, inBam, exclude, readList, outBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if os.path.getsize(readList) == 0:
            # Picard FilterSamReads cannot deal with an empty READ_LIST_FILE
            if exclude:
                shutil.copyfile(inBam, outBam)
            else:
                tmpf = util.file.mkstempfname('.sam')
                if inBam.endswith('.sam'):
                    # output format (sam/bam) is inferred by samtools based on file extension
                    header = pysam.view('-o', tmpf, '-H', '-S', inBam, catch_stdout=False)
                else:
                    header = pysam.view('-o', tmpf, '-H', inBam, catch_stdout=False)
                # pysam.AlignmentFile cannot write an empty file
                # samtools cannot convert SAM -> BAM on an empty file
                # but Picard SamFormatConverter can deal with empty files
                opts = ['INPUT=' + tmpf, 'OUTPUT=' + outBam, 'VERBOSITY=ERROR']
                PicardTools.execute(self, 'SamFormatConverter', opts, JVMmemory='50m')
        else:
            opts = [
                'INPUT=' + inBam, 'OUTPUT=' + outBam, 'READ_LIST_FILE=' + readList,
                'FILTER=' + (exclude and 'excludeReadList' or 'includeReadList'), 'WRITE_READS_FILES=false'
            ]
            PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class CreateSequenceDictionaryTool(PicardTools):
    subtoolName = 'CreateSequenceDictionary'
    jvmMemDefault = '512m'

    def execute(self,
                inFasta,
                outDict=None,
                overwrite=False,
                picardOptions=None,
                JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if not outDict:
            if inFasta.lower().endswith('.fa'):
                outDict = inFasta[:-3] + '.dict'
            elif inFasta.lower().endswith('.fasta'):
                outDict = inFasta[:-6] + '.dict'
            else:
                raise Exception("bad input")
        if os.path.isfile(outDict):
            if overwrite:
                os.unlink(outDict)
            else:
                return
        opts = ['REFERENCE=' + inFasta, 'OUTPUT=' + outDict]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class BuildBamIndexTool(PicardTools):
    subtoolName = 'BuildBamIndex'
    jvmMemDefault = '512m'

    def execute(self, inBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        opts = ['INPUT=' + inBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class ExtractIlluminaBarcodesTool(PicardTools):
    subtoolName = 'ExtractIlluminaBarcodes'
    jvmMemDefault = '8g'
    defaults = {'read_structure': '101T8B8B101T', 'max_mismatches': 0, 'minimum_base_quality': 25, 'num_processors': 0}
    option_list = (
        'read_structure', 'max_mismatches', 'minimum_base_quality', 'min_mismatch_delta', 'max_no_calls',
        'minimum_quality', 'compress_outputs', 'num_processors'
    )

    def execute(self,
                basecalls_dir,
                lane,
                barcode_file,
                output_dir,
                metrics,
                picardOptions=None,
                JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or {}

        opts_dict = self.defaults.copy()
        for k, v in picardOptions.items():
            opts_dict[k] = v
        opts = []
        for k, v in opts_dict.items():
            if v is not None:
                if type(v) in (list, tuple):
                    for x in v:
                        opts.append('='.join((k.upper(), str(x))))
                else:
                    opts.append('='.join((k.upper(), str(v))))
        opts += [
            'BASECALLS_DIR=' + basecalls_dir, 'LANE=' + str(lane), 'BARCODE_FILE=' + barcode_file,
            'METRICS_FILE=' + metrics
        ]
        if output_dir is not None:
            opts += ['OUTPUT_DIR=' + output_dir]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory)


class IlluminaBasecallsToSamTool(PicardTools):
    subtoolName = 'IlluminaBasecallsToSam'
    jvmMemDefault = '54g'
    defaults = {
        'read_structure': '101T8B8B101T',
        'adapters_to_check': ('PAIRED_END', 'NEXTERA_V1', 'NEXTERA_V2'),
        'max_reads_in_ram_per_tile': 100000,
        'max_records_in_ram': 100000,
        'num_processors': 4,
        'force_gc': False,
        'include_non_pf_reads': False,
    }
    option_list = (
        'read_structure', 'sequencing_center', 'adapters_to_check', 'platform', 'max_reads_in_ram_per_tile',
        'max_records_in_ram', 'num_processors', 'apply_eamss_filter', 'force_gc', 'first_tile', 'tile_limit',
        'include_non_pf_reads', 'run_start_date', 'read_group_id'
    )

    # pylint: disable=W0221
    def execute(
        self,
        basecalls_dir,
        barcodes_dir,
        run_barcode,
        lane,
        library_params,
        picardOptions=None,
        JVMmemory=None
    ):
        picardOptions = picardOptions or {}

        opts_dict = self.defaults.copy()
        for k, v in picardOptions.items():
            opts_dict[k] = v
        opts = []
        for k, v in opts_dict.items():
            if v is not None:
                if type(v) in (list, tuple):
                    for x in v:
                        opts.append('='.join((k.upper(), str(x))))
                else:
                    opts.append('='.join((k.upper(), str(v))))
        opts += [
            'BASECALLS_DIR=' + basecalls_dir, 'BARCODES_DIR=' + barcodes_dir, 'LANE=' + str(lane),
            'RUN_BARCODE=' + run_barcode, 'LIBRARY_PARAMS=' + library_params
        ]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory)
    # pylint: enable=W0221

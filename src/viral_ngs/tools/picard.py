'''
    Tools in the Picard suite
'''

import logging
import os
import os.path
import tempfile
import shutil
import subprocess
import contextlib
from decimal import *

import pysam

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = "picard"
TOOL_VERSION = '2.21.1'
TOOL_URL = 'https://github.com/broadinstitute/picard/releases/download/' \
    + '{ver}/picard-tools-{ver}.zip'.format(ver=TOOL_VERSION)
# Note: /seq/software/picard/{versionnumber}/ does not correspond with github release numbers!

_log = logging.getLogger(__name__)


class PicardTools(tools.Tool):
    """Base class for tools in the picard suite."""
    jvmMemDefault = '2g'

    def install(self):
        pass

    def is_installed(self):
        return True

    def install_and_get_path(self):
        # the conda version wraps the jar file with a shell script
        return 'picard'

    def __init__(self, install_methods=None):
        self.subtool_name = self.subtool_name if hasattr(self, "subtool_name") else None
        self.installed_method = True

    def version(self):
        return TOOL_VERSION

    def execute(self, command, picardOptions=None, JVMmemory=None, background=False, **kwargs):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if JVMmemory is None:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        path = self.install_and_get_path()
        tool_cmd = [path, '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), command] + picardOptions + ['USE_JDK_DEFLATER=true','USE_JDK_INFLATER=true']
        _log.debug(' '.join(tool_cmd))

        env = os.environ.copy()
        env.pop('JAVA_HOME', None)
        if background:
            return subprocess.Popen(tool_cmd, env=env, **kwargs)
        else:
            return subprocess.check_call(tool_cmd, env=env, **kwargs)

    @staticmethod
    def dict_to_picard_opts(options):
        return ["%s=%s" % (k, v) for k, v in options.items()]


class RevertSamTool(PicardTools):
    subtoolName = 'RevertSam'

    def execute(self, inBam, outBam, picardOptions=None, JVMmemory=None, background=False):    # pylint: disable=W0221
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            shutil.copyfile(inBam, outBam)
        else:
            picardOptions = picardOptions or []
            opts = ['INPUT=' + inBam, 'OUTPUT=' + outBam]
            PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory, background=background)


class CheckIlluminaDirectoryTool(PicardTools):
    subtoolName = 'CheckIlluminaDirectory'

    def execute(self, basecalls_dir, lanes,  read_structure, data_types=None, fake_files=False, tile_numbers=None, link_locs=False, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []
        opts = [
            'BASECALLS_DIR=' + basecalls_dir,
            'READ_STRUCTURE=' + read_structure
        ]

        if fake_files:
            opts += ['FAKE_FILES=true']

        if tile_numbers is not None:
            if type(tile_numbers)==int:
                tile_numbers = [tile_numbers]
            for tile_number in set(tile_numbers):
                 opts += ['TILE_NUMBERS=' + str(tile_number)]

        if data_types is not None:
            if isinstance(arg, str):
                data_types = [data_types]
            for data_type in set(data_types):
                opts += ['DATA_TYPES=' + data_type]

        # if lanes is a single int, cast it to a list
        if type(lanes)==int:
            lanes = [lanes]

        assert type(lanes)==list, "Lanes must be a list specifying the lanes"
        for lane in set(lanes):
             opts += ['LANES=' + str(lane)]

        if link_locs:
            opts += ['LINK_LOCS=true']

        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class MarkDuplicatesTool(PicardTools):
    subtoolName = 'MarkDuplicates'

    def execute(
        self, inBams, outBam, outMetrics=None, taggingPolicy=None, picardOptions=None, JVMmemory=None
    ):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if not outMetrics:
            outMetrics = util.file.mkstempfname('.metrics')
        opts = ['INPUT=' + bam for bam in inBams] + ['OUTPUT=' + outBam, 'METRICS_FILE=' + outMetrics]
        if taggingPolicy:
            opts += ['TAGGING_POLICY={}'.format(taggingPolicy)]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class SplitSamByLibraryTool(PicardTools):
    subtoolName = 'SplitSamByLibrary'

    def execute(
        self, in_bam, out_dir, picardOptions=None, JVMmemory=None
    ):    # pylint: disable=W0221

        if tools.samtools.SamtoolsTool().isEmpty(in_bam):
            # Picard SplitSamByLibrary cannot deal with an empty input BAM file
            shutil.copyfile(in_bam, os.path.join(out_dir, 'output.bam'))
            return

        opts = ['INPUT=' + in_bam, 'OUTPUT=' + out_dir]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory=JVMmemory)


class SamToFastqTool(PicardTools):
    subtoolName = 'SamToFastq'
    illumina_clipping_attribute = 'XT'

    def execute(self, inBam, outFastq1, outFastq2=None, outFastq0=None,
                illuminaClipping=False, interleave=False,
                picardOptions=None, JVMmemory=None, background=None, **kwargs):    # pylint: disable=W0221
        '''Write paired reads from `inBam` to `outFastq1` and `outFastq1`.  If `outFastq0` is given, write
        any unpaired reads from `inBam` there, else ignore them.  If `illuminaClipping` is True,
        trim reads at the clipping position specified by the Illumina clipping attribute
        (which is defined by the class variable SamToFastqTool.illumina_clipping_attribute).'''

        picardOptions = picardOptions or []

        opts = [
            'INPUT=' + inBam, 'VALIDATION_STRINGENCY=SILENT', 'FASTQ=' + outFastq1,
        ]
        if outFastq2:
            opts.append('SECOND_END_FASTQ=' + outFastq2)
        else:
            if interleave:
                opts.append('INTERLEAVE=true')

        if outFastq0:
            assert outFastq2, "outFastq0 option only applies in paired-end output mode"
            opts.append('UNPAIRED_FASTQ=' + outFastq0)

        if illuminaClipping:
            opts += PicardTools.dict_to_picard_opts({
                'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
                'CLIPPING_ACTION': 'X'
            })

        return PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory,
                                   background=background, **kwargs)

    @contextlib.contextmanager
    def execute_tmp(self, inBam, sfx='', includeUnpaired=False, **kwargs):
        '''Output reads from `inBam` to temp fastq files, and yield their filenames as a tuple.
        If `includeUnpaired` is True, output unpaired reads to a third temp fastq file and yield it as a third
        element of the tuple.
        '''
        if not includeUnpaired:
            with util.file.tempfnames((sfx+'.1.fq', sfx+'.2.fq')) as (outFastq1, outFastq2):
                self.execute(inBam, outFastq1, outFastq2, **kwargs)
                yield outFastq1, outFastq2
        else:
            with util.file.tempfnames((sfx+'.1.fq', sfx+'.2.fq', sfx+'.0.fq')) as (outFastq1, outFastq2, outFastq0):
                self.execute(inBam, outFastq1, outFastq2, outFastq0=outFastq0, **kwargs)
                yield outFastq1, outFastq2, outFastq0

    def per_read_group(self, inBam, outDir, picardOptions=None, JVMmemory=None):
        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # Picard SamToFastq cannot deal with an empty input BAM file
            if not os.path.isdir(outDir):
                os.mkdir(outDir)
        else:
            picardOptions = picardOptions or []
            opts = ['INPUT=' + inBam, 'OUTPUT_DIR=' + outDir, 'OUTPUT_PER_RG=true', 'RG_TAG=ID']
            PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class FastqToSamTool(PicardTools):
    subtoolName = 'FastqToSam'

    def execute(
        self, inFastq1, inFastq2, sampleName,
        outBam, picardOptions=None, JVMmemory=None
    ):    # pylint: disable=W0221
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

    def execute(
        self, inBam, outBam, sort_order=default_sort_order,
        picardOptions=None, JVMmemory=None
    ):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if sort_order not in self.valid_sort_orders:
            raise Exception("invalid sort order")
        opts = ['INPUT=' + inBam, 'OUTPUT=' + outBam, 'SORT_ORDER=' + sort_order]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class DownsampleSamTool(PicardTools):
    subtoolName = 'DownsampleSam'
    valid_strategies = ['HighAccuracy', 'ConstantMemory', 'Chained']
    default_strategy = 'Chained'    # ConstantMemory first, then HighAccuracy to get closer to the target probability
    default_random_seed = 1    # Set to constant int for deterministic results
    jvmMemDefault = '4g'

    def execute(self,
                inBam,
                outBam,
                probability,
                accuracy=None, #Picard default is 1.0E-4
                strategy=default_strategy,
                random_seed=default_random_seed,
                picardOptions=None,
                JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []
        JVMmemory = JVMmemory or self.jvmMemDefault

        if strategy not in self.valid_strategies:
            raise Exception("invalid subsample strategy: %s" % strategy)
        if not probability:
            raise Exception("Probability must be defined")
        if float(probability) <= 0 or float(probability) > 1:
            raise Exception("Probability must be in range (0,1]. This value was given: %s" % probability)

        opts = ['INPUT=' + inBam, 'OUTPUT=' + outBam, 'PROBABILITY=' + str(probability)]

        if accuracy:
            opts.extend(['ACCURACY=' + str(accuracy)])

        if strategy:
            opts.extend(['STRATEGY=' + strategy])

        if not random_seed:
            _log.info("No random seed is set for subsample operation; results will be non-deterministic")
            opts.extend(["RANDOM_SEED=null"])
            raise Exception(
                "Setting RANDOM_SEED=null crashes Picard 1.141, though it may be available when viral-ngs updates to a later version of Picard."
            )
        else:
            _log.info("Random seed is set for subsample operation; results will be deterministic")
            opts.extend(['RANDOM_SEED=' + str(random_seed)])

        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


    def downsample_to_approx_count(
        self, inBam, outBam, read_count, picardOptions=None,
        JVMmemory=None
    ):    # pylint: disable=W0221):

        samtools = tools.samtools.SamtoolsTool()
        total_read_count = samtools.count(inBam)

        if total_read_count == 0:
            _log.info("Input BAM has no reads. Copying to output.")
            shutil.copyfile(inBam, outBam)

        probability = Decimal(int(read_count)) / Decimal(total_read_count)
        probability = 1 if probability > 1 else probability

        assert probability >= 0

        if probability < 1:
            # per the Picard docs, HighAccuracy is recommended for read counts <50k
            strategy = "HighAccuracy" if total_read_count < 50000 else "Chained"
            _log.info("Setting downsample accuracy to %s based on read count of %s" % (strategy, total_read_count))

            self.execute(inBam, outBam, probability, strategy=strategy, accuracy=0.00001, picardOptions=picardOptions, JVMmemory=JVMmemory)
        else:
            _log.info("Requested downsample count exceeds number of reads. Including all reads in output.")
            shutil.copyfile(inBam, outBam)


class MergeSamFilesTool(PicardTools):
    subtoolName = 'MergeSamFiles'

    def execute(self, inBams, outBam, picardOptions=None, JVMmemory=None, background=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if not any(opt.startswith('USE_THREADING') for opt in picardOptions):
            picardOptions.append('USE_THREADING=true')

        opts = ['INPUT=' + bam for bam in inBams] + ['OUTPUT=' + outBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory=JVMmemory, background=background)


class ReplaceSamHeaderTool(PicardTools):
    subtoolName = 'ReplaceSamHeader'

    def execute(self, inBam, headerBam, outBam, picardOptions=None, JVMmemory=None, background=None):    # pylint: disable=W0221

        opts = ['INPUT=' + inBam,
                'HEADER=' + headerBam,
                'OUTPUT=' + outBam]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory=JVMmemory, background=background)


class FilterSamReadsTool(PicardTools):
    ''' TO DO: it might be desirable to replace this tool with a
        non-Picard/non-Java approach that uses samtools/pysam, sqlite,
        and O(1) memory.
    '''
    subtoolName = 'FilterSamReads'
    jvmMemDefault = '4g'

    def execute(self, inBam, exclude, readList, outBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        if tools.samtools.SamtoolsTool().isEmpty(inBam):
            # Picard FilterSamReads cannot deal with an empty input BAM file
            shutil.copyfile(inBam, outBam)
        elif os.path.getsize(readList) == 0:
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

    def execute(
        self, inFasta, outDict=None, overwrite=False,
        picardOptions=None, JVMmemory=None
    ):    # pylint: disable=W0221
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
        with util.file.fastas_with_sanitized_ids(inFasta, use_tmp=False) as sanitized_fastas:
            opts = ['REFERENCE=' + sanitized_fastas[0], 'OUTPUT=' + outDict]
            PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)


class BuildBamIndexTool(PicardTools):
    subtoolName = 'BuildBamIndex'
    jvmMemDefault = '512m'

    def execute(self, inBam, picardOptions=None, JVMmemory=None):    # pylint: disable=W0221
        picardOptions = picardOptions or []

        opts = ['INPUT=' + inBam]
        PicardTools.execute(self, self.subtoolName, opts + picardOptions, JVMmemory)

class CollectIlluminaLaneMetricsTool(PicardTools):
    subtoolName = 'CollectIlluminaLaneMetrics'
    jvmMemDefault = '8g'
    defaults = {'read_structure': '101T8B8B101T'}
    option_list = (
        'read_structure',
    )

    def execute(
        self, run_dir,
        output_dir, output_prefix,
        picardOptions=None,
        JVMmemory=None
    ):    # pylint: disable=W0221
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

        opts += ['RUN_DIRECTORY=' + run_dir]
        opts += ['OUTPUT_DIRECTORY=' + output_dir]
        opts += ['OUTPUT_PREFIX=' + output_prefix]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory)

class ExtractIlluminaBarcodesTool(PicardTools):
    subtoolName = 'ExtractIlluminaBarcodes'
    jvmMemDefault = '8g'
    # minimum_base_quality=20 used to accommodate NovaSeq, which with RTA3 writes only four Q-score values: 2, 12, 23, and 37
    defaults = {'read_structure': '101T8B8B101T', 'max_mismatches': 0, 'minimum_base_quality': 20, 'num_processors': 0}
    option_list = (
        'read_structure', 'max_mismatches', 'minimum_base_quality', 'min_mismatch_delta', 'max_no_calls',
        'minimum_quality', 'compress_outputs', 'num_processors'
    )

    def execute(
        self, basecalls_dir,
        lane, barcode_file,
        output_dir, metrics,
        picardOptions=None,
        JVMmemory=None
    ):    # pylint: disable=W0221
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
    jvmMemDefault = '7g'
    defaults = {
        'read_structure': '101T8B8B101T',
        'adapters_to_check': ('PAIRED_END', 'NEXTERA_V1', 'NEXTERA_V2'),
        'max_reads_in_ram_per_tile': 1000000,
        'max_records_in_ram': 2000000,
        'num_processors': 0,
        'include_non_pf_reads': False,
        'compression_level': 7,
    }
    option_list = (
        'read_structure', 'sequencing_center', 'adapters_to_check', 'platform', 'max_reads_in_ram_per_tile',
        'max_records_in_ram', 'num_processors', 'apply_eamss_filter', 'force_gc', 'first_tile', 'tile_limit',
        'include_non_pf_reads', 'run_start_date', 'read_group_id', 'compression_level'
    )

    # pylint: disable=W0221
    def execute(self, 
        basecalls_dir,
        barcodes_dir,
        read_group_id_precursor,
        lane, library_params,
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
            'RUN_BARCODE=' + read_group_id_precursor, 'LIBRARY_PARAMS=' + library_params
        ]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory)

    def execute_single_sample(self, 
        basecalls_dir,
        output_file,
        read_group_id_precursor, # this is a run-specific ID, ex. the flowcell ID
        lane,
        sample_alias,
        picardOptions=None,
        JVMmemory=None
    ):
        picardOptions = picardOptions or {}

        assert len(read_group_id_precursor) >= 5, "read_group_id_precursor must be >=5 chars per https://github.com/broadinstitute/picard/blob/2.17.6/src/main/java/picard/illumina/IlluminaBasecallsToSam.java#L510"

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
            'OUTPUT=' + output_file,
            'BASECALLS_DIR=' + basecalls_dir, 
            'LANE=' + str(lane),
            'RUN_BARCODE=' + read_group_id_precursor, #
            'SAMPLE_ALIAS=' + sample_alias,
            'LIBRARY_NAME=' + sample_alias
        ]
        PicardTools.execute(self, self.subtoolName, opts, JVMmemory)
    # pylint: enable=W0221

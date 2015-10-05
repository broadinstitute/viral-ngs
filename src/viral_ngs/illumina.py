#!/usr/bin/env python
"""
Utilities for demultiplexing Illumina data.
"""

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import os
import os.path
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree

import util.cmd
import util.file
import tools.picard

log = logging.getLogger(__name__)

# =========================
# ***  illumina_demux   ***
# =========================


def parser_illumina_demux(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Illumina BCL directory (or tar.gz of BCL directory).')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('outDir', help='Output directory for BAM files.')

    parser.add_argument('--outMetrics',
                        help='Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.',
                        default=None)
    parser.add_argument('--sampleSheet',
                        default=None,
                        help='''Override SampleSheet. Input tab or CSV file w/header and four named columns:
                                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2.
                                Default is to look for a SampleSheet.csv in the inDir.''')
    parser.add_argument('--flowcell', help='Override flowcell ID (default: read from RunInfo.xml).', default=None)
    parser.add_argument('--read_structure',
                        help='Override read structure (default: read from RunInfo.xml).',
                        default=None)

    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        if opt not in ('read_structure', 'num_processors'):
            parser.add_argument('--' + opt,
                                help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt))
    for opt in tools.picard.IlluminaBasecallsToSamTool.option_list:
        if opt == 'adapters_to_check':
            parser.add_argument('--' + opt,
                                nargs='*',
                                help='Picard IlluminaBasecallsToSam ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
        elif opt == 'read_structure':
            pass
        else:
            parser.add_argument('--' + opt,
                                help='Picard IlluminaBasecallsToSam ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))

    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.IlluminaBasecallsToSamTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_illumina_demux)
    return parser


def main_illumina_demux(args):
    ''' Demultiplex Illumina runs & produce BAM files, one per sample.
        Wraps together Picard's ExtractBarcodes and IlluminaBasecallsToSam
        while handling the various required input formats. Also can
        read Illumina BCL directories, tar.gz BCL directories.
        TO DO: read BCL or tar.gz BCL directories from S3 / object store.
    '''

    # prepare
    illumina = IlluminaDirectory(args.inDir)
    illumina.load()
    if args.flowcell:
        flowcell = args.flowcell
    else:
        flowcell = illumina.get_RunInfo().get_flowcell()
    if args.run_start_date:
        run_date = args.run_start_date
    else:
        run_date = illumina.get_RunInfo().get_rundate_american()
    if args.read_structure:
        read_structure = args.read_structure
    else:
        read_structure = illumina.get_RunInfo().get_read_structure()
    if args.sampleSheet:
        samples = SampleSheet(args.sampleSheet, only_lane=args.lane)
    else:
        samples = illumina.get_SampleSheet(only_lane=args.lane)

    # Picard ExtractIlluminaBarcodes
    extract_input = util.file.mkstempfname('.txt', prefix='.'.join(['barcodeData', flowcell, str(args.lane)]))
    barcodes_tmpdir = tempfile.mkdtemp(prefix='extracted_barcodes-')
    samples.make_barcodes_file(extract_input)
    out_metrics = (args.outMetrics is None) and util.file.mkstempfname('.metrics.txt') or args.outMetrics
    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['read_structure'] = read_structure
    tools.picard.ExtractIlluminaBarcodesTool().execute(
        illumina.get_BCLdir(),
        args.lane,
        extract_input,
        barcodes_tmpdir,
        out_metrics,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory)

    # Picard IlluminaBasecallsToSam
    basecalls_input = util.file.mkstempfname('.txt', prefix='.'.join(['library_params', flowcell, str(args.lane)]))
    samples.make_params_file(args.outDir, basecalls_input)
    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['run_start_date'] = run_date
    picardOpts['read_structure'] = read_structure
    if not picardOpts.get('sequencing_center') and illumina.get_RunInfo():
        picardOpts['sequencing_center'] = illumina.get_RunInfo().get_machine()
    tools.picard.IlluminaBasecallsToSamTool().execute(
        illumina.get_BCLdir(),
        barcodes_tmpdir,
        flowcell,
        args.lane,
        basecalls_input,
        picardOptions=picardOpts,
        JVMmemory=args.JVMmemory)

    # clean up
    os.unlink(extract_input)
    os.unlink(basecalls_input)
    shutil.rmtree(barcodes_tmpdir)
    illumina.close()
    return 0


__commands__.append(('illumina_demux', parser_illumina_demux))

# ============================
# ***  IlluminaDirectory   ***
# ============================


class IlluminaDirectory(object):
    ''' A class that handles Illumina data directories
    '''

    def __init__(self, uri):
        self.uri = uri
        self.path = None
        self.tempDir = None
        self.runinfo = None
        self.samplesheet = None

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return 0

    def load(self):
        if self.path is None:
            if '://' in self.uri:
                raise NotImplementedError('boto s3 download here uri -> tarball')
                # tarball = util.file.mkstempfname('.tar.gz')
                # # TODO: download here, uri -> tarball
                # self._extract_tarball(tarball)
                # os.unlink(tarball)
            else:
                if os.path.isdir(self.uri):
                    self.path = self.uri
                else:
                    self._extract_tarball(self.uri)
            self._fix_path()

    def _fix_path(self):
        assert self.path is not None
        if not os.path.isdir(os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')):
            # this is not the correct root-level directory
            # sometimes this points to one level up
            subdirs = list(os.path.join(self.path, x) for x in os.listdir(self.path)
                           if os.path.isdir(os.path.join(self.path, x)))
            if len(subdirs) != 1:
                raise Exception('cannot find Data/Intensities/BaseCalls/ inside %s (subdirectories: %s)' %
                                (self.uri, str(subdirs)))
            self.path = subdirs[0]
            if not os.path.isdir(os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')):
                raise Exception('cannot find Data/Intensities/BaseCalls/ inside %s (%s)' % (self.uri, self.path))

    def _extract_tarball(self, tarfile):
        if not os.path.isfile(tarfile):
            raise Exception('file does not exist: %s' % tarfile)
        if tarfile.endswith('.tar.gz') or tarfile.endswith('.tgz'):
            compression_option = 'z'
        elif tarfile.endswith('.tar.bz2'):
            compression_option = 'j'
        elif tarfile.endwith('.tar'):
            compression_option = ''
        else:
            raise Exception("unsupported file type: %s" % tarfile)
        self.tempDir = tempfile.mkdtemp(prefix='IlluminaDirectory-')
        untar_cmd = ['tar', '-C', self.tempDir, '-x{}pf'.format(compression_option), tarfile]
        log.debug(' '.join(untar_cmd))
        with open(os.devnull, 'w') as fnull:
            subprocess.check_call(untar_cmd, stderr=fnull)
        self.path = self.tempDir

    def close(self):
        if self.tempDir:
            shutil.rmtree(self.tempDir)
            self.tempDir = None

    def get_RunInfo(self):
        if self.runinfo is None and os.path.isfile(os.path.join(self.path, 'RunInfo.xml')):
            self.runinfo = RunInfo(os.path.join(self.path, 'RunInfo.xml'))
        return self.runinfo

    def get_SampleSheet(self, only_lane=None):
        if self.samplesheet is None and os.path.isfile(os.path.join(self.path, 'SampleSheet.csv')):
            self.samplesheet = SampleSheet(os.path.join(self.path, 'SampleSheet.csv'), only_lane=only_lane)
        return self.samplesheet

    def get_BCLdir(self):
        return os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')

# ==================
# ***  RunInfo   ***
# ==================


class RunInfo(object):
    ''' A class that reads the RunInfo.xml file emitted by Illumina
        MiSeq and HiSeq machines.
    '''

    def __init__(self, xml_fname):
        self.fname = xml_fname
        self.root = xml.etree.ElementTree.parse(xml_fname).getroot()

    def get_fname(self):
        return self.fname

    def get_flowcell(self):
        fc = self.root[0].find('Flowcell').text
        if '-' in fc:
            # miseq often adds a bunch of leading zeros and a dash in front
            fc = fc.split('-')[1]
        assert 4 <= len(fc) <= 9
        return fc

    def get_rundate_american(self):
        rundate = self.root[0].find('Date').text
        if len(rundate) == 6:
            y, m, d = (rundate[0:2], rundate[2:4], rundate[4:6])
            y = '20' + y
        elif len(rundate) == 8:
            y, m, d = (rundate[0:4], rundate[4:6], rundate[6:8])
        else:
            raise Exception()
        return '%s/%s/%s' % (m, d, y)

    def get_rundate_iso(self):
        rundate = self.root[0].find('Date').text
        if len(rundate) == 6:
            y, m, d = (rundate[0:2], rundate[2:4], rundate[4:6])
            y = '20' + y
        elif len(rundate) == 8:
            y, m, d = (rundate[0:4], rundate[4:6], rundate[6:8])
        else:
            raise Exception()
        return '%s-%s-%s' % (y, m, d)

    def get_machine(self):
        return self.root[0].find('Instrument').text

    def get_read_structure(self):
        reads = []
        for x in self.root[0].find('Reads').findall('Read'):
            order = int(x.attrib['Number'])
            read = x.attrib['NumCycles'] + (x.attrib['IsIndexedRead'] == 'Y' and 'B' or 'T')
            reads.append((order, read))
        return ''.join([r for _, r in sorted(reads)])

    def num_reads(self):
        return sum(1 for x in self.root[0].find('Reads').findall('Read') if x.attrib['IsIndexedRead'] == 'N')

# ======================
# ***  SampleSheet   ***
# ======================


class SampleSheet(object):
    ''' A class that reads an Illumina SampleSheet.csv or alternative/simplified
        tab-delimited versions as well.
    '''

    def __init__(self, infile, use_sample_name=True, only_lane=None, allow_non_unique=False):
        self.fname = infile
        self.use_sample_name = use_sample_name
        if only_lane is not None:
            only_lane = str(only_lane)
        self.only_lane = only_lane
        self.allow_non_unique = allow_non_unique
        self.rows = []
        self._detect_and_load_sheet(infile)

    def _detect_and_load_sheet(self, infile):
        if infile.endswith('.csv'):
            # one of a few possible CSV formats (watch out for line endings from other OSes)
            with util.file.open_or_gzopen(infile, 'rU') as inf:
                header = None
                miseq_skip = False
                row_num = 0
                for line in inf:
                    row = line.rstrip('\n').split(',')
                    if miseq_skip:
                        if line.startswith('[Data]'):
                            # start paying attention *after* this line
                            miseq_skip = False
                        # otherwise, skip all the miseq headers
                    elif line.startswith('['):
                        # miseq: ignore all lines until we see "[Data]"
                        miseq_skip = True
                    elif header is None:
                        header = row
                        if 'Sample_ID' in header:
                            # this is a MiSeq-generated SampleSheet.csv
                            keymapper = {
                                'Sample_ID': 'sample',
                                'index': 'barcode_1',
                                'index2': 'barcode_2',
                                'Sample_Name': 'sample_name'
                            }
                            header = list(map(keymapper.get, header))
                        elif 'SampleID' in header:
                            # this is a Broad Platform generated SampleSheet.csv
                            keymapper = {
                                'SampleID': 'sample',
                                'Index': 'barcode_1',
                                'Index2': 'barcode_2',
                                'libraryName': 'library_id_per_sample',
                                'FCID': 'flowcell',
                                'Lane': 'lane'
                            }
                            header = list(map(keymapper.get, header))
                        elif len(row) == 3:
                            # hopefully this is a Broad walk-up submission sheet (_web_iww_htdocs_seq...)
                            header = ['sample', 'barcode_1', 'barcode_2']
                            if 'sample' not in row[0].lower():
                                # this is an actual data row! (no header exists in this file)
                                row_num += 1
                                self.rows.append({
                                    'sample': row[0],
                                    'barcode_1': row[1],
                                    'barcode_2': row[2],
                                    'row_num': str(row_num)
                                })
                        else:
                            raise Exception('unrecognized filetype: %s' % infile)
                        for h in ('sample', 'barcode_1'):
                            assert h in header
                    else:
                        # data rows
                        row_num += 1
                        assert len(header) == len(row)
                        row = dict((k, v) for k, v in zip(header, row) if k and v)
                        row['row_num'] = str(row_num)
                        if (self.only_lane is not None and row.get('lane') and self.only_lane != row['lane']):
                            continue
                        if row['sample'] and row['barcode_1']:
                            self.rows.append(row)
            # go back and re-shuffle miseq columns if use_sample_name applies
            if (self.use_sample_name and 'sample_name' in header and all(row.get('sample_name') for row in self.rows)):
                for row in self.rows:
                    row['library_id_per_sample'] = row['sample']
                    row['sample'] = row['sample_name']
            for row in self.rows:
                if 'sample_name' in row:
                    del row['sample_name']
        elif infile.endswith('.txt'):
            # our custom tab file format: sample, barcode_1, barcode_2, library_id_per_sample
            self.rows = []
            row_num = 0
            for row in util.file.read_tabfile_dict(infile):
                assert row.get('sample') and row.get('barcode_1')
                row_num += 1
                row['row_num'] = str(row_num)
                self.rows.append(row)
        else:
            raise Exception('unrecognized filetype: %s' % infile)

        if not self.rows:
            raise Exception('empty file')

        # populate library IDs, run IDs (ie BAM filenames)
        for row in self.rows:
            row['library'] = row['sample']
            if row.get('library_id_per_sample'):
                row['library'] += '.l' + row['library_id_per_sample']
            row['run'] = row['library']
        if len(set(row['run'] for row in self.rows)) != len(self.rows):
            if self.allow_non_unique:
                log.warn("non-unique library IDs in this lane")
                unique_count = {}
                for row in self.rows:
                    unique_count.setdefault(row['library'], 0)
                    unique_count[row['library']] += 1
                    row['run'] += '.r' + str(unique_count[row['library']])
            else:
                raise Exception('non-unique library IDs in this lane')

        # are we single or double indexed?
        if all(row.get('barcode_2') for row in self.rows):
            self.indexes = 2
        elif any(row.get('barcode_2') for row in self.rows):
            raise Exception('inconsistent single/double barcoding in sample sheet')
        else:
            self.indexes = 1

    def make_barcodes_file(self, outFile):
        ''' Create input file for Picard ExtractBarcodes '''
        if self.num_indexes() == 2:
            header = ['barcode_name', 'library_name', 'barcode_sequence_1', 'barcode_sequence_2']
        else:
            header = ['barcode_name', 'library_name', 'barcode_sequence_1']
        with open(outFile, 'wt') as outf:
            outf.write('\t'.join(header) + '\n')
            for row in self.rows:
                out = {
                    'barcode_sequence_1': row['barcode_1'],
                    'barcode_sequence_2': row.get('barcode_2', ''),
                    'barcode_name': row['sample'],
                    'library_name': row['library']
                }
                outf.write('\t'.join(out[h] for h in header) + '\n')

    def make_params_file(self, bamDir, outFile):
        ''' Create input file for Picard IlluminaBasecallsToXXX '''
        if self.num_indexes() == 2:
            header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1', 'BARCODE_2']
        else:
            header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1']
        with open(outFile, 'wt') as outf:
            outf.write('\t'.join(header) + '\n')
            # add one catchall entry at the end called Unmatched
            rows = self.rows + [{
                'barcode_1': 'N',
                'barcode_2': 'N',
                'sample': 'Unmatched',
                'library': 'Unmatched',
                'run': 'Unmatched'
            }]
            for row in rows:
                out = {
                    'BARCODE_1': row['barcode_1'],
                    'BARCODE_2': row.get('barcode_2', ''),
                    'SAMPLE_ALIAS': row['sample'],
                    'LIBRARY_NAME': row['library']
                }
                out['OUTPUT'] = os.path.join(bamDir, row['run'] + ".bam")
                outf.write('\t'.join(out[h] for h in header) + '\n')

    def get_fname(self):
        return self.fname

    def get_rows(self):
        return self.rows

    def num_indexes(self):
        ''' Return 1 or 2 depending on whether pools are single or double indexed '''
        return self.indexes

    def fetch_by_index(self, idx):
        idx = str(idx)
        for row in self.rows:
            if idx == row['row_num']:
                return row
        return None

# =============================
# ***  miseq_fastq_to_bam   ***
# =============================


def miseq_fastq_to_bam(outBam, sampleSheet, inFastq1, inFastq2=None, runInfo=None,
                       sequencing_center=None,
                       JVMmemory=tools.picard.FastqToSamTool.jvmMemDefault):
    ''' Convert fastq read files to a single bam file. Fastq file names must conform
        to patterns emitted by Miseq machines. Sample metadata must be provided
        in a SampleSheet.csv that corresponds to the fastq filename. Specifically,
        the _S##_ index in the fastq file name will be used to find the corresponding
        row in the SampleSheet
    '''

    # match miseq based on fastq filenames
    mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq1)
    assert mo, "fastq filename %s does not match the patterns used by an Illumina Miseq machine" % inFastq1
    assert mo.group(2) == '1', "fastq1 must correspond to read 1, not read %s" % mo.group(2)
    sample_num = mo.group(1)
    if inFastq2:
        mo = re.match(r"^\S+_S(\d+)_L001_R(\d)_001.fastq(?:.gz|)$", inFastq2)
        assert mo, "fastq filename %s does not match the patterns used by an Illumina Miseq machine" % inFastq2
        assert mo.group(2) == '2', "fastq2 must correspond to read 2, not read %s" % mo.group(2)
        assert mo.group(1) == sample_num, "fastq1 (%s) and fastq2 (%s) must have the same sample number" % (
            sample_num, mo.group(1))

    # load metadata
    samples = SampleSheet(sampleSheet, allow_non_unique=True)
    sample_info = samples.fetch_by_index(sample_num)
    assert sample_info, "sample %s not found in %s" % (sample_num, sampleSheet)
    sampleName = sample_info['sample']
    log.info("Using sample name: %s", sampleName)
    if sample_info.get('barcode_2'):
        barcode = '-'.join((sample_info['barcode_1'], sample_info['barcode_2']))
    else:
        barcode = sample_info['barcode_1']
    picardOpts = {
        'LIBRARY_NAME': sample_info['library'],
        'PLATFORM': 'illumina',
        'VERBOSITY': 'WARNING',
        'QUIET': 'TRUE',
    }
    if runInfo:
        runInfo = RunInfo(runInfo)
        flowcell = runInfo.get_flowcell()
        picardOpts['RUN_DATE'] = runInfo.get_rundate_iso()
        if inFastq2:
            assert runInfo.num_reads() == 2, "paired fastqs given for a single-end RunInfo.xml"
        else:
            assert runInfo.num_reads() == 1, "second fastq missing for a paired-end RunInfo.xml"
    else:
        flowcell = 'A'
    if sequencing_center is None and runInfo:
        sequencing_center = runInfo.get_machine()
    if sequencing_center:
        picardOpts['SEQUENCING_CENTER'] = sequencing_center
    picardOpts['PLATFORM_UNIT'] = '.'.join((flowcell, '1', barcode))
    if len(flowcell) > 5:
        flowcell = flowcell[:5]
    picardOpts['READ_GROUP_NAME'] = flowcell

    # run Picard
    picard = tools.picard.FastqToSamTool()
    picard.execute(inFastq1,
                   inFastq2,
                   sampleName,
                   outBam,
                   picardOptions=picard.dict_to_picard_opts(picardOpts),
                   JVMmemory=JVMmemory)
    return 0


def parser_miseq_fastq_to_bam(parser=argparse.ArgumentParser()):
    parser.add_argument('outBam', help='Output BAM file.')
    parser.add_argument('sampleSheet', help='Input SampleSheet.csv file.')
    parser.add_argument('inFastq1', help='Input fastq file; 1st end of paired-end reads if paired.')
    parser.add_argument('--inFastq2', help='Input fastq file; 2nd end of paired-end reads.', default=None)
    parser.add_argument('--runInfo', help='Input RunInfo.xml file.', default=None)
    parser.add_argument(
        '--sequencing_center',
        default=None,
        help='Name of your sequencing center (default is the sequencing machine ID from the RunInfo.xml)')
    parser.add_argument('--JVMmemory',
                        default=tools.picard.FastqToSamTool.jvmMemDefault,
                        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, miseq_fastq_to_bam, split_args=True)
    return parser


__commands__.append(('miseq_fastq_to_bam', parser_miseq_fastq_to_bam))


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

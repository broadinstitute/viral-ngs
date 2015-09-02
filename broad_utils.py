#!/usr/bin/env python
"""
Utilities for getting sequences out of the Broad walk-up sequencing pipeline.
These utilities are probably not of much use outside the Broad.
"""

__author__ = "dpark@broadinstitute.org"
__commands__ = []

import argparse
import logging
import os
import os.path
import shutil
import tempfile
import time
import json
import glob
import hashlib
import base64
import xml.etree.ElementTree
import util.cmd
import util.file
import tools.picard

log = logging.getLogger(__name__)

# ==========================================
# ***  get stuff from Picard json file   ***
# ==========================================


def get_json_from_picard(picardDir):
    ''' for example, /seq/walkup/picard/{flowcell_minus_first_char} '''
    analysisDir = max(
        (os.path.getmtime(os.path.join(picardDir, d)), d) for d in os.listdir(picardDir)
        if os.path.isdir(os.path.join(picardDir, d)))[1]
    jsonfile = list(glob.glob(os.path.join(picardDir, analysisDir, 'info', 'logs', '*.json')))
    if len(jsonfile) != 1:
        raise Exception("error")
    return jsonfile[0]

def get_run_date(jsonfile):
    with open(jsonfile, 'rt') as inf:
        runDate = json.load(inf)['workflow']['runDate']
    return runDate

def get_bustard_dir(jsonfile):
    with open(jsonfile, 'rt') as inf:
        bustard = json.load(inf)['workflow']['runFolder']
    return bustard


def parser_get_bustard_dir(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Picard directory')
    util.cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util.cmd.attach_main(parser, main_get_bustard_dir)
    return parser
def main_get_bustard_dir(args):
    'Find the basecalls directory from a Picard directory'
    print(get_bustard_dir(get_json_from_picard(args.inDir)))
    return 0
__commands__.append(('get_bustard_dir', parser_get_bustard_dir))


def parser_get_run_date(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Picard directory')
    util.cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util.cmd.attach_main(parser, main_get_run_date)
    return parser
def main_get_run_date(args):
    'Find the sequencing run date from a Picard directory'
    print(get_run_date(get_json_from_picard(args.inDir)))
    return 0
__commands__.append(('get_run_date', parser_get_run_date))


# ===============
# ***  misc   ***
# ===============

def iterate_lanes(runfile):
    for flowcellfile in util.file.read_tabfile(runfile):
        for lane in util.file.read_tabfile_dict(flowcellfile[0]):
            yield lane


def iterate_wells(runfile):
    for lane in iterate_lanes(runfile):
        for well in util.file.read_tabfile_dict(lane['barcode_file']):
            yield (lane, well)


def get_all_samples(runfile):
    return list(sorted(set(well['sample'] for lane, well in iterate_wells(runfile))))


def get_all_libraries(runfile):
    return list(sorted(set(well['sample'] + '.l' + well['library_id_per_sample'] for lane, well in iterate_wells(
        runfile))))


def get_run_id(well):
    run_id = well['sample']
    if well.get('library_id_per_sample'):
        run_id += '.l' + well['library_id_per_sample']
    if well.get('run_id_per_library'):
        run_id += '.r' + well['run_id_per_library']
    return run_id


def get_all_runs(runfile):
    return list(sorted(get_run_id(well) + '.' + lane['flowcell'] + '.' + lane['lane'] for lane, well in iterate_wells(
        runfile)))


def parser_get_all_names(parser=argparse.ArgumentParser()):
    parser.add_argument('type', help='Type of name', choices=['samples', 'libraries', 'runs'])
    parser.add_argument('runfile', help='File with seq run information')
    util.cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util.cmd.attach_main(parser, main_get_all_names)
    return parser
def main_get_all_names(args):
    'Get all samples'
    if args.type == 'samples':
        method = get_all_samples
    elif args.type == 'libraries':
        method = get_all_libraries
    elif args.type == 'runs':
        method = get_all_runs
    for s in method(args.runfile):
        print(s)
    return 0
__commands__.append(('get_all_names', parser_get_all_names))



# =========================
# ***  illumina_demux   ***
# =========================

def make_barcodes_file(inFile, outFile):
    'Create input file for extract_barcodes'
    if any(row.get('barcode_2') for row in util.file.read_tabfile_dict(inFile)):
        header = ['barcode_name', 'library_name', 'barcode_sequence_1', 'barcode_sequence_2']
    else:
        header = ['barcode_name', 'library_name', 'barcode_sequence_1']
    with open(outFile, 'wt') as outf:
        outf.write('\t'.join(header) + '\n')
        for row in util.file.read_tabfile_dict(inFile):
            out = {
                'barcode_sequence_1': row['barcode_1'],
                'barcode_sequence_2': row.get('barcode_2', ''),
                'barcode_name': row['sample'],
                'library_name': row['sample']
            }
            if row.get('library_id_per_sample'):
                out['library_name'] += '.l' + row['library_id_per_sample']
            outf.write('\t'.join(out[h] for h in header) + '\n')

def make_params_file(inFile, bamDir, outFile):
    'Create input file for illumina_basecalls'
    if any(row.get('barcode_2') for row in util.file.read_tabfile_dict(inFile)):
        header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1', 'BARCODE_2']
    else:
        header = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1']
    with open(outFile, 'wt') as outf:
        outf.write('\t'.join(header) + '\n')
        rows = list(util.file.read_tabfile_dict(inFile))
        rows.append({'barcode_1': 'N', 'barcode_2': 'N', 'sample': 'Unmatched'})
        for row in rows:
            out = {
                'BARCODE_1': row['barcode_1'],
                'BARCODE_2': row.get('barcode_2', ''),
                'SAMPLE_ALIAS': row['sample'],
                'LIBRARY_NAME': row['sample']
            }
            if row.get('library_id_per_sample'):
                out['LIBRARY_NAME'] += '.l' + row['library_id_per_sample']
            run_id = out['LIBRARY_NAME']
            if row.get('run_id_per_library'):
                run_id += '.r' + row['run_id_per_library']
            out['OUTPUT'] = os.path.join(bamDir, run_id + ".bam")
            outf.write('\t'.join(out[h] for h in header) + '\n')

def get_earliest_date(inDir):
    ''' Looks at the dates of all first-level members of a directory, plus the directory
        itself, and returns the earliest date seen.
    '''
    fnames = [inDir] + [os.path.join(inDir, x) for x in os.listdir(inDir)]
    earliest = min(os.path.getmtime(fn) for fn in fnames)
    # return time.strftime("%Y-%m-%d", time.localtime(earliest))
    # what?? http://sourceforge.net/p/samtools/mailman/message/27441767/
    return time.strftime("%m/%d/%Y", time.localtime(earliest))


def parser_illumina_demux(parser=argparse.ArgumentParser()):
    parser.add_argument('sampleSheet',
                        help='''Input tab file w/header and four named columns:
                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2''')
    parser.add_argument('bclDir', help='Illumina BCL directory (or tar.gz of BCL directory).')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('outDir', help='Output directory for BAM files.')
    
    parser.add_argument('--outMetrics', help='Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file.', default=None)
    parser.add_argument('--flowcell', help='Override flowcell ID (default: read from RunInfo.xml).', default=None)
    parser.add_argument('--read_structure', help='Override read structure (default: read from RunInfo.xml).', default=None)
    
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
    illumina = IlluminaDirectory(args.bclDir)
    illumina.load()
    if args.flowcell:
        flowcell = args.flowcell
    else:
        flowcell = illumina.get_RunInfo().get_flowcell()
    if args.run_start_date:
        run_date = args.run_start_date
    else:
        run_date = illumina.get_RunInfo().get_rundate()
    if args.read_structure:
        read_structure = args.read_structure
    else:
        read_struture = illumina.get_RunInfo().get_read_structure()
    
    # Picard ExtractIlluminaBarcodes
    extract_input = util.file.mkstempfname('.txt', prefix='.'.join(['barcodeData', flowcell, str(args.lane)]))
    barcodes_tmpdir = tempfile.mkdtemp(prefix='extracted_barcodes-')
    make_barcodes_file(args.sampleSheet, extract_input)
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
    make_params_file(args.sampleSheet, args.outDir, basecalls_input)
    picardOpts = dict((opt, getattr(args, opt)) for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    picardOpts['run_start_date'] = run_date
    picardOpts['read_structure'] = read_structure
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
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return 0
    
    def load(self):
        if self.path is None:
            if '://' in self.uri:
                tarball = util.file.mkstempfname('.tar.gz')
                raise NotImplementedError('boto s3 download here uri -> tarball')
                self._extract_tarball(tarball)
                os.unlink(tarball)
            else:
                if os.path.isdir(self.uri):
                    self.path = self.uri
                else:
                    self._extract_tarball(self.uri)
            self._fix_path()
    
    def _fix_path(self):
        if not os.path.isdir(os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')):
            # this is not the correct root-level directory
            # sometimes this points to one level up
            subdirs = list(x for x in os.listdir(self.path) if os.path.isdir(x))
            if len(subdirs) != 1:
                raise Exception('cannot find Data/Intensities/BaseCalls/ inside %s' % self.uri)
            self.path = subdirs[0]
            if not os.path.isdir(os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')):
                raise Exception('cannot find Data/Intensities/BaseCalls/ inside %s' % self.uri)
    
    def _extract_tarball(self, tarfile):
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
        subprocess.check_call(untar_cmd)
    
    def close(self):
        if self.tempDir:
            shutil.rmtree(self.tempDir)
            self.tempDir = None
    
    def get_RunInfo(self):
        if self.runinfo is None:
            self.runinfo = RunInfo(os.path.join(self.path, 'RunInfo.xml'))
        return self.runinfo
    
    def get_SampleSheet(self):
        if self.samplesheet is None:
            self.samplesheet = SampleSheet(os.path.join(self.path, 'SampleSheet.csv'))
        return self.samplesheet
    
    def get_BCLdir(self):
        return os.path.join(self.path, 'Data', 'Intensities', 'BaseCalls')


class RunInfo(object):
    ''' A class that reads the RunInfo.xml file emitted by Illumina
        MiSeq and HiSeq machines.
    '''
    def __init__(self, xml_fname):
        self.root = xml.etree.ElementTree.parse(xml_fname).getroot()
    
    def get_flowcell(self):
        fc = self.root[0].find('Flowcell').text
        if '-' in fc:
            # miseq often adds a bunch of leading zeros and a dash in front
            fc = fc.split('-')[1]
        assert 4 <= len(fc) <= 9
        return fc
    
    def get_rundate(self):
        rundate = self.root[0].find('Date').text
        assert len(rundate) == 6
        return '%s/%s/20%s' % (rundate[2:4], rundate[4:6], rundate[0:2])
        
    def get_machine(self):
        return self.root[0].find('Instrument').text
    
    def get_read_structure(self):
        reads = []
        for x in self.root[0].find('Reads').findall('Read'):
            order  = int(x.attrib['Number'])
            read  = x.attrib['NumCycles'] + (x.attrib['IsIndexedRead'] == 'Y' and 'B' or 'T')
            reads.append((order, read))
        return ''.join([r for o,r in sorted(reads)])

class SampleSheet(object):
    ''' A class that reads an Illumina SampleSheet.csv or alternative/simplified
        tab-delimited versions as well.
    '''
    def __init__(self, infile):
        pass


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

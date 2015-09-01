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
import time
import json
import glob
import hashlib
import base64
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
        (os.path.getmtime(os.path.join(picardDir, d)), d)
        for d in os.listdir(picardDir)
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
    return list(sorted(set(well['sample']
                           for lane, well in iterate_wells(runfile))))


def get_all_libraries(runfile):
    return list(sorted(set(well['sample'] + '.l' + well['library_id_per_sample']
                           for lane, well in iterate_wells(runfile))))


def get_run_id(well):
    run_id = well['sample']
    if well.get('library_id_per_sample'):
        run_id += '.l' + well['library_id_per_sample']
    if well.get('run_id_per_library'):
        run_id += '.r' + well['run_id_per_library']
    return run_id


def get_all_runs(runfile):
    return list(sorted(get_run_id(well) + '.' + lane['flowcell'] + '.' + lane['lane']
                       for lane, well in iterate_wells(runfile)))


def parser_get_all_names(parser=argparse.ArgumentParser()):
    parser.add_argument('type', help='Type of name',
                        choices=['samples', 'libraries', 'runs'])
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


# =============================
# ***  make_barcodes_file   ***
# =============================

def make_barcodes_file(inFile, outFile):
    'Create input file for extract_barcodes'
    if any(row.get('barcode_2') for row in util.file.read_tabfile_dict(inFile)):
        header = ['barcode_name', 'library_name', 'barcode_sequence_1', 'barcode_sequence_2']
    else:
        header = ['barcode_name', 'library_name', 'barcode_sequence_1']
    with open(outFile, 'wt') as outf:
        outf.write('\t'.join(header) + '\n')
        for row in util.file.read_tabfile_dict(inFile):
            out = {'barcode_sequence_1': row['barcode_1'],
                   'barcode_sequence_2': row.get('barcode_2', ''),
                   'barcode_name': row['sample'],
                   'library_name': row['sample']}
            if row.get('library_id_per_sample'):
                out['library_name'] += '.l' + row['library_id_per_sample']
            outf.write('\t'.join(out[h] for h in header) + '\n')


def parser_make_barcodes_file(parser=argparse.ArgumentParser()):
    parser.add_argument('inFile',
                        help='''Input tab file w/header and 2-5 named columns (last three are optional):
                sample, barcode_1, barcode_2, library_id_per_sample, run_id_per_library''')
    parser.add_argument('outFile', help='Output BARCODE_FILE file for Picard.')
    util.cmd.attach_main(parser, make_barcodes_file, split_args=True)
    return parser
__commands__.append(('make_barcodes_file', parser_make_barcodes_file))

# ===========================
# ***  extract_barcodes   ***
# ===========================


def parser_extract_barcodes(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Bustard directory.')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('barcodeFile',
                        help='''Input tab file w/header and four named columns:
                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2''')
    parser.add_argument('outDir', help='Output directory for barcodes.')
    parser.add_argument('--outMetrics',
                        help='Output metrics file. Default is to dump to a temp file.',
                        default=None)
    for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list:
        parser.add_argument('--' + opt,
                            help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                            default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt))
    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_extract_barcodes)
    return parser


def main_extract_barcodes(args):
    'Match every read in a lane against their barcode.'
    out_metrics = (args.outMetrics is None) and util.file.mkstempfname('.metrics.txt') or args.outMetrics
    picardOpts = dict((opt, getattr(args, opt))
                      for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    tools.picard.ExtractIlluminaBarcodesTool().execute(
        os.path.join(args.inDir, 'Data', 'Intensities', 'BaseCalls'), args.lane, args.barcodeFile,
        args.outDir, out_metrics,
        picardOptions=picardOpts, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('extract_barcodes', parser_extract_barcodes))

# ==============================
# ***  make_library_params   ***
# ==============================


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
            out = {'BARCODE_1': row['barcode_1'],
                   'BARCODE_2': row.get('barcode_2', ''),
                   'SAMPLE_ALIAS': row['sample'],
                   'LIBRARY_NAME': row['sample']}
            if row.get('library_id_per_sample'):
                out['LIBRARY_NAME'] += '.l' + row['library_id_per_sample']
            run_id = out['LIBRARY_NAME']
            if row.get('run_id_per_library'):
                run_id += '.r' + row['run_id_per_library']
            out['OUTPUT'] = os.path.join(bamDir, run_id + ".bam")
            outf.write('\t'.join(out[h] for h in header) + '\n')


def parser_make_params_file(parser=argparse.ArgumentParser()):
    parser.add_argument('inFile',
                        help='''Input tab file w/header and four named columns:
                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2''')
    parser.add_argument('bamDir', help='Directory for output bams')
    parser.add_argument('outFile', help='Output LIBRARY_PARAMS file for Picard')
    util.cmd.attach_main(parser, make_params_file, split_args=True)
    return parser
__commands__.append(('make_params_file', parser_make_params_file))

# =============================
# ***  illumina_basecalls   ***
# =============================


def get_earliest_date(inDir):
    ''' Looks at the dates of all first-level members of a directory, plus the directory
        itself, and returns the earliest date seen.
    '''
    fnames = [inDir] + [os.path.join(inDir, x) for x in os.listdir(inDir)]
    earliest = min(os.path.getmtime(fn) for fn in fnames)
    # return time.strftime("%Y-%m-%d", time.localtime(earliest))
    # what?? http://sourceforge.net/p/samtools/mailman/message/27441767/
    return time.strftime("%m/%d/%Y", time.localtime(earliest))


def short_hash(inString, length=None):
    ''' Returns a base32 encoding of a SHA1 hash of the inString, optionally truncated
        to a maximum length.  The base32 encoding is uppercase A-Z and 2-7.
    '''
    hash_obj = hashlib.sha1(inString.encode('utf-8'))
    b32_str = base64.b32encode(bytes(hash_obj.digest())).decode('utf-8')
    if length > 0 and len(b32_str) > length:
        b32_str = b32_str[:length]
    return b32_str


def parser_illumina_basecalls(parser=argparse.ArgumentParser()):
    parser.add_argument('inBustardDir', help='Bustard directory.')
    parser.add_argument('inBarcodesDir', help='Barcodes directory.')
    parser.add_argument('flowcell', help='Flowcell ID')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('paramsFile',
                        help='''Input tab file w/header and five named columns:
                BARCODE_1, BARCODE_2, OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME''')
    for opt in tools.picard.IlluminaBasecallsToSamTool.option_list:
        if opt == 'adapters_to_check':
            parser.add_argument('--' + opt, nargs='*',
                                help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
        else:
            parser.add_argument('--' + opt,
                                help='Picard ExtractIlluminaBarcodes ' + opt.upper() + ' (default: %(default)s)',
                                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
    parser.add_argument('--JVMmemory',
                        help='JVM virtual memory size (default: %(default)s)',
                        default=tools.picard.IlluminaBasecallsToSamTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    util.cmd.attach_main(parser, main_illumina_basecalls)
    return parser


def main_illumina_basecalls(args):
    'Demultiplex Illumina runs & produce BAM files, one per sample'
    picardOpts = dict((opt, getattr(args, opt))
                      for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
                      if hasattr(args, opt) and getattr(args, opt) != None)
    if not picardOpts.get('run_start_date'):
        picardOpts['run_start_date'] = get_earliest_date(args.inBustardDir)
    # if not picardOpts.get('read_group_id'):
    #    picardOpts['read_group_id'] = short_hash('{}.{}'.format(args.flowcell,args.lane), 6)
    tools.picard.IlluminaBasecallsToSamTool().execute(
        os.path.join(args.inBustardDir, 'Data/Intensities/BaseCalls'),
        args.inBarcodesDir, args.flowcell, args.lane, args.paramsFile,
        picardOptions=picardOpts, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('illumina_basecalls', parser_illumina_basecalls))


# =======================

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

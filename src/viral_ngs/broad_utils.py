#!/usr/bin/env python
"""
Utilities for getting sequences out of the Broad walk-up sequencing pipeline.
"""

__author__ = "dpark@broadinstitute.org"
__version__ = "PLACEHOLDER"
__commands__ = []

import argparse, logging, os, os.path, time, json, glob, hashlib
import util.cmd, util.file
import tools.picard

log = logging.getLogger(__name__)


# ==========================
# ***  get_bustard_dir   ***
# ==========================

def get_bustard_dir(picardDir):
    ''' for example, /seq/walkup/picard/{flowcell_minus_first_char} '''
    analysisDir = max(
        (os.path.getmtime(os.path.join(picardDir, d)), d)
        for d in os.listdir(picardDir)
        if os.path.isdir(os.path.join(picardDir, d)))[1]
    jsonfile = list(glob.glob(os.path.join(picardDir, analysisDir, 'info', 'logs', '*.json')))
    if len(jsonfile) != 1:
        raise Exception("error")
    jsonfile = jsonfile[0]
    with open(jsonfile, 'rt') as inf:
        bustard = json.load(inf)['workflow']['runFolder']
    return bustard
def parser_get_bustard_dir() :
    parser = argparse.ArgumentParser(
        description='Find the Bustard directory from a Picard directory')
    parser.add_argument('inDir',  help='Picard directory')
    util.cmd.common_args(parser, (('loglevel', 'ERROR'),))
    return parser
def main_get_bustard_dir(args) :
    print(get_bustard_dir(args.inDir))
    return 0
__commands__.append(('get_bustard_dir', main_get_bustard_dir, parser_get_bustard_dir))


# =============================
# ***  make_barcodes_file   ***
# =============================

def make_barcodes_file(inFile, outFile):
    header = ['barcode_name', 'library_name', 'barcode_sequence_1', 'barcode_sequence_2']
    with open(outFile, 'wt') as outf:
        outf.write('\t'.join(header)+'\n')
        for row in util.file.read_tabfile_dict(inFile):
            out  = {'barcode_sequence_1':row['barcode_1'],
                    'barcode_sequence_2':row['barcode_2'],
                    'barcode_name':row['sample'],
                    'library_name':row['sample']}
            if row.get('library_id_per_sample'):
                out['library_name'] += '.l' + row['library_id_per_sample']
            outf.write('\t'.join(out[h] for h in header)+'\n')
def parser_make_barcodes_file() :
    parser = argparse.ArgumentParser(description='Create input file for extract_barcodes')
    parser.add_argument('inFile',
        help='''Input tab file w/header and 3-5 named columns (last two are optional):
                sample, barcode_1, barcode_2, library_id_per_sample, run_id_per_library''')
    parser.add_argument('outFile', help='Output BARCODE_FILE file for Picard.')
    return parser
def main_make_barcodes_file(args) :
    make_barcodes_file(args.inFile, args.outFile)
    return 0
__commands__.append(('make_barcodes_file', main_make_barcodes_file, parser_make_barcodes_file))

# ===========================
# ***  extract_barcodes   ***
# ===========================

def parser_extract_barcodes():
    parser = argparse.ArgumentParser(
        description='''Match every read in a lane against their barcode.''')
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
        parser.add_argument('--'+opt,
            help='Picard ExtractIlluminaBarcodes '+opt.upper()+' (default: %(default)s)',
            default=tools.picard.ExtractIlluminaBarcodesTool.defaults.get(opt))
    parser.add_argument('--JVMmemory',
        help='JVM virtual memory size (default: %(default)s)',
        default = tools.picard.ExtractIlluminaBarcodesTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_extract_barcodes(args):
    out_metrics = (args.outMetrics==None) and util.file.mkstempfname('.metrics.txt') or args.outMetrics
    picardOpts = dict((opt, getattr(args, opt))
        for opt in tools.picard.ExtractIlluminaBarcodesTool.option_list
        if hasattr(args, opt) and getattr(args, opt)!=None)
    tools.picard.ExtractIlluminaBarcodesTool().execute(
        os.path.join(args.inDir, 'Data/Intensities/BaseCalls'), args.lane, args.barcodeFile,
        args.outDir, out_metrics,
        picardOptions=picardOpts, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('extract_barcodes', main_extract_barcodes, parser_extract_barcodes))

# ==============================
# ***  make_library_params   ***
# ==============================

def make_params_file(inFile, flowcell, lane, bamDir, outFile):
    header = ['BARCODE_1', 'BARCODE_2', 'OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'ID']
    with open(outFile, 'wt') as outf:
        outf.write('\t'.join(header)+'\n')
        rows = list(util.file.read_tabfile_dict(inFile))
        rows.append({'barcode_1':'N','barcode_2':'N','sample':'Unmatched'})
        for row in rows:
            out  = {'BARCODE_1':row['barcode_1'],
                    'BARCODE_2':row['barcode_2'],
                    'SAMPLE_ALIAS':row['sample'],
                    'LIBRARY_NAME':row['sample']}
            if row.get('library_id_per_sample'):
                out['LIBRARY_NAME'] += '.l' + row['library_id_per_sample']
            run_id = out['LIBRARY_NAME']
            if row.get('run_id_per_library'):
                run_id += '.r' + row['run_id_per_library']
            out['OUTPUT'] = os.path.join(bamDir, run_id + ".bam")
            out['ID'] = hashlib.sha1('{}.{}.{}'.format(run_id,flowcell,lane).encode('utf-8')).hexdigest().upper()[:8]
            outf.write('\t'.join(out[h] for h in header)+'\n')
def parser_make_params_file():
    parser = argparse.ArgumentParser(description='Create input file for illumina_basecalls')
    parser.add_argument('barcodeFile',
        help='''Input tab file w/header and four named columns:
                barcode_name, library_name, barcode_sequence_1, barcode_sequence_2''')
    parser.add_argument('flowcell', help='Flowcell ID')
    parser.add_argument('lane', help='Lane number', type=int)
    parser.add_argument('bamDir', help='Directory for output bams')
    parser.add_argument('outFile', help='Output LIBRARY_PARAMS file for Picard')
    return parser
def main_make_params_file(args):
    make_params_file(args.barcodeFile, args.flowcell, args.lane, args.bamDir, args.outFile)
    return 0
__commands__.append(('make_params_file', main_make_params_file, parser_make_params_file))

# =============================
# ***  illumina_basecalls   ***
# =============================

def get_earliest_date(inDir):
    ''' Looks at the dates of all first-level members of a directory, plus the directory
        itself, and returns the earliest date seen.
    '''
    fnames = [inDir] + [os.path.join(inDir, x) for x in os.listdir(inDir)]
    earliest = min(os.path.getmtime(fn) for fn in fnames)
    return time.strftime("%Y-%m-%d", time.localtime(earliest))
    
def parser_illumina_basecalls():
    parser = argparse.ArgumentParser(
        description='''Demultiplex Illumina runs & produce BAM files, one per sample''')
    parser.add_argument('inBustardDir', help='Bustard directory.')
    parser.add_argument('inBarcodesDir', help='Barcodes directory.')
    parser.add_argument('flowcell', help='Flowcell ID')
    parser.add_argument('lane', help='Lane number.', type=int)
    parser.add_argument('paramsFile',
        help='''Input tab file w/header and five named columns:
                BARCODE_1, BARCODE_2, OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME''')
    for opt in tools.picard.IlluminaBasecallsToSamTool.option_list:
        if opt=='adapters_to_check':
            parser.add_argument('--'+opt, nargs='*',
                help='Picard ExtractIlluminaBarcodes '+opt.upper()+' (default: %(default)s)',
                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
        else:
            parser.add_argument('--'+opt,
                help='Picard ExtractIlluminaBarcodes '+opt.upper()+' (default: %(default)s)',
                default=tools.picard.IlluminaBasecallsToSamTool.defaults.get(opt))
    parser.add_argument('--JVMmemory',
        help='JVM virtual memory size (default: %(default)s)',
        default = tools.picard.IlluminaBasecallsToSamTool.jvmMemDefault)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmpDir', None)))
    return parser
def main_illumina_basecalls(args):
    picardOpts = dict((opt, getattr(args, opt))
        for opt in tools.picard.IlluminaBasecallsToSamTool.option_list
        if hasattr(args, opt) and getattr(args, opt)!=None and opt != 'read_structure')
    params_file = util.file.mkstempfname('library_params.txt')
    if not picardOpts.get('run_start_date'):
        picardOpts['run_start_date'] = get_earliest_date(args.inIlluminaDir)
    tools.picard.IlluminaBasecallsToSamTool().execute(
        os.path.join(args.inBustardDir, 'Data/Intensities/BaseCalls'),
        args.lane, args.read_structure, args.inBarcodesDir,
        args.flowcell, args.paramsFile,
        picardOptions=picardOpts, JVMmemory=args.JVMmemory)
    return 0
__commands__.append(('illumina_basecalls',
    main_illumina_basecalls, parser_illumina_basecalls))



# =======================

if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

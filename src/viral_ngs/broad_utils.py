#!/usr/bin/env python3
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
import json
import glob
from .core import file as util_file, misc as util_misc, cmd as util_cmd

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
    util_cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util_cmd.attach_main(parser, main_get_bustard_dir)
    return parser
def main_get_bustard_dir(args):
    'Find the basecalls directory from a Picard directory'
    print(get_bustard_dir(get_json_from_picard(args.inDir)))
    return 0
__commands__.append(('get_bustard_dir', parser_get_bustard_dir))


def parser_get_run_date(parser=argparse.ArgumentParser()):
    parser.add_argument('inDir', help='Picard directory')
    util_cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util_cmd.attach_main(parser, main_get_run_date)
    return parser
def main_get_run_date(args):
    'Find the sequencing run date from a Picard directory'
    print(get_run_date(get_json_from_picard(args.inDir)))
    return 0
__commands__.append(('get_run_date', parser_get_run_date))


# ===============
# ***  misc   ***
# ===============

def iterate_wells(runfile):
    for lane in util_file.read_tabfile_dict(runfile):
        for well in util_file.read_tabfile_dict(lane['barcode_file']):
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
    util_cmd.common_args(parser, (('loglevel', 'ERROR'),))
    util_cmd.attach_main(parser, main_get_all_names)
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


# =======================
def full_parser():
    return util_cmd.make_parser(__commands__, __doc__)


def main():
    util_cmd.main_argparse(__commands__, __doc__)


if __name__ == '__main__':
    main()

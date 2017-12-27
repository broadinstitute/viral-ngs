#!/usr/bin/env python
''' Utilities for managing data.
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import logging
import random
import numpy
import os
import os.path
import shutil
import subprocess
import functools
import operator
import concurrent.futures
import csv
import json

try:
    from itertools import zip_longest    # pylint: disable=E0611
except ImportError:
    from itertools import izip_longest as zip_longest    # pylint: disable=E0611

# intra-module
import util.cmd
import util.file
import util.misc
from util.metadata import InFile, OutFile

# third-party
import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.IUPACData
import networkx
import networkx.algorithms.operators.binary

log = logging.getLogger(__name__)


def _load_provenance_graph():
    """Read the provenance graph"""


    pgraph = networkx.DiGraph()
    
    for step_graph_fname in os.listdir(util.metadata.metadata_dir()):
        print(step_graph_fname)
        if step_graph_fname.endswith('.json'):
            step_graph = json.loads(util.file.slurp_file(os.path.join(util.metadata.metadata_dir(), step_graph_fname)))
            step_id = step_graph['step']['step_id']
            pgraph.add_node(step_id, **step_graph['step'])
            pgraph.add_node(step_id, node_kind='step')
            for file_kind in ('in_files', 'out_files'):
                for port_name, file_info in step_graph[file_kind].items():
                    port_id = step_id + '__' + port_name
                    pgraph.add_node(port_id, **file_info)
                    pgraph.add_node(port_id, node_kind=('input' if file_kind=='in_files' else 'output'))
                    pgraph.add_edge(*((port_id, step_id) if file_kind=='in_files' else (step_id, port_id)))
                    if 'hash' in file_info:
                        file_id = file_info['hash']
                        if file_id:
                            pgraph.add_node(file_id, node_kind='file')
                            pgraph.add_edge(*((file_id, port_id) if file_kind=='in_files' else (port_id, file_id)))

    return pgraph


def plot_provenance(prov_plot):
    """Plot the provenance graph"""

    print('saving to ', prov_plot)

    pgraph = _load_provenance_graph()
    
    with open(prov_plot, 'wt') as out:
        out.write('digraph G {\n')
        def fix_n(n):
            return 'n'+n.replace('-','_')
        for n, d in pgraph.nodes(data=True):
            print('node n={} data d={}'.format(n, d))
            if d['node_kind'] == 'file':
                out.write('{} [shape=point];\n'.format(fix_n(n)))
            elif d['node_kind'] in ('input, output'):
                out.write('{} [shape=ellipse, label="{}"];\n'.format(fix_n(n), d['fname']))
            else:
                out.write('{} [shape=invtrapezium, label="{}"];\n'.format(fix_n(n), d['cmd_name']))
        for n1, n2, d in pgraph.edges(data=True):
            out.write('{} -> {};\n'.format(fix_n(n1), fix_n(n2)))
                
        out.write('}\n')
            
def parser_plot_provenance(parser=argparse.ArgumentParser()):
    parser.add_argument('prov_plot', help='Where to save the provenance plot')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, plot_provenance, split_args=True)
    return parser

__commands__.append(('plot_provenance', parser_plot_provenance))

###############################################


def compute_paths(paths_file):
    """Compute computation paths, and their attributes."""

    pgraph = networkx.DiGraph()
    
    for step_graph_fname in os.listdir(util.metadata.metadata_dir()):
        print(step_graph_fname)
        if step_graph_fname.endswith('.json'):
            step_graph = json.loads(util.file.slurp_file(os.path.join(util.metadata.metadata_dir(), step_graph_fname)))
            step_id = step_graph['step']['step_id']
            pgraph.add_node(step_id, **step_graph['step'])
            pgraph.add_node(step_id, node_kind='step')
            for file_kind in ('in_files', 'out_files'):
                for file_hash, file_attrs in step_graph[file_kind].items():
                    pgraph.add_node(file_hash, **file_attrs)
                    pgraph.add_node(file_hash, node_kind='file')
                    pgraph.add_edge(*((file_hash, step_id) if file_kind=='in_files' else (step_id, file_hash)), arg=file_attrs['arg'])

    with open(prov_plot, 'wt') as out:
        out.write('digraph G {\n')
        def fix_n(n):
            return 'n'+n.replace('-','_')
        for n, d in pgraph.nodes(data=True):
            if d['node_kind'] == 'file':
                out.write('{} [shape=ellipse, label="{}"];\n'.format(fix_n(n), d['fname']))
            else:
                out.write('{} [shape=invtrapezium, label="{}"];\n'.format(fix_n(n), d['cmd_name']))
        for n1, n2, d in pgraph.edges(data=True):
            out.write('{} -> {} [label="{}"];\n'.format(fix_n(n1), fix_n(n2), d['arg']))
                
        out.write('}\n')
            
def parser_compute_paths(parser=argparse.ArgumentParser()):
    parser.add_argument('paths_filel', help='Where to save the paths')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, compute_paths, split_args=True)
    return parser

__commands__.append(('compute_paths', parser_compute_paths))


###############################################

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)



if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

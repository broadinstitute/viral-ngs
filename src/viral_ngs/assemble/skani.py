'''
    SKANI - accurate, fast nucleotide identity calculation for MAGs and databases
    https://github.com/bluenote-1577/skani
'''

__author__ = "dpark@broadinstitute.org"

import logging
import tools
import util.file
import util.misc
import csv
import os
import os.path
import shutil
import subprocess

TOOL_NAME = "skani"

_log = logging.getLogger(__name__)

class UndirectedGraph:
    ''' Simple utility class for finding clusters from pairwise relationships
    '''
    def __init__(self):
        self.edges = {}

    def add_node(self, node):
        self.edges.setdefault(node, set())

    def add_edge(self, node1, node2):
        self.edges.setdefault(node1, set()).add(node2)
        self.edges.setdefault(node2, set()).add(node1)

    def _dfs(self, node, visited):
        visited.add(node)
        cluster = set()
        cluster.add(node)
        for neighbor in self.edges[node]:
            if neighbor not in visited:
                cluster.update(self._dfs(neighbor, visited))
        return cluster

    def get_clusters(self):
        visited = set()
        clusters = []
        for node in self.edges.keys():
            if node not in visited:
                clusters.append(self._dfs(node, visited))
        return clusters


class SkaniTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(SkaniTool, self).__init__(install_methods=install_methods)

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '--version']).decode('UTF-8').strip().split()[1]

    def execute(self, subcommand, args, outfile, threads=None):
        ''' generic execution of skani
        '''

        # build the skani command
        tool_cmd = [self.install_and_get_path(), subcommand]
        tool_cmd.extend(args)
        tool_cmd.extend(["-t", "{}".format(util.misc.sanitize_thread_count(threads))])

        # run the command
        _log.debug(' '.join(tool_cmd) + ' > ' + outfile)
        with open(outfile, 'w') as outf:
            util.misc.run_and_save(tool_cmd, outf=outf)

    def triangle(self, ref_fastas, outfile_ani, outfile_af, other_args = (), threads=None):
        ''' skani triangle computes an all-to-all ANI distance matrix for a set of sequences
        '''
        self.execute('triangle', list(ref_fastas) + list(other_args), outfile_ani, threads=threads)
        shutil.copyfile('skani_matrix.af', outfile_af)

    def dist(self, query_fasta, ref_fastas, outfile, other_args = (), threads=None):
        ''' skani dist computes ANI distance between a specified query set of
            sequences (MAGs) and reference genomes (database)
        '''
        self.execute('dist', ['-q', query_fasta, '-r'] + list(ref_fastas) + list(other_args), outfile, threads=threads)

    def find_reference_clusters(self, ref_fastas,
                                other_args = ('-m', 50, '--no-learned-ani', '--slow', '--robust', '--detailed', '--ci', '--sparse'),
                                threads=None):
        ''' use skani triangle to define clusters of highly-related genomes
            (default settings here are for viral genomes)
        '''
        g = UndirectedGraph()
        for ref_fasta in ref_fastas:
            g.add_node(ref_fasta)

        with util.file.tempfnames(('.skani_matrix.ani', '.skani_matrix.af')) \
            as (tmp_matrix_ani, tmp_matrix_af):
            # run skani triangle
            self.triangle(ref_fastas, 'skani_matrix.ani', 'skani_matrix.af', other_args, threads=threads)

            # parse the skani triangle results and define clusters
            with open(tmp_matrix_ani, 'r') as inf:
                for row in csv.DictReader(inf, delimiter='\t'):
                    g.add_edge(row['Ref_file'], row['Query_file'])

        return g.get_clusters()

    def find_closest_reference(self, contigs_fasta, ref_fastas, out_file,
                                other_args = ('-m', 50, '--no-learned-ani', '--slow', '--robust', '--detailed', '--ci', '-s', 85, '-n', 10, '--no-marker-index'),
                                threads=None):
        ''' use skani dist to find the closest reference genome for each contig
            (default settings here are for viral genomes)
        '''
        self.dist(contigs_fasta, ref_fastas, out_file, other_args, threads=threads)
        with open(out_file, 'r') as inf:
            top_row = None
            for row in csv.DictReader(inf, delimiter='\t'):
                top_row = row
                break
        return (top_row['Ref_file'], top_row['ANI'], top_row['Align_fraction_ref'], top_row['Total_bases_covered']) if top_row is not None else None

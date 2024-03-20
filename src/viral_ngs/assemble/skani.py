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
import Bio.SeqIO

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

    def _is_fasta_basically_empty(self, inFasta, min_length=500):
        ''' Check if a fasta file contains no sequences longer than min_length.
            If the fasta contains only sequences < min_length, return True
        '''
        with open(inFasta, 'r') as inf:
            for seq in Bio.SeqIO.parse(inf, 'fasta'):
                if len(seq.seq.replace("-","").replace("N","")) >= min_length:
                    return False
        return True

    def execute(self, subcommand, args, outfile, threads=None):
        ''' generic execution of skani
        '''

        # build the skani command
        tool_cmd = [self.install_and_get_path(), subcommand]
        tool_cmd.extend(map(str, args))
        tool_cmd.extend(["-t", "{}".format(util.misc.sanitize_thread_count(threads))])

        # run the command
        _log.debug(' '.join(tool_cmd) + ' > ' + outfile)
        with open(outfile, 'w') as outf:
            util.misc.run_and_save(tool_cmd, outf=outf)

    def triangle(self, ref_fastas, outfile_ani, other_args = (), threads=None):
        ''' skani triangle computes an all-to-all ANI distance matrix for a set of sequences
        '''
        self.execute('triangle', list(ref_fastas) + list(other_args), outfile_ani, threads=threads)

    def dist(self, query_fasta, ref_fastas, outfile, other_args = (), threads=None):
        ''' skani dist computes ANI distance between a specified query set of
            sequences (MAGs) and reference genomes (database)
        '''
        if self._is_fasta_basically_empty(query_fasta):
            _log.warning("Query fasta file is empty or contains only very short sequences. Skipping skani dist (which will fail in this scenario).")
            with open(outfile, 'w') as outf:
                outf.write('\t'.join(['Ref_file','Query_file','ANI','Align_fraction_ref','Align_fraction_query','Ref_name','Query_name','Num_ref_contigs','Num_query_contigs','ANI_5_percentile','ANI_95_percentile','Standard_deviation','Ref_90_ctg_len','Ref_50_ctg_len','Ref_10_ctg_len','Query_90_ctg_len','Query_50_ctg_len','Query_10_ctg_len','Avg_chain_len','Total_bases_covered']) + '\n')
        else:
            self.execute('dist', ['-q', query_fasta, '-r'] + list(ref_fastas) + list(other_args), outfile, threads=threads)

    def find_reference_clusters(self, ref_fastas,
                                m=50, s=50, c=20, min_af=15,
                                other_args = ('--no-learned-ani', '--robust', '--detailed', '--ci', '--sparse'),
                                threads=None):
        ''' use skani triangle to define clusters of highly-related genomes
            (default settings here are for viral genomes)
        '''
        g = UndirectedGraph()
        for ref_fasta in ref_fastas:
            g.add_node(ref_fasta)

        with util.file.tempfname('.skani_matrix.ani') as tmp_matrix_ani:
            # run skani triangle
            self.triangle(ref_fastas, tmp_matrix_ani,
                          ['-m', m, '-c', c, '-s', s, '--min-af', min_af] + list(other_args), threads=threads)

            # parse the skani triangle results and define clusters
            with open(tmp_matrix_ani, 'r') as inf:
                for row in csv.DictReader(inf, delimiter='\t'):
                    g.add_edge(row['Ref_file'], row['Query_file'])

        return g.get_clusters()

    def find_closest_reference(self, contigs_fasta, ref_fastas, out_file,
                                m=50, s=50, c=20, min_af=15,
                                other_args = ('--no-learned-ani', '--robust', '--detailed', '--ci', '-n', 10, '--no-marker-index'),
                                threads=None):
        ''' use skani dist to find the closest reference genome for each contig
            (default settings here are for viral genomes)
        '''
        self.dist(contigs_fasta, ref_fastas, out_file,
                  ['-m', m, '-c', c, '-s', s, '--min-af', min_af] + list(other_args), threads=threads)
        with open(out_file, 'r') as inf:
            top_row = None
            for row in csv.DictReader(inf, delimiter='\t'):
                top_row = row
                break
        return (top_row['Ref_file'], top_row['ANI'], top_row['Align_fraction_ref'], top_row['Total_bases_covered']) if top_row is not None else None

#!/usr/bin/env python3
''' This script contains a number of utilities for metagenomic analyses.
'''
from __future__ import print_function
from __future__ import division

__author__ = "yesimon@broadinstitute.org"

import argparse
import collections
import csv
import glob
import gzip
import io
import itertools
import logging
import os.path
import operator
import queue
import re
import shutil
import sys
import tempfile
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import anndata
import numpy as np
import pandas as pd
import pysam

import util.cmd
import util.file
import util.misc
import tools.picard
import tools.samtools

import classify.kaiju
import classify.kraken
import classify.kraken2
import classify.krona
import classify.kb

__commands__ = []

log = logging.getLogger(__name__)


class TaxIdError(ValueError):
    '''Taxonomy ID couldn't be determined.'''


def maybe_compressed(fn):
    fn_gz = fn + '.gz'
    if os.path.exists(fn):
        return fn
    elif os.path.exists(fn_gz):
        return fn_gz
    else:
        raise FileNotFoundError(fn)


class TaxonomyDb(object):
    """
        This class loads NCBI taxonomy information from:
        ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
    """

    def __init__(
        self,
        tax_dir=None,
        gis=None,
        nodes=None,
        names=None,
        gis_paths=None,
        nodes_path=None,
        names_path=None,
        load_gis=False,
        load_nodes=False,
        load_names=False
    ):
        if tax_dir:
            if load_gis and not gis:
                gis_paths = [maybe_compressed(os.path.join(tax_dir, 'gi_taxid_nucl.dmp')),
                            maybe_compressed(os.path.join(tax_dir, 'gi_taxid_prot.dmp'))]
            else:
                gis_paths = []
            nodes_path = load_nodes and maybe_compressed(os.path.join(tax_dir, 'nodes.dmp')) or None
            names_path = load_names and maybe_compressed(os.path.join(tax_dir, 'names.dmp')) or None
        self.tax_dir = tax_dir
        self.gis_paths = gis_paths
        self.nodes_path = nodes_path
        self.names_path = names_path
        if load_gis:
            if gis:
                self.gis = gis
            elif gis_paths:
                self.gis = {}
                for gi_path in gis_paths:
                    log.info('Loading taxonomy gis: %s', gi_path)
                    self.gis.update(self.load_gi_single_dmp(gi_path))
        if load_nodes:
            if nodes:
                self.ranks, self.parents = nodes
            elif nodes_path:
                log.info('Loading taxonomy nodes: %s', nodes_path)
                self.ranks, self.parents = self.load_nodes(nodes_path)
        if load_names:
            if names:
                self.names = names
            elif names_path:
                log.info('Loading taxonomy names: %s', names_path)
                self.names = self.load_names(names_path)

    def load_gi_single_dmp(self, dmp_path):
        '''Load a gi->taxid dmp file from NCBI taxonomy.'''
        gi_array = {}
        with util.file.open_or_gzopen(dmp_path) as f:
            for i, line in enumerate(f):
                gi, taxid = line.strip().split('\t')
                gi = int(gi)
                taxid = int(taxid)
                gi_array[gi] = taxid
                if (i + 1) % 1000000 == 0:
                    log.info('Loaded %s gis', i)
        return gi_array

    def load_names(self, names_db, scientific_only=True):
        '''Load the names.dmp file from NCBI taxonomy.'''
        if scientific_only:
            names = {}
        else:
            names = collections.defaultdict(list)
        for line in util.file.open_or_gzopen(names_db):
            parts = line.strip().split('|')
            taxid = int(parts[0])
            name = parts[1].strip()
            #unique_name = parts[2].strip()
            class_ = parts[3].strip()
            if scientific_only:
                if class_ == 'scientific name':
                    names[taxid] = name
            else:
                names[taxid].append(name)
        return names

    def load_nodes(self, nodes_db):
        '''Load ranks and parents arrays from NCBI taxonomy.'''
        ranks = {}
        parents = {}
        with util.file.open_or_gzopen(nodes_db) as f:
            for line in f:
                parts = line.strip().split('|')
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
                rank = parts[2].strip()
                #embl_code = parts[3].strip()
                #division_id = parts[4].strip()
                parents[taxid] = parent_taxid
                ranks[taxid] = rank
        return ranks, parents

    def get_ordered_ancestors(self, taxid):
        ''' returns all ancestors of a taxid in proximity order: [parent, grandparent, greatgrandparent, etc] '''
        if taxid in self.parents and taxid != self.parents[taxid]:
            return [self.parents[taxid]] + self.get_ordered_ancestors(self.parents[taxid])
        else:
            return []

    def process_blast_hits(self, hits, top_percent):
        '''Filter groups of blast hits and perform lca.

        Args:
        hits: []BlastRecord groups of hits.
        top_percent: (float) Only consider hits within this percent of top bit score.

        Return:
        (int) Tax id of LCA.
        '''
        hits = (self.translate_gi_to_tax_id(hit) for hit in hits)

        hits = [hit for hit in hits if hit.subject_id != 0]
        if len(hits) == 0:
            return

        best_score = max(hit.bit_score for hit in hits)
        cutoff_bit_score = (100 - top_percent) / 100 * best_score
        valid_hits = (hit for hit in hits if hit.bit_score >= cutoff_bit_score)
        valid_hits = list(valid_hits)
        # Sort requires realized list
        valid_hits.sort(key=operator.attrgetter('bit_score'), reverse=True)
        if valid_hits:
            tax_ids = tuple(itertools.chain(*(blast_m8_taxids(hit) for hit in valid_hits)))
            return coverage_lca(tax_ids, self.parents)

    def process_sam_hits(self, sam_hits, top_percent):
        '''Filter groups of blast hits and perform lca.

        Args:
        sam_hits: []Sam groups of hits.
        top_percent: (float) Only consider hits within this percent of top bit score.

        Return:
        (int) Tax id of LCA.
        '''
        best_score = max(hit.get_tag('AS') for hit in sam_hits)
        cutoff_alignment_score = (100 - top_percent) / 100 * best_score
        valid_hits = (hit for hit in sam_hits if hit.get_tag('AS') >= cutoff_alignment_score)
        valid_hits = list(valid_hits)
        # Sort requires realized list
        valid_hits.sort(key=lambda sam1: sam1.get_tag('AS'), reverse=True)

        tax_ids = [extract_tax_id(hit) for hit in valid_hits]
        return coverage_lca(tax_ids, self.parents)

    def translate_gi_to_tax_id(self, record):
        '''Replace gi headers in subject ids to int taxonomy ids.'''
        gi = int(record.subject_id.split('|')[1])
        tax_id = self.gis[gi]
        rec_list = list(record)
        rec_list[1] = tax_id
        return BlastRecord(*rec_list)

    def sam_lca(self, sam_file, output=None, top_percent=10, unique_only=True):
        ''' Calculate the LCA taxonomy id for multi-mapped reads in a samfile.

        Assumes the sam is sorted by query name. Writes tsv output: query_id \t tax_id.

        Args:
        sam_file: (path) Sam file.
        output: (io) Output file.
        top_percent: (float) Only this percent within top hit are used.
        unique_only: (bool) If true, only output assignments for unique, mapped reads. If False, set unmapped or duplicate reads as unclassified.

        Return:
        (collections.Counter) Counter of taxid hits
        '''

        c = collections.Counter()
        with pysam.AlignmentFile(sam_file, 'rb') as sam:
            seg_groups = (v for k, v in itertools.groupby(sam, operator.attrgetter('query_name')))
            for seg_group in seg_groups:
                segs = list(seg_group)
                query_name = segs[0].query_name
                # 0x4 is unmapped, 0x400 is duplicate
                mapped_segs = [seg for seg in segs if seg.flag & 0x4 == 0 and seg.flag & 0x400 == 0]
                if unique_only and not mapped_segs:
                    continue

                if mapped_segs:
                    tax_id = self.process_sam_hits(mapped_segs, top_percent)
                    if tax_id is None:
                        log.warning('Query: {} has no valid taxonomy paths.'.format(query_name))
                        if unique_only:
                            continue
                        else:
                            tax_id = 0
                else:
                    tax_id = 0

                if output:
                    classified = 'C' if tax_id else 'U'
                    output.write('{}\t{}\t{}\n'.format(classified, query_name, tax_id))
                c[tax_id] += 1
        return c

    def sam_lca_report(self, bam_aligned, outReport, outReads=None, unique_only=None):
        if outReads:
            lca_tsv = outReads
        else:
            lca_tsv = util.file.mkstempfname('.tsv')
        with util.file.open_or_gzopen(lca_tsv, 'wt') as lca:
            hits = self.sam_lca(bam_aligned, lca, top_percent=10, unique_only=unique_only)
        with open(outReport, 'w') as f:
            for line in self.kraken_dfs_report(hits):
                print(line, file=f)

    def blast_lca(self,
                m8_file,
                output,
                paired=False,
                min_bit_score=50,
                max_expected_value=0.01,
                top_percent=10,):
        '''Calculate the LCA taxonomy id for groups of blast hits.

        Writes tsv output: query_id \t tax_id

        Args:
        m8_file: (io) Blast m8 file to read.
        output: (io) Output file.
        paired: (bool) Whether to count paired suffixes /1,/2 as one group.
        min_bit_score: (float) Minimum bit score or discard.
        max_expected_value: (float) Maximum e-val or discard.
        top_percent: (float) Only this percent within top hit are used.
        '''
        records = blast_records(m8_file)
        records = (r for r in records if r.e_val <= max_expected_value)
        records = (r for r in records if r.bit_score >= min_bit_score)
        if paired:
            records = (paired_query_id(rec) for rec in records)
        blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
        for blast_group in blast_groups:
            blast_group = list(blast_group)
            tax_id = self.process_blast_hits(blast_group, top_percent)
            query_id = blast_group[0].query_id
            if not tax_id:
                log.debug('Query: {} has no valid taxonomy paths.'.format(query_id))
            classified = 'C' if tax_id else 'U'
            output.write('{}\t{}\t{}\n'.format(classified, query_id, tax_id))

    def kraken_dfs(self, lines, taxa_hits, total_hits, taxid, level):
        '''Recursively do DFS for number of hits per taxa.'''
        cum_hits = num_hits = taxa_hits.get(taxid, 0)
        for child_taxid in self.children[taxid]:
            cum_hits += self.kraken_dfs(lines, taxa_hits, total_hits, child_taxid, level + 1)
        percent_covered = '%.2f' % (cum_hits / total_hits * 100)
        rank = rank_code(self.ranks[taxid])
        name = self.names[taxid]
        if cum_hits > 0:
            lines.append('\t'.join([percent_covered, str(cum_hits), str(num_hits), rank, str(taxid), '  ' * level + name]))
        return cum_hits

    def kraken_dfs_report(self, taxa_hits):
        '''Return a kraken compatible DFS report of taxa hits.

        Warning: this potentially mutates the self.children variable

        Args:
        taxa_hits: (collections.Counter) # of hits per tax id.

        Return:
        []str lines of the report
        '''

        self.children = parents_to_children(self.parents)  # huh??? are we sure we want to mutate?
        total_hits = sum(taxa_hits.values())
        if total_hits == 0:
            return ['\t'.join(['100.00', '0', '0', 'U', '0', 'unclassified'])]

        lines = []
        self.kraken_dfs(lines, taxa_hits, total_hits, 1, 0)
        unclassified_hits = taxa_hits.get(0, 0)
        unclassified_hits += taxa_hits.get(-1, 0)

        if unclassified_hits > 0:
            percent_covered = '%.2f' % (unclassified_hits / total_hits * 100)
            lines.append(
                '\t'.join([
                    str(percent_covered), str(unclassified_hits), str(unclassified_hits), 'U', '0', 'unclassified'
                ])
            )
        return reversed(lines)



BlastRecord = collections.namedtuple(
    'BlastRecord', [
        'query_id', 'subject_id', 'percent_identity', 'aln_length', 'mismatch_count', 'gap_open_count', 'query_start',
        'query_end', 'subject_start', 'subject_end', 'e_val', 'bit_score', 'extra'
    ]
)

def blast_records(f):
    '''Yield blast m8 records line by line'''
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        for field in range(3, 10):
            parts[field] = int(parts[field])
        for field in (2, 10, 11):
            parts[field] = float(parts[field])
        args = parts[:12]
        extra = parts[12:]
        args.append(extra)

        yield BlastRecord(*args)

def paired_query_id(record):
    '''Replace paired suffixes in query ids.'''
    suffixes = ('/1', '/2')
    for suffix in suffixes:
        if record.query_id.endswith(suffix):
            rec_list = list(record)
            rec_list[0] = record.query_id[:-len(suffix)]
            return BlastRecord(*rec_list)
    return record

def blast_m8_taxids(record):
    return [int(record.subject_id)]

def extract_tax_id(sam1):
    '''Replace gi headers in subject ids to int taxonomy ids.'''
    parts = sam1.reference_name.split('|')
    if parts[0] == 'taxid':
        return int(parts[1])
    else:
        raise TaxIdError(parts)



def coverage_lca(query_ids, parents, lca_percent=100):
    '''Calculate the lca that will cover at least this percent of queries.

    Args:
      query_ids: []int list of nodes.
      parents: []int array of parents.
      lca_percent: (float) Cover at least this percent of queries.

    Return:
      (int) LCA
    '''
    lca_needed = lca_percent / 100 * len(query_ids)
    paths = []
    for query_id in query_ids:
        path = []
        while query_id != 1:
            path.append(query_id)
            if parents.get(query_id, 0) == 0:
                log.warning('Parent for query id: {} missing'.format(query_id))
                break
            query_id = parents[query_id]
        if query_id == 1:
            path.append(1)
            path = list(reversed(path))
            paths.append(path)
    if not paths:
        return

    last_common = 1
    max_path_length = max(len(path) for path in paths)
    for level in range(max_path_length):
        valid_paths = (path for path in paths if len(path) > level)
        max_query_id, hits_covered = collections.Counter(path[level] for path in valid_paths).most_common(1)[0]
        if hits_covered >= lca_needed:
            last_common = max_query_id
        else:
            break
    return last_common


def tree_level_lookup(parents, node, level_cache):
    '''Get the node level/depth.

    Args:
      parents: Array of node parents.
      node: Node to get level (root == 1).
      level_cache: Cache of previously found levels.

    Returns:
      (int) level of node
    '''
    path = []
    while True:
        level = level_cache.get(node)
        if level:
            for i, node in enumerate(reversed(path)):
                level_cache[node] = level + i + 1
            return level + len(path)
        path.append(node)
        node = parents[node]


def push_up_tree_hits(parents, hits, min_support_percent=None, min_support=None, update_assignments=False):
    '''Push up hits on nodes until min support is reached.

    Args:
      parents: Array of node parents.
      hits: Counter of hits on each node.
      min_support_percent: Push up hits until each node has
        this percent of the sum of all hits.
      min_support: Push up hits until each node has this number of hits.

    Returns:
      (counter) Hits mutated pushed up the tree.
    '''
    assert min_support_percent or min_support

    if update_assignments:
        pass

    total_hits = sum(hits.values())
    if not min_support:
        min_support = round(min_support_percent * 0.01 * total_hits)
    pq_level = queue.PriorityQueue()
    level_cache = {1: 1}
    for hit_id, num_hits in hits.items():
        if num_hits < min_support:
            pq_level.put((-tree_level_lookup(parents, hit_id, level_cache), hit_id))

    while not pq_level.empty() > 0:
        level, hit_id = pq_level.get()
        level = -level

        if hits[hit_id] >= min_support:
            continue
        if hit_id == 1:
            del hits[1]
            break

        parent_hit_id = parents[hit_id]
        num_hits = hits[hit_id]
        hits[parent_hit_id] += num_hits
        # Can't pop directly from hits because hit_id might not be stored in counter
        if hit_id in hits:
            del hits[hit_id]
        if hits[parent_hit_id] < min_support:
            pq_level.put((-tree_level_lookup(parents, parent_hit_id, level_cache), parent_hit_id))
    return hits


def parents_to_children(parents):
    '''Convert an array of parents to lists of children for each parent.

    Returns:
      (dict[list]) Lists of children
    '''
    children = collections.defaultdict(list)
    for node, parent in parents.items():
        if node == 1:
            continue
        if parent != 0:
            children[parent].append(node)
    return children



def file_lines(filename):
    if filename is not None:
        with open(filename) as f:
            for line in f:
                yield line


def collect_children(children, original_taxids):
    '''Collect nodes with all children recursively.'''
    taxids = original_taxids
    while taxids:
        taxid = taxids.pop()
        yield taxid
        for child_taxid in children[taxid]:
            taxids.add(child_taxid)


def collect_parents(parents, taxids):
    '''Collect nodes with all parents recursively.'''
    # The root taxid node is 1
    yield 1
    taxids_with_parents = set([1])
    for taxid in taxids:
        while taxid not in taxids_with_parents:
            yield taxid
            taxids_with_parents.add(taxid)
            taxid = parents[taxid]


def parser_subset_taxonomy(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "taxDb",
        help="Taxonomy database directory (containing nodes.dmp, parents.dmp etc.)",
    )
    parser.add_argument(
        "outputDb",
        help="Output taxonomy database directory",
    )
    parser.add_argument(
        "--whitelistTaxids",
        help="List of taxids to add to taxonomy (with parents)",
        nargs='+', type=int
    )
    parser.add_argument(
        "--whitelistTaxidFile",
        help="File containing taxids - one per line - to add to taxonomy with parents.",
    )
    parser.add_argument(
        "--whitelistTreeTaxids",
        help="List of taxids to add to taxonomy (with parents and children)",
        nargs='+', type=int
    )
    parser.add_argument(
        "--whitelistTreeTaxidFile",
        help="File containing taxids - one per line - to add to taxonomy with parents and children.",
    )
    parser.add_argument(
        "--whitelistGiFile",
        help="File containing GIs - one per line - to add to taxonomy with nodes.",
    )
    parser.add_argument(
        "--whitelistAccessionFile",
        help="File containing accessions - one per line - to add to taxonomy with nodes.",
    )
    parser.add_argument(
        "--skipGi", action='store_true',
        help="Skip GI to taxid mapping files"
    )
    parser.add_argument(
        "--skipAccession", action='store_true',
        help="Skip accession to taxid mapping files"
    )
    parser.add_argument(
        "--skipDeadAccession", action='store_true',
        help="Skip dead accession to taxid mapping files"
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, subset_taxonomy, split_args=True)
    return parser
def subset_taxonomy(taxDb, outputDb, whitelistTaxids=None, whitelistTaxidFile=None,
                    whitelistTreeTaxids=None, whitelistTreeTaxidFile=None,
                    whitelistGiFile=None, whitelistAccessionFile=None,
                    skipGi=None, skipAccession=None, skipDeadAccession=None,
                    stripVersion=True):
    '''
    Generate a subset of the taxonomy db files filtered by the whitelist. The
    whitelist taxids indicate specific taxids plus their parents to add to
    taxonomy while whitelistTreeTaxids indicate specific taxids plus both
    parents and all children taxa. Whitelist GI and accessions can only be
    provided in file form and the resulting gi/accession2taxid files will be
    filtered to only include those in the whitelist files. Finally, taxids +
    parents for the gis/accessions will also be included.
    '''
    util.file.mkdir_p(os.path.join(outputDb, 'accession2taxid'))
    db = TaxonomyDb(tax_dir=taxDb, load_nodes=True)

    taxids = set()
    if whitelistTaxids is not None:
        taxids.update(set(whitelistTaxids))
    taxids.update((int(x) for x in file_lines(whitelistTaxidFile)))

    tree_taxids = set()
    if whitelistTreeTaxids is not None:
        tree_taxids.update(set(whitelistTreeTaxids))
    taxids.update((int(x) for x in file_lines(whitelistTreeTaxidFile)))
    keep_taxids = set(collect_parents(db.parents, taxids))

    if tree_taxids:
        db.children = parents_to_children(db.parents)
        children_taxids = collect_children(db.children, tree_taxids)
        keep_taxids.update(children_taxids)

    # Taxids kept based on GI or Accession. Get parents afterwards to not pull in all GIs/accessions.
    keep_seq_taxids = set()
    def filter_file(path, sep='\t', taxid_column=0, gi_column=None, a2t=False, header=False):
        input_path = os.path.join(db.tax_dir, path)
        output_path = os.path.join(outputDb, path)

        input_path = maybe_compressed(input_path)
        with util.file.open_or_gzopen(input_path, 'rt') as f, \
             util.file.open_or_gzopen(output_path, 'wt') as out_f:
            if header:
                out_f.write(f.readline())  # Cannot use next(f) for python2
            for line in f:
                parts = line.split(sep)
                taxid = int(parts[taxid_column])
                if gi_column is not None:
                    gi = int(parts[gi_column])
                    if gi in gis:
                        keep_seq_taxids.add(taxid)
                        out_f.write(line)
                        continue
                if a2t:
                    accession = parts[accession_column_i]
                    if stripVersion:
                        accession = accession.split('.', 1)[0]
                    if accession in accessions:
                        keep_seq_taxids.add(taxid)
                        out_f.write(line)
                        continue
                if taxid in keep_taxids:
                    out_f.write(line)

    if not skipGi:
        gis = set(int(x) for x in file_lines(whitelistGiFile))

        filter_file('gi_taxid_nucl.dmp', taxid_column=1, gi_column=0)
        filter_file('gi_taxid_prot.dmp', taxid_column=1, gi_column=0)

    if not skipAccession:
        if stripVersion:
            accessions = set(x.strip().split('.', 1)[0] for x in file_lines(whitelistAccessionFile))
            accession_column_i = 0
        else:
            accessions = set(file_lines(whitelistAccessionFile))
            accession_column_i = 1

        acc_dir = os.path.join(db.tax_dir, 'accession2taxid')
        acc_paths = []
        for fn in os.listdir(acc_dir):
            if fn.endswith('.accession2taxid') or fn.endswith('.accession2taxid.gz'):
                if skipDeadAccession and fn.startswith('dead_'):
                    continue
                acc_paths.append(os.path.join(acc_dir, fn))
        for acc_path in acc_paths:
            filter_file(os.path.relpath(acc_path, db.tax_dir), taxid_column=2, header=True, a2t=True)


    # Add in taxids found from processing GI/accession
    keep_seq_taxids = collect_parents(db.parents, keep_seq_taxids)
    keep_taxids.update(keep_seq_taxids)

    filter_file('nodes.dmp', sep='|')
    filter_file('names.dmp', sep='|')
    filter_file('merged.dmp')
    filter_file('delnodes.dmp')
__commands__.append(('subset_taxonomy', parser_subset_taxonomy))


def parser_filter_taxids_to_focal_hits(parser=argparse.ArgumentParser()):
    parser.add_argument("taxids_tsv",   help="TSV file where first column is a taxid")
    parser.add_argument("focal_report_tsv", help="TSV produced by taxlevel_plurality")
    parser.add_argument("taxdb_dir", help="Taxonomy database directory (containing nodes.dmp, parents.dmp etc.)")
    parser.add_argument("min_read_count", type=int, help="ignore focal_report_tsv entries below this read count")
    parser.add_argument("output_tsv",  help="Output TSV file where first column is a taxid")

    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_taxids_to_focal_hits, split_args=True)
    return parser
def filter_taxids_to_focal_hits(taxids_tsv, focal_report_tsv, taxdb_dir, min_read_count, output_tsv):
    '''
    Generate a subset of the taxids_tsv file filtered by the focal_report_tsv.
    We will only emit rows from the taxids_tsv that contain taxids that are either
    contained within or are a child/descendant of nodes contained within the
    focal_report_tsv
    '''

    # load taxonomy database structure
    taxdb = TaxonomyDb(tax_dir=taxdb_dir, load_nodes=True, load_gis=False)

    # load focal hits
    hits = set()
    with util.file.open_or_gzopen(focal_report_tsv, "rt") as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
            if int(row['reads_excl_children']) >= min_read_count:
                hits.add(int(row['taxon_id']))

    # filter taxids_tsv -> output_tsv
    with util.file.open_or_gzopen(taxids_tsv, "rt") as inf:
        with util.file.open_or_gzopen(output_tsv, "wt") as outf:
            for line in inf:
                taxid = int(line.rstrip('\r\n').split('\t')[0])
                ancestors = taxdb.get_ordered_ancestors(taxid)
                if any(node in hits for node in [taxid] + ancestors):
                    outf.write(line)

__commands__.append(('filter_taxids_to_focal_hits', parser_filter_taxids_to_focal_hits))


def rank_code(rank):
    '''Get the short 1 letter rank code for named ranks.'''
    if rank == "species":
        return "S"
    elif rank == "genus":
        return "G"
    elif rank == "family":
        return "F"
    elif rank == "order":
        return "O"
    elif rank == "class":
        return "C"
    elif rank == "phylum":
        return "P"
    elif rank == "kingdom":
        return "K"
    elif rank == "superkingdom":
        return "D"
    elif rank == "unclassified":
        return "U"
    else:
        return "-"


def taxa_hits_from_tsv(f, taxid_column=2):
    '''Return a counter of hits from tsv.'''
    c = collections.Counter()
    for row in csv.reader(f, delimiter='\t'):
        tax_id = int(row[taxid_column - 1])
        c[tax_id] += 1
    return c


def parser_kraken2(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('inBams', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outReports', nargs='+', help='Kraken2 summary report output file. Multiple filenames space separated.')
    parser.add_argument('--outReads', nargs='+', help='Kraken2 per read classification output file. Multiple filenames space separated.')
    parser.add_argument(
        '--minimum_hit_groups', default=None, type=int, help='Minimum hit groups (Kraken2 default: 2)'
    )
    parser.add_argument(
        '--min_base_qual', default=None, type=int, help='Minimum base quality (default %(default)s)'
    )
    parser.add_argument(
        '--confidence', default=None, type=float, help='Kraken2 confidence score threshold (default %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken2, split_args=True)
    return parser
def kraken2(db, inBams, outReports=None, outReads=None, min_base_qual=None, confidence=None, minimum_hit_groups=None, threads=None):
    '''
        Classify reads by taxon using Kraken2
    '''

    assert outReads or outReports, ('Either --outReads or --outReport must be specified.')
    kraken_tool = classify.kraken2.Kraken2()
    kraken_tool.pipeline(db, inBams, out_reports=outReports, out_reads=outReads,
                         min_base_qual=min_base_qual, confidence=confidence,
                         minimum_hit_groups=minimum_hit_groups, num_threads=threads)
__commands__.append(('kraken2', parser_kraken2))


def parser_kb(parser=argparse.ArgumentParser()):
    """Argument parser for the kb_python wrapper.

    Args:
        parser (_type_, optional): _description_. Defaults to argparse.ArgumentParser().

    Returns:
        argparse.ArgumentParser: The parser with arguments added.
    """
    parser.add_argument('inBam', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--index', help='kb index file.')
    parser.add_argument('--t2g', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--kmerSize', type=int, help='k-mer size (default: 31bp)', default=31)
    parser.add_argument('--technology', choices=['10xv2', '10xv3', '10xv3-3prime', '10xv3-5prime', 'dropseq', 
                                                 'indrop', 'celseq', 'celseq2', 'smartseq2', 'bulk'], 
                        help='Technology used to generate the data (default: bulk)', default='bulk')
    parser.add_argument('--h5ad', action='store_true', help='Output HDF5 file (default: False)', default=False)
    parser.add_argument('--loom', action='store_true', help='Output Loom file (default: False)', default=False)
    parser.add_argument('--outDir', help='Output directory (default: kb_out)', default='kb_out')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken2, split_args=True)
    return parser
def kb_python(inBams, index_file, t2g_file, kmer_size=31, technology='bulk', h5ad=False, loom=False, outDir=None, threads=None):
    """Runs kb count on the input BAM files.

    Args:
        inBams (list): List of input BAM files.
        outDir (str): Output directory. Defaults to None.
        index_file (str): Path to the kb index file.
        t2g_file (str): Path to the transcript-to-gene mapping file.
        kmer_size (int, optional): K-mer size for the alignment. Defaults to 31.
        technology (str, optional): Sequencing technology used. Defaults to 'bulk'.
        h5ad (bool, optional): Whether to output HDF5 file. Defaults to False.
        loom (bool, optional): Whether to output Loom file. Defaults to False.
        threads (int, optional): Number of threads to use. Defaults to None.
    """

    assert outDir, ('Output directory must be specified.')
    kb_tool = classify.kb.kb()
    kb_tool.pipeline(
        inBams=inBams,
        outDir=outDir,
        index_file=index_file,
        t2g_file=t2g_file,
        kmer_size=kmer_size,
        technology=technology,
        h5ad=h5ad,
        loom=loom,
        threads=threads
    )
__commands__.append(('kb', parser_kb))

def parser_krakenuniq(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('inBams', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outReports', nargs='+', help='Kraken summary report output file. Multiple filenames space separated.')
    parser.add_argument('--outReads', nargs='+', help='Kraken per read classification output file. Multiple filenames space separated.')
    parser.add_argument(
        '--filterThreshold', default=0.05, type=float, help='Kraken filter threshold (default %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, krakenuniq, split_args=True)
    return parser
def krakenuniq(db, inBams, outReports=None, outReads=None, lockMemory=False, filterThreshold=None, threads=None):
    '''
        Classify reads by taxon using KrakenUniq
    '''

    assert outReads or outReports, ('Either --outReads or --outReport must be specified.')
    kuniq_tool = classify.kraken.KrakenUniq()
    kuniq_tool.pipeline(db, inBams, out_reports=outReports, out_reads=outReads,
                        filter_threshold=filterThreshold, num_threads=threads)
# Disable KrakenUniq from future versions for now, pending conda rebuild of @yesimon's fork and its dependencies
#__commands__.append(('krakenuniq', parser_krakenuniq))


def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inReports', nargs='+', help='Input report file (default: tsv)')
    parser.add_argument('db', help='Krona taxonomy database directory.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--sample_name', help='Title of dataset (default basename(inReport))', default=None)
    parser.add_argument('--queryColumn', help='Column of query id. (default %(default)s)', type=int, default=2)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id. (default %(default)s)', type=int, default=3)
    parser.add_argument('--scoreColumn', help='Column of score. (default %(default)s)', type=int, default=None)
    parser.add_argument('--magnitudeColumn', help='Column of magnitude. (default %(default)s)', type=int, default=None)
    parser.add_argument('--noHits', help='Include wedge for no hits.', action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.', action='store_true')
    parser.add_argument('--inputType', help='Handling for specialized report types.', default='tsv', choices=['tsv', 'kraken2', 'krakenuniq', 'kaiju'])
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser
def krona(inReports, db, outHtml, queryColumn=None, taxidColumn=None, scoreColumn=None, magnitudeColumn=None, noHits=None, noRank=None,
          inputType=None, sample_name=None):
    '''
        Create an interactive HTML report from a tabular metagenomic report
    '''

    krona_tool = classify.krona.Krona()
    if sample_name is not None:
        dataset_names = list(sample_name for fn in inReports)
    else:
        dataset_names = list(os.path.basename(fn) for fn in inReports)

    with util.file.tmp_dir() as tmp_dir:
        to_import = []

        if inputType == 'tsv':
            for inReport in inReports:
                if inReport.endswith('.gz'):
                    tmp_tsv = util.file.mkstempfname('.tsv', directory=tmp_dir)
                    with gzip.open(inReport, 'rb') as f_in:
                        with open(tmp_tsv, 'w') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                            to_import.append(tmp_tsv)
                else:
                    to_import.append(inReport)

        elif inputType == 'kraken2':
            queryColumn=None
            taxidColumn=5
            scoreColumn=None
            magnitudeColumn=3
            to_import = inReports

        elif inputType == 'krakenuniq':
            queryColumn=None
            taxidColumn=1
            scoreColumn=3
            magnitudeColumn=2
            for inReport in inReports:
                tmp_tsv = util.file.mkstempfname('.tsv', directory=tmp_dir)
                with open(tmp_tsv, 'w') as f_out:
                    report = classify.kraken.KrakenUniq().read_report(inReport)
                    for taxid, (tax_reads, tax_kmers) in report.items():
                        f_out.write('{}\t{}\t{}\n'.format(taxid, tax_reads, tax_kmers))
                to_import.append(tmp_tsv)

        elif inputType == 'kaiju':
            queryColumn=None
            taxidColumn=1
            scoreColumn=None
            magnitudeColumn=2
            for inReport in inReports:
                tmp_tsv = util.file.mkstempfname('.tsv', directory=tmp_dir)
                with open(tmp_tsv, 'w') as f_out:
                    report = classify.kaiju.Kaiju().read_report(inReport)
                    for taxid, reads in report.items():
                        f_out.write('{}\t{}\n'.format(taxid, reads))
                to_import.append(tmp_tsv)

        else:
            raise NotImplementedError

        # run krona
        krona_tool.import_taxonomy(
            db,
            list(','.join(pair) for pair in zip(to_import, dataset_names)),
            outHtml,
            query_column=queryColumn,
            taxid_column=taxidColumn,
            score_column=scoreColumn,
            magnitude_column=magnitudeColumn,
            no_hits=noHits,
            no_rank=noRank
        )

        if inputType == 'krakenuniq':
            # Rename "Avg. score" to "Est. genome coverage"
            html_lines = util.file.slurp_file(outHtml).split('\n')
            tmp_tsv = util.file.mkstempfname('.tsv', directory=tmp_dir)
            with open(tmp_tsv, 'w') as new_report:
                for line in html_lines:
                    if '<attribute display="Avg. score">score</attribute>' in line:
                        line = line.replace('Avg. score', 'Est. unique kmers')
                    print(line, file=new_report)
            shutil.copyfile(tmp_tsv, outHtml)

__commands__.append(('krona', parser_krona))


def parser_kaiju(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Kaiju database .fmi file.')
    parser.add_argument('taxDb', help='Taxonomy database directory.')
    parser.add_argument('outReport', help='Output taxonomy report.')
    parser.add_argument('--outReads', help='Output LCA assignments for each read.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kaiju, split_args=True)
    return parser
def kaiju(inBam, db, taxDb, outReport, outReads=None, threads=None):
    '''
        Classify reads by the taxon of the Lowest Common Ancestor (LCA)
    '''

    kaiju_tool = classify.kaiju.Kaiju()
    kaiju_tool.classify(db, taxDb, inBam, output_report=outReport, output_reads=outReads, num_threads=threads)
__commands__.append(('kaiju', parser_kaiju))


def parser_metagenomic_report_merge(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "metagenomic_reports",
        help="Input metagenomic reports with the query ID and taxon ID in the 2nd and 3rd columns (Kraken format)",
        nargs='+',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        "--outSummaryReport",
        dest="out_kraken_summary",
        help="Path of human-readable metagenomic summary report, created by kraken-report"
    )
    parser.add_argument(
        "--krakenDB",
        dest="kraken_db",
        help="Kraken database (needed for outSummaryReport)",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        "--outByQueryToTaxonID", dest="out_krona_input", help="Output metagenomic report suitable for Krona input. "
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, metagenomic_report_merge, split_args=True)
    return parser
def metagenomic_report_merge(metagenomic_reports, out_kraken_summary, kraken_db, out_krona_input):
    '''
        Merge multiple metegenomic reports into a single metagenomic report.
        Any Krona input files created by this
    '''
    assert out_kraken_summary or out_krona_input, (
        "Either --outSummaryReport or --outByQueryToTaxonID must be specified"
    )
    assert kraken_db if out_kraken_summary else True, (
        'A Kraken db must be provided via --krakenDB if outSummaryReport is specified'
    )

    # column numbers containing the query (sequence) ID and taxonomic ID
    # these are one-indexed
    # See: http://ccb.jhu.edu/software/kraken/MANUAL.html#output-format
    # tool_data_columns = {
    #     "kraken": (2, 3)
    # }

    # if we're creating a Krona input file
    if out_krona_input:
        # open the output file (as gz if necessary)
        with util.file.open_or_gzopen(out_krona_input, "wt") as outf:
            # create a TSV writer for the output file
            output_writer = csv.writer(outf, delimiter='\t', lineterminator='\n')

            if metagenomic_reports:
                # for each Kraken-format metag file specified, pull out the appropriate columns
                # and write them to the TSV output
                for metag_file in metagenomic_reports:
                    with util.file.open_or_gzopen(metag_file.name, "rt") as inf:
                        file_reader = csv.reader(inf, delimiter='\t')
                        for row in file_reader:
                            # for only the two relevant columns
                            output_writer.writerow([f for f in row])

    # create a human-readable summary of the Kraken reports
    # kraken-report can only be used on kraken reports since it depends on queries being in its database
    if out_kraken_summary:
        # create temporary file to hold combined kraken report
        tmp_metag_combined_txt = util.file.mkstempfname('.txt')

        util.file.cat(tmp_metag_combined_txt, [metag_file.name for metag_file in metagenomic_reports])

        kraken_tool = classify.kraken.Kraken()
        kraken_tool.report(tmp_metag_combined_txt, kraken_db.name, out_kraken_summary)
__commands__.append(('report_merge', parser_metagenomic_report_merge))



def fasta_library_accessions(library):
    '''Parse accession from ids of fasta files in library directory. '''
    library_accessions = set()
    for dirpath, dirnames, filenames in os.walk(library, followlinks=True):
        for filename in filenames:
            if not filename.endswith('.fna') and not filename.endswith('.fa') and not filename.endswith('.ffn'):
                continue
            filepath = os.path.join(dirpath, filename)
            for seqr in SeqIO.parse(filepath, 'fasta'):
                name = seqr.name
                # Search for accession
                mo = re.search(r'([A-Z]+_?\d+\.\d+)', name)
                if mo:
                    accession = mo.group(1)
                    library_accessions.add(accession)
    return library_accessions


class KrakenBuildError(Exception):
    '''Error while building kraken database.'''

class KrakenUniqBuildError(KrakenBuildError):
    '''Error while building KrakenUniq database.'''

def parser_filter_bam_to_taxa(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input bam file.')
    parser.add_argument('read_IDs_to_tax_IDs', help='TSV file mapping read IDs to taxIDs, Kraken-format by default. Assumes bijective mapping of read ID to tax ID.')
    parser.add_argument('out_bam', help='Output bam file, filtered to the taxa specified')
    parser.add_argument('nodes_dmp', help='nodes.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('names_dmp', help='names.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('--exclude', action='store_true',  dest="exclude", help='Switch filtration to remove all reads falling under matching taxa (and keep all non-matching). Default is the inverse: keep all reads falling under matching taxa (and remove all non-matching).')
    parser.add_argument('--taxNames', nargs="+", dest="tax_names", help='The taxonomic names to include. More than one can be specified. Mapped to Tax IDs by lowercase exact match only. Ex. "Viruses" This is in addition to any taxonomic IDs provided.')
    parser.add_argument('--taxIDs', nargs="+", type=int, dest="tax_ids", help='The NCBI taxonomy IDs to include. More than one can be specified. This is in addition to any taxonomic names provided.')
    parser.add_argument('--without-children', action='store_true', dest="omit_children", help='Omit reads classified more specifically than each taxon specified (without this a taxon and its children are included).')
    parser.add_argument('--read_id_col', type=int, dest="read_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing read IDs. (default: %(default)s)', default=1)
    parser.add_argument('--tax_id_col', type=int, dest="tax_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing Taxonomy IDs. (default: %(default)s)', default=2)
    parser.add_argument('--out_count', help='Write a file with the number of reads matching the specified taxa.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_bam_to_taxa, split_args=True)
    return parser

def filter_bam_to_taxa(in_bam, read_IDs_to_tax_IDs, out_bam,
                       nodes_dmp, names_dmp,
                       exclude=False,
                       tax_names=None, tax_ids=None,
                       omit_children=False,
                       read_id_col=1, tax_id_col=2,
                       out_count=None,
                       JVMmemory=None):
    """
        Filter an (already classified) input bam file to only include reads that have been mapped to specified
        taxonomic IDs or scientific names. This requires a classification file, as produced
        by tools such as Kraken, as well as the NCBI taxonomy database.
    """
    tax_ids = set(tax_ids) if tax_ids else set()
    tax_names = tax_names or []
    # use TaxonomyDb() class above and tree traversal/collection functions above
    db = TaxonomyDb(nodes_path=nodes_dmp, names_path=names_dmp, load_nodes=True, load_names=True)
    db.children = parents_to_children(db.parents)

    paired_read_base_pattern = re.compile(r'^(.*?)(/[1-2])?$')

    # get taxIDs for each of the heading values specifed (exact matches only)
    tax_ids_from_headings = set()
    for heading in tax_names:
        # map heading to taxID
        name_pattern = re.compile(heading, flags=re.IGNORECASE)
        for row_tax_id, names in db.names.items():
            found_heading = False
            if type(names) != list:
                # if taxID->str, create list (of names) with cardinality=1
                names = [names]
            for name in names:
                if name_pattern.match(name):
                    tax_ids_from_headings.add(row_tax_id)
                    log.debug("Found taxName match: %s -> %s" % (row_tax_id,name))
                    found_heading = True
                    break
            if found_heading:
                break

    tax_ids |= tax_ids_from_headings

    log.debug("tax_ids %s", tax_ids)
    log.debug("tax_names %s", tax_names)

    # extend tax_ids to include IDs of children
    tax_ids_to_include = set()
    for tax_id in tax_ids:
        tax_ids_to_include.add(tax_id)
        if not omit_children:
            child_ids = collect_children(db.children, set([tax_id]))
            tax_ids_to_include |= set(child_ids)

    tax_ids_to_include = frozenset(tax_ids_to_include) # frozenset membership check slightly faster
    log.info("matching against {} taxa".format(len(tax_ids_to_include)))

    # perform the actual filtering to return a list of read IDs, writeen to a temp file
    with util.file.tempfname(".txt.gz") as temp_read_list:
        with util.file.open_or_gzopen(temp_read_list, "wt") as read_IDs_file:
            read_ids_written = 0
            for row in util.file.read_tabfile(read_IDs_to_tax_IDs):
                assert tax_id_col<len(row), "tax_id_col does not appear to be in range for number of columns present in mapping file"
                assert read_id_col<len(row), "read_id_col does not appear to be in range for number of columns present in mapping file"
                read_id = row[read_id_col]
                read_tax_id = int(row[tax_id_col])

                # transform read ID to take read pairs into account
                read_id_match = re.match(paired_read_base_pattern,read_id)
                if (read_id_match and
                    read_tax_id in tax_ids_to_include):
                    read_IDs_file.write(read_id_match.group(1)+"\n")
                    read_ids_written+=1
        log.info("matched {} reads".format(read_ids_written))

        # report count if desired
        if out_count:
            with open(out_count, 'wt') as outf:
                outf.write("{}\n".format(read_ids_written))

        # if we found reads matching the taxNames requested,
        if (read_ids_written > 0) or exclude:
            # filter the input bam to include only these
            tools.picard.FilterSamReadsTool().execute(in_bam,
                                                        exclude,
                                                        temp_read_list,
                                                        out_bam,
                                                        JVMmemory=JVMmemory)
        else:
            # otherwise, "touch" the output bam to contain the
            # header of the input bam (no matching reads)
            tools.samtools.SamtoolsTool().dumpHeader(in_bam,out_bam)
__commands__.append(('filter_bam_to_taxa', parser_filter_bam_to_taxa))



def parser_kraken_taxlevel_summary(parser=argparse.ArgumentParser()):
    parser.add_argument('summary_files_in', nargs="+", help='Kraken-format summary text file with tab-delimited taxonomic levels.')
    parser.add_argument('--jsonOut', dest="json_out", type=argparse.FileType('w'), help='The path to a json file containing the relevant parsed summary data in json format.')
    parser.add_argument('--csvOut', dest="csv_out", type=argparse.FileType('w'), help='The path to a csv file containing sample-specific counts.')
    parser.add_argument('--taxHeading', nargs="+", dest="tax_headings", help='The taxonomic heading to analyze (default: %(default)s). More than one can be specified.', default="Viruses")
    parser.add_argument('--taxlevelFocus', dest="taxlevel_focus", help='The taxonomic heading to summarize (totals by Genus, etc.) (default: %(default)s).', default="species")#,
                        #choices=["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"])
    parser.add_argument('--topN', type=int, dest="top_n_entries", help='Only include the top N most abundant taxa by read count (default: %(default)s)', default=100)
    parser.add_argument('--countThreshold', type=int, dest="count_threshold", help='Minimum number of reads to be included (default: %(default)s)', default=1)
    parser.add_argument('--zeroFill', action='store_true', dest="zero_fill", help='When absent from a sample, write zeroes (rather than leaving blank).')
    parser.add_argument('--noHist', action='store_true', dest="no_hist", help='Write out a report by-sample rather than a histogram.')
    parser.add_argument('--includeRoot', action='store_true', dest="include_root", help='Include the count of reads at the root level and the unclassified bin.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, taxlevel_summary, split_args=True)
    return parser

def taxlevel_summary(summary_files_in, json_out, csv_out, tax_headings, taxlevel_focus, top_n_entries, count_threshold, no_hist, zero_fill, include_root):
    """
        Aggregates taxonomic abundance data from multiple Kraken-format summary files.
        It is intended to report information on a particular taxonomic level (--taxlevelFocus; ex. 'species'),
        within a higher-level grouping (--taxHeading; ex. 'Viruses'). By default, when --taxHeading
        is at the same level as --taxlevelFocus a summary with lines for each sample is emitted.
        Otherwise, a histogram is returned. If per-sample information is desired, --noHist can be specified.
        In per-sample data, the suffix "-pt" indicates percentage, so a value of 0.02 is 0.0002 of the total number of reads for the sample.
        If --topN is specified, only the top N most abundant taxa are included in the histogram count or per-sample output.
        If a number is specified for --countThreshold, only taxa with that number of reads (or greater) are included.
        Full data returned via --jsonOut (filtered by --topN and --countThreshold), whereas -csvOut returns a summary.
    """

    samples = {}
    same_level = False

    Abundance = collections.namedtuple("Abundance", "percent,count,kmers,dup,cov")

    def indent_len(in_string):
        return len(in_string)-len(in_string.lstrip())

    for f in list(summary_files_in):
        sample_name, extension = os.path.splitext(f)
        sample_summary = {}
        sample_root_summary = {}
        tax_headings_copy = [s.lower() for s in tax_headings]

        # -----------------------------------------------------------------
        # KrakenUniq has two lines prefixed by '#', a blank line,
        # and then a TSV header beginning with "%". The column fields are:
        # (NB:field names accurate, but space-separated in this comment 
        # for readability here)
        #   %        reads  taxReads  kmers  dup   cov  taxID  rank          taxName
        #   0.05591  2      0         13     1.85  NA   10239  superkingdom  Viruses
        #
        # Where the fields are:
        #   %:
        #   reads:
        #   taxReads:
        #   kmers: number of unique k-mers
        #   dup: average number of times each unique k-mer has been seen
        #   cov: coverage of the k-mers of the clade in the database
        #   taxID: 
        #   rank: row["rank"]; A rank code (see list below)
        #   taxName: row["sci_name"]; indented scientific name
        #
        # Taxonomic ranks used by KrakenUniq include:
        #   unknown, no rank, sequence, assembly, subspecies, 
        #   species, species subgroup, species group, subgenus, 
        #   genus, tribe, subfamily, family, superfamily, parvorder, 
        #   infraorder, suborder, order, superorder, parvclass, 
        #   infraclass, subclass, class, superclass, subphylum, 
        #   phylum, kingdom, superkingdom, root
        #
        #   via: https://github.com/fbreitwieser/krakenuniq/blob/a8b4a2dbf50553e02d3cab3c32f93f91958aa575/src/taxdb.hpp#L96-L131
        # -----------------------------------------------------------------
        # Kraken (standard) reports lack header lines. 
        # (NB:field names below are only for reference. Space-separated for 
        # readability here)
        #   %     reads  taxReads  rank  taxID      taxName
        #   0.00  16     0         D     10239      Viruses
        #
        # Where the fields are:
        #   %:        row["pct_of_reads"]; Percentage of reads covered by the clade rooted at this taxon
        #   reads:    row["num_reads"]; Number of reads covered by the clade rooted at this taxon
        #   taxReads: row["reads_exc_children"]; Number of reads assigned directly to this taxon
        #   rank:     row["rank"]; A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
        #   taxID:    row["NCBI_tax_ID"]; NCBI taxonomy ID
        #   taxName:  row["sci_name"]; indented scientific name
        # -----------------------------------------------------------------


        with util.file.open_or_gzopen(f, 'rt') as inf:
            report_type=None
            should_process = False
            indent_of_selection = -1
            currently_being_processed = ""
            for lineno, line in enumerate(inf):
                if len(line.rstrip('\r\n').strip()) == 0 or ( report_type != None and line.startswith("#") or line.startswith("%")):
                    continue

                # KrakenUniq is mentioned on the first line of
                # summary reports created by KrakenUniq
                if not report_type and "KrakenUniq" in line:
                    report_type="krakenuniq"
                    continue
                elif not report_type:
                    report_type="kraken"                

                csv.register_dialect('kraken_report', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
                if report_type == "kraken":
                    fieldnames = [ "pct_of_reads",
                                    "num_reads",
                                    "reads_exc_children",
                                    "rank",
                                    "NCBI_tax_ID",
                                    "sci_name"
                                ]
                elif report_type == "krakenuniq":
                    fieldnames = [ 
                                    "pct_of_reads",
                                    "num_reads",
                                    "reads_exc_children",
                                    "uniq_kmers",
                                    "kmer_dups",
                                    "cov_of_clade_kmers",
                                    "NCBI_tax_ID",
                                    "rank",
                                    "sci_name"
                                ]
                else:
                    continue #never reached since we fall back to kraken above

                row = next(csv.DictReader([line.strip().rstrip('\n')], fieldnames=fieldnames, dialect="kraken_report"))

                try:
                    indent_of_line = indent_len(row["sci_name"])
                except AttributeError as e:
                    log.warning("Report type: '{}'".format(report_type))
                    log.warning("Issue with line {}: '{}'".format(lineno,line.strip().rstrip('\n')))
                    log.warning("From file: {}".format(f))
                # remove leading/trailing whitespace from each item
                row = { k:v.strip() for k, v in row.items()}

                # rows are formatted as described above. 
                # Kraken:
                #   0.00  16  0   D   10239     Viruses
                # KrakenUniq:
                #   0.05591  2      0         13     1.85  NA   10239  superkingdom  Viruses

                # if the root-level bins (root, unclassified) should be included, do so, but bypass normal
                # stateful parsing logic since root does not have a distinct rank level
                if row["sci_name"].lower() in ["root","unclassified"] and include_root:
                    sample_root_summary[row["sci_name"]] = collections.OrderedDict()
                    sample_root_summary[row["sci_name"]][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]),row.get("kmers",None),row.get("dup",None),row.get("cov",None))
                    continue

                if indent_of_line <= indent_of_selection:
                    should_process = False
                    indent_of_selection=-1

                if indent_of_selection == -1:
                    if row["sci_name"].lower() in tax_headings_copy:
                        tax_headings_copy.remove(row["sci_name"].lower())

                        should_process = True
                        indent_of_selection = indent_of_line
                        currently_being_processed = row["sci_name"]
                        sample_summary[currently_being_processed] = collections.OrderedDict()
                        if row["rank"] == rank_code(taxlevel_focus) or row["rank"].lower().replace(" ","") == taxlevel_focus.lower().replace(" ",""):
                            same_level = True
                        if row["rank"] in ("-","no rank"):
                            log.warning("Non-taxonomic parent level selected")

                if should_process:
                    # skip "-" rank levels since they do not occur at the sample level
                    # otherwise include the taxon row if the rank matches the desired level of focus
                    if (row["rank"] not in ("-","no rank") and (rank_code(taxlevel_focus) == row["rank"] or row["rank"].lower().replace(" ","") == taxlevel_focus.lower().replace(" ","")) ):
                        if int(row["num_reads"])>=count_threshold:
                            sample_summary[currently_being_processed][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]),row.get("kmers",None),row.get("dup",None),row.get("cov",None))


        for k,taxa in sample_summary.items():
            sample_summary[k] = collections.OrderedDict(sorted(taxa.items(), key=lambda item: (item[1][1]) , reverse=True)[:top_n_entries])

            if len(list(sample_summary[k].items()))>0:
                log.info("{f}: most abundant among {heading} at the {level} level: "
                            "\"{name}\" with {reads} reads ({percent:.2%} of total); "
                            "included since >{threshold} read{plural}".format(
                                                                          f=f,
                                                                          heading=k,
                                                                          level=taxlevel_focus,
                                                                          name=list(sample_summary[k].items())[0][0],
                                                                          reads=list(sample_summary[k].items())[0][1].count,
                                                                          percent=list(sample_summary[k].items())[0][1].percent/100.0,
                                                                          threshold=count_threshold,
                                                                          plural="s" if count_threshold>1 else "" )
                )

        if include_root:
            # include root-level bins (root, unclassified) in the returned data
            for k,taxa in sample_root_summary.items():
                assert (k not in sample_summary), "{k} already in sample summary".format(k=k)
                sample_summary[k] = taxa
        samples[sample_name] = sample_summary

    if json_out != None:
        json_summary = json.dumps(samples, sort_keys=True, indent=4, separators=(',', ': '))
        json_out.write(json_summary)
        json_out.close()


    if csv_out != None:

        # if we're writing out at the same level as the query header
        # write out the fractions and counts
        if same_level or no_hist:

            fieldnames = set()
            for sample, taxa in samples.items():
                for heading,taxon in taxa.items():
                    if len(taxon):
                        for k in taxon.keys():
                            fieldnames |= set([k+"-pt",k+"-ct"])

            heading_columns = ["sample"]
            if include_root:
                root_fields = ["root-pt","root-ct","unclassified-pt","unclassified-ct"]
                fieldnames -= set(root_fields)
                heading_columns += root_fields

            writer = csv.DictWriter(csv_out, restval=0 if zero_fill else '', fieldnames=heading_columns+sorted(list(fieldnames)))
            writer.writeheader()

            for sample, taxa in samples.items():
                sample_dict = {}
                sample_dict["sample"] = sample
                for heading,taxon in taxa.items():
                    for entry in taxon.keys():
                        sample_dict[entry+"-pt"] = taxon[entry].percent
                        sample_dict[entry+"-ct"] = taxon[entry].count
                writer.writerow(sample_dict)


            csv_out.close()

        # otherwise write out a histogram
        else:
            count = 0
            summary_counts = collections.defaultdict(dict)
            for sample, totals in samples.items():
                for heading,taxa in totals.items():
                    for taxon in taxa.keys():
                        if taxon not in summary_counts[heading].keys():
                            summary_counts[heading][taxon] = 1
                        else:
                            summary_counts[heading][taxon] += 1

            for k,taxa in summary_counts.items():
                summary_counts[k] = collections.OrderedDict(sorted(taxa.items(), key=lambda item: (item[1]) , reverse=True))


            fieldnames = ["heading","taxon","num_samples"]
            writer = csv.DictWriter(csv_out, restval=0 if zero_fill else '', fieldnames=fieldnames)
            writer.writeheader()

            for heading,taxa_counts in summary_counts.items():
                writer.writerows([{"heading":heading,"taxon":taxon,"num_samples":count} for taxon,count in taxa_counts.items()])

            csv_out.close()


__commands__.append(('taxlevel_summary', parser_kraken_taxlevel_summary))


def parser_taxlevel_plurality(parser=argparse.ArgumentParser()):
    parser.add_argument('summary_file', help='input Kraken-format summary text file with tab-delimited taxonomic levels.')
    parser.add_argument('tax_heading', help='The taxonomic heading to analyze.')
    parser.add_argument('out_report', help='tab-delimited output file.')
    parser.add_argument('--min_reads', type=int, dest="min_reads", help='Only include hits with more than min_reads (default: %(default)s)', default=1)
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, taxlevel_plurality, split_args=True)
    return parser

def taxlevel_plurality(summary_file, tax_heading, out_report, min_reads):
    """
        Identifies the most abundant taxon (of any rank) contributing to a node of interest in the taxonomic tree.
        It is intended to highlight the primary contributor of taxonomic signal within a taxonomic category of interest,
        for example, the most abundant virus among all viruses.
    """

    keeper_rows = []
    with util.file.open_or_gzopen(summary_file, 'rt') as inf:
        report_type=None
        should_process = False
        indent_of_selection = -1
        def indent_len(in_string):
            return len(in_string)-len(in_string.lstrip())
        for lineno, line in enumerate(inf):
            if len(line.rstrip('\r\n').strip()) == 0 or ( report_type != None and line.startswith("#") or line.startswith("%")):
                continue

            if not report_type and "KrakenUniq" in line:
                report_type="krakenuniq"
                continue
            elif not report_type:
                report_type="kraken"

            csv.register_dialect('kraken_report', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
            if report_type == "kraken":
                fieldnames = [  "pct_of_total",
                                "reads_cumulative",
                                "reads_excl_children",
                                "taxon_rank",
                                "taxon_id",
                                "taxon_sci_name"
                            ]
            elif report_type == "krakenuniq":
                fieldnames = [  "pct_of_total",
                                "reads_cumulative",
                                "reads_excl_children",
                                "uniq_kmers",
                                "kmer_dups",
                                "cov_of_clade_kmers",
                                "taxon_id",
                                "taxon_rank",
                                "taxon_sci_name"
                            ]
            else:
                continue #never reached since we fall back to kraken above

            row = next(csv.DictReader([line.strip().rstrip('\n')], fieldnames=fieldnames, dialect="kraken_report"))

            try:
                indent_of_line = indent_len(row["taxon_sci_name"])
            except AttributeError as e:
                log.warning("Report type: '{}'".format(report_type))
                log.warning("Issue with line {}: '{}'".format(lineno,line.strip().rstrip('\n')))
                log.warning("From file: {}".format(f))
            # remove leading/trailing whitespace from each item
            row = { k:v.strip() for k, v in row.items()}

            if indent_of_line <= indent_of_selection:
                should_process = False
                indent_of_selection=-1

            if indent_of_selection == -1:
                if row["taxon_sci_name"].lower() == tax_heading.lower():
                    should_process = True
                    indent_of_selection = indent_of_line

            if should_process:
                keeper_rows.append(row)

    # find the top hits within the focal taxon (if any)
    out = []
    if not keeper_rows:
        # we didn't find anything under the tax_heading of interest
        out.append({'order_within_focal':0, 'focal_taxon_name':tax_heading, 'focal_taxon_count':0, 'pct_of_focal':0.0, 'pct_of_total':0.0, 'reads_cumulative':0, 'reads_excl_children':0, 'taxon_rank':'', 'taxon_id':'', 'taxon_sci_name':''})
    else:
        # find the most abundant taxon classified underneath (but not including) the tax_heading 
        assert keeper_rows[0]["taxon_sci_name"].lower() == tax_heading.lower()
        keepers_sorted = enumerate(sorted(keeper_rows[1:], key=lambda row:int(row["reads_excl_children"]), reverse=True))
        for i, hit in keepers_sorted:
            outrow = {  'order_within_focal': i+1,
                        'focal_taxon_name':keeper_rows[0]["taxon_sci_name"],
                        'focal_taxon_count':keeper_rows[0]["reads_cumulative"],
                        'pct_of_focal': 100.0 * float(hit["reads_excl_children"]) / float(keeper_rows[0]["reads_cumulative"])}
            for k in ("pct_of_total", "reads_cumulative", "reads_excl_children", "taxon_rank", "taxon_id", "taxon_sci_name"):
                outrow[k] = hit[k]
            if int(outrow['reads_excl_children']) >= min_reads:
                out.append(outrow)

    # write outputs
    with util.file.open_or_gzopen(out_report, 'wt') as outf:
        header = ('focal_taxon_name', 'focal_taxon_count', 'order_within_focal', 'pct_of_focal', 'pct_of_total', 'reads_cumulative', 'reads_excl_children', 'taxon_rank', 'taxon_id', 'taxon_sci_name')
        writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(out)

__commands__.append(('taxlevel_plurality', parser_kraken_taxlevel_plurality))

def parser_kb_extract(parser=argparse.ArgumentParser()):
    """Argument parser for the kb_python extract command.

    Args:
        parser (argparse.ArgumentParser): Argument parser instance. Defaults to argparse.ArgumentParser().

    Returns:
        argparse.ArgumentParser: The parser with arguments added.
    """
    parser.add_argument('inBam', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--index', help='kb index file.')
    parser.add_argument('--t2g', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outDir', help='Output directory (default: kb_out)', default='kb_out')
    parser.add_argument('--aa', action='store_true', help='True if sequence contains amino acids(default: False).')
    parser.add_argument('--targets', nargs='+', help='List of target sequences to extract from input sequences.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken2, split_args=True)
    return parser
def kb_extract(inBams, index_file, t2g_file, target_ids, aa=False, outDir=None, threads=None):
    """Runs kb count on the input BAM files.

    Args:
        inBams (list): List of input BAM files.
        outDir (str): Output directory. Defaults to None.
        index_file (str): Path to the kb index file.
        t2g_file (str): Path to the transcript-to-gene mapping file.
        kmer_size (int, optional): K-mer size for the alignment. Defaults to 31.
        technology (str, optional): Sequencing technology used. Defaults to 'bulk'.
        h5ad (bool, optional): Whether to output HDF5 file. Defaults to False.
        loom (bool, optional): Whether to output Loom file. Defaults to False.
        threads (int, optional): Number of threads to use. Defaults to None.
    """

    assert outDir, ('Output directory must be specified.')
    kb_tool = classify.kb.kb()
    kb_tool.extract(
        inBams=inBams,
        outDir=outDir,
        index_file=index_file,
        t2g_file=t2g_file,
        target_ids=target_ids,
        aa=aa,
        threads=threads
    )
__commands__.append(('kb_extract', parser_kb_extract))

def parser_kb_top_taxa(parser=argparse.ArgumentParser()):
    parser.add_argument('counts_tar', help='input Kraken-format summary text file with tab-delimited taxonomic levels.')
    parser.add_argument('out_report', help='tab-delimited output file.')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kb_top_taxa, split_args=True)
    return parser

def kb_top_taxa(counts_tar, id_to_tax_map, target_taxon, out_report):
    """Identifies the most abundant taxon (of any rank) contributing to a taxa node of interest in kb count output. This function
    is currently only set to parse h5ad files and further functionality will be added in the future.
        
    It is intended to highlight the primary contributor of taxonomic signal within a taxonomic category of interest,
    for example, the most abundant virus among all viruses.
        
    Args:
        counts_tar (str): Path to the input counts file (HDF5 format).
        id_to_tax_map (str): Path to the ID to taxonomy mapping file.
        target_taxon (str): The taxonomic heading to analyze.
        out_report (str): Path to the output report file.
    """
    sotu_counts = pd.DataFrame(columns=['sotu_id', 'count', 'tax_name', 'tax_id', 'tax_category'])
    
    tax_map_df = pd.read_csv(id_to_tax_map)
    tax_map_df = tax_map_df.replace('.', np.nan)
    
    tax_map_df = tax_map_df.apply(lambda row: row[row.first_valid_index():].dropna().iloc[-2], axis=1)
    tax_map_df['farthest_taxa'] = tax_map_df.apply(lambda r: r.dropna().iloc[-2] if r.notna().any() else pd.NA, axis=1)
    tax_map_df['farthest_taxa'] = tax_map_df['farthest_taxa'].fillna('root')
    
    pattern = r'^(?:\d+|[A-Za-z]{1,3}\d{2,}|[A-Za-z]_?\d{3,})$'
    mask = (
        tax_map_df['farthest_taxa']
        .astype(str)
        .str.fullmatch(pattern, na=False)
    )
    tax_map_df.loc[mask, 'farthest_taxa'] = "Unknown"
    
    ## Start by extracting our the contents of our tarball
    with util.file.tmp_dir() as tmp_dir:
        util.file.untar(counts_tar, tmp_dir)
        h5ad_files = glob.glob(os.path.join(tmp_dir, "counts_unfiltered", '*.h5ad'))
        
        assert len(h5ad_files) == 1, "Expected exactly one .h5ad file in the counts tarball, found {}".format(len(h5ad_files))
        h5ad_file = h5ad_files[0]
        
        adata = anndata.read_h5ad(h5ad_file)
        
        counts_mtx = adata.X.to_array()
        sotu_labels = adata.var.index
        sotu_totals = list(zip(sotu_labels, counts_mtx.sum(axis=0)))
        
        tax_mapping = [(tax_map_df.loc[tax_map_df['id'] == sotu_id, 'root'].values[0], 
                        tax_map_df.loc[tax_map_df['id'] == sotu_id, 'furthest_taxa'].values[0])
                       for (sotu_id, _count) in sotu_totals]
        
        sotu_counts = pd.concat([sotu_counts, list(zip(*sotu_totals, *zip(*tax_mapping)))])
        sotu_counts.query('tax_category == @target_taxon', inplace=True)
        
        if sotu_counts.empty:
            out = [{'focal_taxon_name':target_taxon, 'focal_taxon_count':0, 'pct_of_focal':0.0, 'hit_id':'', 'hit_lowest_taxa_name':'', 'hit_reads': 0}]
        else:
            target_taxon_count = sotu_counts.groupby('tax_category')['count'].sum().values[0]
            
            sotu_counts = sotu_counts.sort_values(by='count', ascending=False).reset_index(drop=True)
            out = []
            for i, row in sotu_counts.iterrows():
                outrow = {  'focal_taxon_name': target_taxon,
                            'focal_taxon_count': target_taxon_count,
                            'hit_id': row['sotu_id'],
                            'hit_lowest_taxa_name': row['furthest_taxa'],
                            'hit_reads': row['count'],
                            'pct_of_focal': 100.0 * float(row['count']) / float(target_taxon_count),
                }
                out.append(outrow)
    
    # write outputs
    with util.file.open_or_gzopen(out_report, 'wt') as outf:
        header = ('focal_taxon_name', 'focal_taxon_count', 'hit_id', 'hit_lowest_taxa_name', 'hit_reads', 'pct_of_focal')
        writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(out)

__commands__.append(('kb_top_taxa', parser_kb_top_taxa))


def parser_krona_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Krona taxonomy database output directory.')
    parser.add_argument('--taxdump_tar_gz', help='NCBI taxdump.tar.gz file', default=None)
    parser.add_argument('--get_accessions', action='store_true', help='Fetch NCBI accession to taxid mappings. This is not required for processing kraken1/2/uniq hits, only for BLAST hits, and adds a significant amount of time and database space (default false).')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, krona_build, split_args=True)
    return parser
def krona_build(db, taxdump_tar_gz=None, get_accessions=False):
    '''
    Builds a Krona taxonomy database
    '''
    classify.krona.Krona().build_db(
        db, taxdump_tar_gz=taxdump_tar_gz, get_accessions=get_accessions)
__commands__.append(('krona_build', parser_krona_build))


def parser_kraken2_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database output directory.')
    parser.add_argument('--tax_db', help='Use pre-existing kraken2 taxonomy db structure', default=None)
    parser.add_argument('--taxdump_out', help='Save ncbi taxdump.tar.gz file', default=None)
    parser.add_argument('--standard_libraries',
        nargs='+',
        choices=["archaea", "bacteria", "plasmid",
             "viral", "human", "fungi", "plant", "protozoa",
             "nr", "nt", "env_nr", "env_nt", "UniVec",
             "UniVec_Core"],
        help='A list of "standard" kraken libraries to download on the fly and add.')
    parser.add_argument('--custom_libraries', nargs='+', help='Custom fasta files with properly formatted headers.')
    parser.add_argument('--kmerLen', type=int, help='(k) k-mer length (kraken2 default: 35nt/15aa)')
    parser.add_argument('--minimizerLen', type=int, help='(l) Minimizer length (kraken2 default: 31nt/12aa)')
    parser.add_argument('--minimizerSpaces', type=int, help='(s) Number of characters in minimizer that are ignored in comparisons (kraken2 default: 7nt/0aa)')
    parser.add_argument('--protein', action='store_true', help='Build protein database (default false/nucleotide).')
    parser.add_argument('--maxDbSize', type=int, help='Maximum db size in GB (default: none)')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken2_build, split_args=True)
    return parser
def kraken2_build(db,
                tax_db=None, taxdump_out=None,
                standard_libraries=(), custom_libraries=(),
                protein=False,
                minimizerLen=None, kmerLen=None, minimizerSpaces=None,
                maxDbSize=None, threads=None):
    '''
    Builds a kraken2 database from library directory of fastas and taxonomy db
    directory. The --subsetTaxonomy option allows shrinking the taxonomy to only
    include taxids associated with the library folders. For this to work, the
    library fastas must have the standard id names such as `>NC1234.1`
    accessions, `>gi|123456789|ref|XXXX||`, or custom kraken name
    `>kraken:taxid|1234|`.
    '''

    kraken_tool = classify.kraken2.Kraken2()
    kraken_tool.build(db,
        tax_db=tax_db,
        standard_libraries=standard_libraries,
        custom_libraries=custom_libraries,
        taxdump_out=taxdump_out,
        k=kmerLen, l=minimizerLen, s=minimizerSpaces,
        protein=protein,
        max_db_size=maxDbSize,
        num_threads=threads)

__commands__.append(('kraken2_build', parser_kraken2_build))


def parser_kb_build(parser=argparse.ArgumentParser()):
    parser.add_argument('ref_fasta', help='Reference sequence fasta file.')
    parser.add_argument('--index', help='kb index output file.')
    parser.add_argument('--workflow', choices=['standard', 'nac', 'kite', 'custom'],
                        default='standard', help='Type of index to create (default: %(default)s).')
    parser.add_argument('--kmer_len', type=int, help='k-mer length (default: 31).')
    parser.add_argument('--aa', action='store_true', help='True if sequence contains amino acids(default: False).')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kb_build, split_args=True)
    return parser
def kb_build(ref_fasta, index, workflow='standard', kmer_len=31, aa_seq=False, threads=None):
    '''
    Builds a kb index from a reference fasta file.
    
    Args:
        ref_fasta (str): Path to the reference sequence fasta file.
        index (str): Path to the output kb index file.
        workflow (str): Type of index to create. Options are 'standard', 'nac', 'kite', 'custom'.
        kmer_len (int): k-mer length (default: 31).
        aa_seq (bool): True if sequence contains amino acids (default: False).
        threads (int): Number of threads to use (default: None).
    '''

    kb_tool = classify.kb.kb()
    kb_tool.build(ref_fasta,
                        index=index,
                        workflow=workflow,
                        kmer_len=kmer_len,
                        aa_seq=aa_seq,
                        threads=threads)
    
__commands__.append(('kb_build', parser_kb_build))


def parser_krakenuniq_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Krakenuniq database output directory.')
    parser.add_argument('--library', help='Input library directory of fasta files. If not specified, it will be read from the "library" subdirectory of "db".')
    parser.add_argument('--taxonomy', help='Taxonomy db directory. If not specified, it will be read from the "taxonomy" subdirectory of "db".')
    parser.add_argument('--subsetTaxonomy', action='store_true', help='Subset taxonomy based on library fastas.')
    parser.add_argument('--minimizerLen', type=int, help='Minimizer length (krakenuniq default: 15)')
    parser.add_argument('--kmerLen', type=int, help='k-mer length (krakenuniq default: 31)')
    parser.add_argument('--maxDbSize', type=int, help='Maximum db size in GB (will shrink if too big)')
    parser.add_argument('--clean', action='store_true', help='Clean by deleting other database files after build')
    parser.add_argument('--workOnDisk', action='store_true', help='Work on disk instead of RAM. This is generally much slower unless the "db" directory lives on a RAM disk.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, krakenuniq_build, split_args=True)
    return parser
def krakenuniq_build(db, library, taxonomy=None, subsetTaxonomy=None,
                     threads=None, workOnDisk=False,
                     minimizerLen=None, kmerLen=None, maxDbSize=None, clean=False):
    '''
    Builds a krakenuniq database from library directory of fastas and taxonomy
    db directory. The --subsetTaxonomy option allows shrinking the taxonomy to
    only include taxids associated with the library folders. For this to work,
    the library fastas must have the standard accession id names such as
    `>NC1234.1` or `>NC_01234.1`.

    Setting the --minimizerLen (default: 16) small, such as 10, will drastically
    shrink the db size for small inputs, which is useful for testing.

    The built db may include symlinks to the original --library / --taxonomy
    directories. If you want to build a static archiveable version of the
    library, simply use the --clean option, which will also remove any
    unnecessary files.
    '''
    util.file.mkdir_p(db)
    library_dir = os.path.join(db, 'library')
    library_exists = os.path.exists(library_dir)
    if library:
        try:
            os.symlink(os.path.abspath(library), os.path.join(db, 'library'))
        except FileExistsError:
            pass
    else:
        if not library_exists:
            raise FileNotFoundError('Library directory {} not found'.format(library_dir))

    taxonomy_dir = os.path.join(db, 'taxonomy')
    taxonomy_exists = os.path.exists(taxonomy_dir)
    if taxonomy:
        if taxonomy_exists:
            raise KrakenUniqBuildError('Output db directory already contains taxonomy directory {}'.format(taxonomy_dir))
        if subsetTaxonomy:
            accessions = fasta_library_accessions(library)

            whitelist_accession_f = util.file.mkstempfname()
            with open(whitelist_accession_f, 'wt') as f:
                for accession in accessions:
                    print(accession, file=f)

            # Context-managerize eventually
            taxonomy_tmp = tempfile.mkdtemp()
            subset_taxonomy(taxonomy, taxonomy_tmp, whitelistAccessionFile=whitelist_accession_f)
            shutil.move(taxonomy_tmp, taxonomy_dir)
        else:
            os.symlink(os.path.abspath(taxonomy), taxonomy_dir)
    else:
        if not taxonomy_exists:
            raise FileNotFoundError('Taxonomy directory {} not found'.format(taxonomy_dir))
        if subsetTaxonomy:
            raise KrakenUniqBuildError('Cannot subset taxonomy if already in db folder')

    krakenuniq_tool = classify.kraken.KrakenUniq()
    options = {'--build': None}
    if threads:
        options['--threads'] = threads
    if minimizerLen:
        options['--minimizer-len'] = minimizerLen
    if kmerLen:
        options['--kmer-len'] = kmerLen
    if maxDbSize:
        options['--max-db-size'] = maxDbSize
    if workOnDisk:
        options['--work-on-disk'] = None
    krakenuniq_tool.build(db, options=options)

    if clean:
        krakenuniq_tool.execute('krakenuniq-build', db, '', options={'--clean': None})
# KrakenUniq disabled from future versions for now, pending conda rebuild of @yesimon's custom fork & dependencies
#__commands__.append(('krakenuniq_build', parser_krakenuniq_build))

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

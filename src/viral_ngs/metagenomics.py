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

from .core import cmd
from .core import file
from .core import misc
from . import read_utils

from .classify import kma
from .classify import kraken2
from .classify import krona
from .classify import kb
from .classify.taxonomy import (
    TaxIdError, TaxonomyDb, BlastRecord, blast_records, paired_query_id,
    blast_m8_taxids, extract_tax_id, coverage_lca, parents_to_children,
    maybe_compressed, rank_code
)

__commands__ = []

log = logging.getLogger(__name__)


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
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, subset_taxonomy, split_args=True)
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
    file.mkdir_p(os.path.join(outputDb, 'accession2taxid'))
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
        with file.open_or_gzopen(input_path, 'rt') as f, \
             file.open_or_gzopen(output_path, 'wt') as out_f:
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

    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, filter_taxids_to_focal_hits, split_args=True)
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
    with file.open_or_gzopen(focal_report_tsv, "rt") as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
            if int(row['reads_excl_children']) >= min_read_count:
                hits.add(int(row['taxon_id']))

    # filter taxids_tsv -> output_tsv
    with file.open_or_gzopen(taxids_tsv, "rt") as inf:
        with file.open_or_gzopen(output_tsv, "wt") as outf:
            for line in inf:
                taxid = int(line.rstrip('\r\n').split('\t')[0])
                ancestors = taxdb.get_ordered_ancestors(taxid)
                if any(node in hits for node in [taxid] + ancestors):
                    outf.write(line)

__commands__.append(('filter_taxids_to_focal_hits', parser_filter_taxids_to_focal_hits))


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
    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kraken2, split_args=True)
    return parser
def kraken2(db, inBams, outReports=None, outReads=None, min_base_qual=None, confidence=None, minimum_hit_groups=None, threads=None):
    '''
        Classify reads by taxon using Kraken2
    '''

    assert outReads or outReports, ('Either --outReads or --outReport must be specified.')
    kraken_tool = kraken2.Kraken2()
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
    parser.add_argument('in_bam', help='Input unaligned reads, BAM format.')
    parser.add_argument('--index', help='kb index file.')
    parser.add_argument('--t2g', help='Input unaligned reads, BAM format.')
    parser.add_argument('--kmer_len', type=int, help='k-mer size (default: 31bp)', default=31)
    parser.add_argument('--parity', choices=['single', 'paired'], help='Library parity (default: single)', default='single')
    parser.add_argument('--technology', choices=['10xv2', '10xv3', '10xv3-3prime', '10xv3-5prime', 'dropseq', 
                                                 'indrop', 'celseq', 'celseq2', 'smartseq2', 'bulk'], 
                        help='Technology used to generate the data (default: bulk)', default='bulk')
    parser.add_argument('--h5ad', action='store_true', help='Output HDF5 file (default: False)', default=False)
    parser.add_argument('--loom', action='store_true', help='Output Loom file (default: False)', default=False)
    parser.add_argument('--protein', action='store_true', help='True if sequence contains amino acids (default: False).')
    parser.add_argument('--out_dir', help='Output directory (default: kb_out)', default='kb_out')
    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kb_python, split_args=True)
    return parser
def kb_python(in_bam, index=None, t2g=None, kmer_len=31, parity='single', technology='bulk', h5ad=False, loom=False, protein=False, out_dir=None, threads=None):
    """Runs kb count on the input BAM files.

    Args:
        in_bam (list): List of input BAM files.
        out_dir (str): Output directory. Defaults to None.
        index (str): Path to the kb index file.
        t2g (list|str): Transcript-to-gene mapping file(s).
        kmer_len (int, optional): K-mer size for the alignment. Defaults to 31.
        parity (str, optional): Library parity (default: single). Defaults to 'single'.
        technology (str, optional): Sequencing technology used. Defaults to 'bulk'.
        h5ad (bool, optional): Whether to output HDF5 file. Defaults to False.
        loom (bool, optional): Whether to output Loom file. Defaults to False.
        protein (bool, optional): Whether the sequence contains amino acids. Defaults to False.
        threads (int, optional): Number of threads to use. Defaults to None.
    """

    assert out_dir, ('Output directory must be specified.')
    kb_tool = kb.kb()
    kb_tool.classify(
        in_bam=in_bam,
        out_dir=out_dir,
        index_file=index,
        t2g_file=t2g,
        k=kmer_len,
        parity=parity,
        technology=technology,
        h5ad=h5ad,
        loom=loom,
        protein=protein,
        num_threads=threads
    )
__commands__.append(('kb', parser_kb))

def parser_kma(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='KMA database prefix.')
    parser.add_argument('inBams', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outPrefixes', nargs='+', help='KMA output prefixes.')
    parser.add_argument('--threads', type=int, help='Number of threads.')
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, main_kma, split_args=True)
    return parser

def main_kma(db, inBams, outPrefixes=None, threads=None):
    if outPrefixes and len(inBams) != len(outPrefixes):
        raise ValueError(f"Number of input BAMs ({len(inBams)}) must match number of output prefixes ({len(outPrefixes)})")
    kma_tool = kma.KMA()
    for in_bam, out_prefix in itertools.zip_longest(inBams, outPrefixes):
        kma_tool.classify(in_bam, db, out_prefix, num_threads=threads)

__commands__.append(('kma', parser_kma))

def parser_kma_build(parser=argparse.ArgumentParser()):
    parser.add_argument('ref_fasta', help='Reference FASTA file.')
    parser.add_argument('db_prefix', help='Output database prefix.')
    parser.add_argument('--threads', type=int, help='Number of threads.')
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, main_kma_build, split_args=True)
    return parser

def main_kma_build(ref_fasta, db_prefix, threads=None):
    kma_tool = kma.KMA()
    kma_tool.build(ref_fasta, db_prefix, num_threads=threads)

__commands__.append(('kma_build', parser_kma_build))


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
    parser.add_argument('--inputType', help='Handling for specialized report types.', default='tsv', choices=['tsv', 'kraken2'])
    cmd.common_args(parser, (('loglevel', None), ('version', None)))
    cmd.attach_main(parser, krona, split_args=True)
    return parser
def krona(inReports, db, outHtml, queryColumn=None, taxidColumn=None, scoreColumn=None, magnitudeColumn=None, noHits=None, noRank=None,
          inputType=None, sample_name=None):
    '''
        Create an interactive HTML report from a tabular metagenomic report
    '''

    krona_tool = krona.Krona()
    if sample_name is not None:
        dataset_names = list(sample_name for fn in inReports)
    else:
        dataset_names = list(os.path.basename(fn) for fn in inReports)

    with file.tmp_dir() as tmp_dir:
        to_import = []

        if inputType == 'tsv':
            for inReport in inReports:
                if inReport.endswith('.gz'):
                    tmp_tsv = file.mkstempfname('.tsv', directory=tmp_dir)
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

__commands__.append(('krona', parser_krona))


def parser_metagenomic_report_merge(parser=argparse.ArgumentParser()):
    parser.add_argument(
        "metagenomic_reports",
        help="Input metagenomic reports with the query ID and taxon ID in the 2nd and 3rd columns (Kraken format)",
        nargs='+',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        "out_krona_input",
        help="Output metagenomic report suitable for Krona input."
    )
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, metagenomic_report_merge, split_args=True)
    return parser
def metagenomic_report_merge(metagenomic_reports, out_krona_input):
    '''
        Merge multiple metagenomic reports into a single metagenomic report suitable for Krona input.
    '''
    with file.open_or_gzopen(out_krona_input, "wt") as outf:
        output_writer = csv.writer(outf, delimiter='\t', lineterminator='\n')
        for metag_file in metagenomic_reports:
            with file.open_or_gzopen(metag_file.name, "rt") as inf:
                file_reader = csv.reader(inf, delimiter='\t')
                for row in file_reader:
                    output_writer.writerow([f for f in row])
__commands__.append(('report_merge', parser_metagenomic_report_merge))


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
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, filter_bam_to_taxa, split_args=True)
    return parser

def filter_bam_to_taxa(in_bam, read_IDs_to_tax_IDs, out_bam,
                       nodes_dmp, names_dmp,
                       exclude=False,
                       tax_names=None, tax_ids=None,
                       omit_children=False,
                       read_id_col=1, tax_id_col=2,
                       out_count=None):
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

    def _matching_read_ids():
        """Generator that yields read IDs matching the specified taxa."""
        for row in file.read_tabfile(read_IDs_to_tax_IDs):
            assert tax_id_col<len(row), "tax_id_col does not appear to be in range for number of columns present in mapping file"
            assert read_id_col<len(row), "read_id_col does not appear to be in range for number of columns present in mapping file"
            read_id = row[read_id_col]
            read_tax_id = int(row[tax_id_col])

            # transform read ID to take read pairs into account
            read_id_match = re.match(paired_read_base_pattern, read_id)
            if (read_id_match and read_tax_id in tax_ids_to_include):
                yield read_id_match.group(1)

    # Stream matching read IDs directly into ReadIdStore
    with file.tmp_dir(suffix='_filter_taxa') as tmpdir:
        db_path = os.path.join(tmpdir, 'read_ids.db')
        with read_utils.ReadIdStore(db_path) as store:
            store.extend(_matching_read_ids())

            log.info("matched {} reads".format(len(store)))

            # report count if desired
            if out_count:
                with open(out_count, 'wt') as outf:
                    outf.write("{}\n".format(len(store)))

            # filter the input bam (include=True keeps matching, include=False removes matching)
            store.filter_bam_by_ids(in_bam, out_bam, include=not exclude)
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
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, taxlevel_summary, split_args=True)
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


        with file.open_or_gzopen(f, 'rt') as inf:
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


def parser_kraken_taxlevel_plurality(parser=argparse.ArgumentParser()):
    parser.add_argument('summary_file', help='input Kraken-format summary text file with tab-delimited taxonomic levels.')
    parser.add_argument('tax_heading', help='The taxonomic heading to analyze.')
    parser.add_argument('out_report', help='tab-delimited output file.')
    parser.add_argument('--min_reads', type=int, dest="min_reads", help='Only include hits with more than min_reads (default: %(default)s)', default=1)
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, taxlevel_plurality, split_args=True)
    return parser

def taxlevel_plurality(summary_file, tax_heading, out_report, min_reads):
    """
        Identifies the most abundant taxon (of any rank) contributing to a node of interest in the taxonomic tree.
        It is intended to highlight the primary contributor of taxonomic signal within a taxonomic category of interest,
        for example, the most abundant virus among all viruses.
    """

    keeper_rows = []
    with file.open_or_gzopen(summary_file, 'rt') as inf:
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
    with file.open_or_gzopen(out_report, 'wt') as outf:
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
    parser.add_argument('in_bam', help='Input unaligned reads, BAM format.')
    parser.add_argument('--index', help='kb index file.')
    parser.add_argument('--t2g', help='Transcript to gene mapping file.')
    parser.add_argument('--out_dir', dest='out_dir', help='Output directory (default: kb_out)', default='kb_out')
    parser.add_argument('--protein', action='store_true', help='True if sequence contains amino acids (default: False).')
    parser.add_argument('--targets', help='Comma-separated list of target sequences to extract from input sequences.', default=None)
    parser.add_argument('--h5ad', help='Path to the output h5ad file. Can pull IDs to extract from this file.', default=None)
    parser.add_argument('--threshold', type=int, help='Minimum read count threshold for a target to be extracted (only used when extractin IDs from h5ad; default: %(default)s)', default=1)
    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kb_extract, split_args=True)
    return parser
def kb_extract(in_bam, index, t2g, targets, protein=False, out_dir=None, h5ad=None, threads=None, threshold=None):
    """Runs kb extract on the input BAM file.

    Args:
        in_bam (str): Input BAM file.
        index (str): Path to the kb index file.
        t2g (str): Path to the transcript-to-gene mapping file.
        targets (str): Comma-separated list of target sequences to extract.
        protein (bool): True if sequence contains amino acids. Defaults to False.
        out_dir (str): Output directory. Defaults to None.
        h5ad (str): Path to the output h5ad file. Can pull IDs to extract from this file. Defaults to None.
        threshold (int, optional): Minimum read count threshold for a target to be extracted. Defaults to 1.
        threads (int, optional): Number of threads to use. Defaults to None.
    """
    assert out_dir, ('Output directory must be specified.')

    kb_tool = kb.kb()
    
    target_ids = targets.split(',') if targets else []
    if not target_ids or len(target_ids) == 0:
        # TODO: This extraction method expects only to have a single row h5ad (i.e. 1 sample). This should be handled more robustly.
        log.warning('No targets specified for extraction. Trying to extract IDs from h5ad.')
        target_ids = kb_tool.extract_hit_ids_from_h5ad(h5ad, threshold=threshold)
        log.info("Target IDs extracted from h5ad: {}".format(target_ids))
        if len(target_ids) == 0:
            raise ValueError('No targets specified for extraction and no IDs found in h5ad.')

    kb_tool.extract(
        in_bam=in_bam,
        index_file=index,
        target_ids=target_ids,
        out_dir=out_dir,
        t2g_file=t2g,
        protein=protein,
        num_threads=threads
    )
__commands__.append(('kb_extract', parser_kb_extract))

def parser_kb_top_taxa(parser=argparse.ArgumentParser()):
    parser.add_argument('counts_tar', help='Input kb count tarball (tar.zst format).')
    parser.add_argument('--id-to-tax-map', dest='id_to_tax_map', help='ID to taxonomy mapping file (CSV format).')
    parser.add_argument('--target-taxon', dest='target_taxon', default='Viruses', help='Target taxonomic category to analyze (default: Viruses).')
    parser.add_argument('out_report', help='Tab-delimited output file.')
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kb_top_taxa, split_args=True)
    return parser

def kb_top_taxa(counts_tar, out_report, id_to_tax_map=None, target_taxon='Viruses'):
    """Identifies the most abundant taxon (of any rank) contributing to a taxa node of interest in kb count output.

    It is intended to highlight the primary contributor of taxonomic signal within a taxonomic category of interest,
    for example, the most abundant virus among all viruses.

    Args:
        counts_tar (str): Path to the input kb count tarball (tar.zst format).
        out_report (str): Path to the output report file.
        id_to_tax_map (str, optional): Path to the ID to taxonomy mapping file (CSV format).
        target_taxon (str): The taxonomic category to analyze (default: 'Viruses').
    """
    kb_tool = kb.kb()

    # Extract and read h5ad file from tarball
    with file.tmp_dir() as tmp_dir:
        file.extract_tarball(counts_tar, tmp_dir)
        h5ad_files = glob.glob(os.path.join(tmp_dir, "counts_unfiltered", '*.h5ad'))

        assert len(h5ad_files) == 1, "Expected exactly one .h5ad file in the counts tarball, found {}".format(len(h5ad_files))
        h5ad_file = h5ad_files[0]

        # Use kb helper to parse h5ad file
        gene_counts = kb_tool.parse_h5ad_counts(h5ad_file)

        # Build results data with taxonomy information
        results_data = []

        if id_to_tax_map:
            # Load taxonomy mapping
            tax_map_df = pd.read_csv(id_to_tax_map)
            tax_map_df = tax_map_df.replace('.', np.nan)

            # Standardize column name: handle both 'ID' and 'id'
            if 'ID' in tax_map_df.columns:
                tax_map_df = tax_map_df.rename(columns={'ID': 'id'})
            elif 'id' not in tax_map_df.columns:
                # If neither exists, use the first column as 'id'
                tax_map_df = tax_map_df.rename(columns={tax_map_df.columns[0]: 'id'})

            # Get the farthest (most specific) non-empty taxonomy for each gene
            # Define taxonomy rank columns in order from broad to specific
            taxonomy_columns = [col for col in ['phylum', 'class', 'order', 'family', 'genus', 'species']
                              if col in tax_map_df.columns]

            def get_farthest_taxa(row):
                # Only consider taxonomy columns, in reverse order (most specific first)
                for col in reversed(taxonomy_columns):
                    if pd.notna(row[col]) and row[col] != '':
                        return row[col]
                return 'Unknown'

            tax_map_df['farthest_taxa'] = tax_map_df.apply(get_farthest_taxa, axis=1)

            # Get the root taxonomy category
            if 'root' in tax_map_df.columns:
                tax_map_df['tax_category'] = tax_map_df['root']
            elif 'phylum' in tax_map_df.columns:
                tax_map_df['tax_category'] = tax_map_df['phylum']
            else:
                # Assume first column after 'id' is the root category
                tax_map_df['tax_category'] = tax_map_df.iloc[:, 1]

            # OPTIMIZATION: Set 'id' as index for O(1) lookups instead of O(n) filtering
            tax_map_df = tax_map_df.set_index('id')

            # Filter by target taxon - optimized with index-based lookup
            for gene_id, count in gene_counts:
                # Skip genes with 0 counts
                if count == 0:
                    continue
                if gene_id in tax_map_df.index:
                    tax_info = tax_map_df.loc[gene_id]
                    tax_category = tax_info['tax_category']
                    if tax_category == target_taxon:
                        # Use species name if available, otherwise use farthest_taxa
                        species_name = tax_info.get('species', tax_info['farthest_taxa']) if 'species' in tax_info.index else tax_info['farthest_taxa']
                        if pd.isna(species_name) or species_name == '':
                            species_name = tax_info['farthest_taxa']

                        results_data.append({
                            'palmdb_id': gene_id,
                            'sotu_id': species_name,
                            'count': int(count),
                            'tax_category': tax_category,
                            'farthest_taxa': tax_info['farthest_taxa']
                        })
        else:
            # No taxonomy mapping - report all genes
            for gene_id, count in gene_counts:
                # Skip genes with 0 counts
                if count == 0:
                    continue
                results_data.append({
                    'palmdb_id': gene_id,
                    'sotu_id': gene_id,
                    'count': int(count),
                    'tax_category': 'Unknown',
                    'farthest_taxa': gene_id
                })

        # Create output
        if not results_data:
            out = [{'focal_taxon_name': target_taxon, 'focal_taxon_count': 0, 'pct_of_focal': 0.0,
                    'palmdb_id': '', 'hit_id': '', 'hit_lowest_taxa_name': '', 'hit_reads': 0}]
        else:
            results_df = pd.DataFrame(results_data)
            total_focal_count = results_df['count'].sum()
            results_df = results_df.sort_values(by='count', ascending=False).reset_index(drop=True)

            # OPTIMIZATION: Use vectorized operations instead of iterrows()
            results_df['focal_taxon_name'] = target_taxon
            results_df['focal_taxon_count'] = total_focal_count
            results_df['pct_of_focal'] = 100.0 * results_df['count'] / total_focal_count

            # Create output dict
            out = results_df.rename(columns={
                'sotu_id': 'hit_id',
                'farthest_taxa': 'hit_lowest_taxa_name',
                'count': 'hit_reads'
            })[['focal_taxon_name', 'focal_taxon_count', 'palmdb_id', 'hit_id', 'hit_lowest_taxa_name', 'hit_reads', 'pct_of_focal']].to_dict('records')

    # Write output
    with file.open_or_gzopen(out_report, 'wt') as outf:
        header = ('focal_taxon_name', 'focal_taxon_count', 'palmdb_id', 'hit_id', 'hit_lowest_taxa_name', 'hit_reads', 'pct_of_focal')
        writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(out)

__commands__.append(('kb_top_taxa', parser_kb_top_taxa))

def parser_kb_merge_h5ads(parser=argparse.ArgumentParser()):
    parser.add_argument('in_count_tars', nargs='+', help='Input kb count tarballs to merge (tar.zst format).')
    parser.add_argument('--out-h5ad', dest='out_h5ad', help='Output merged h5ad file.')
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kb_merge_h5ads, split_args=True)
    return parser
def kb_merge_h5ads(in_count_tars, out_h5ad, tmp_dir=None):
    '''
    Merge multiple kb count output tarballs into a single h5ad file with sample metadata.

    Extracts h5ad files from counts_unfiltered folder and adds sample names from matrix.cells.

    Args:
        in_count_tars (list): List of input kb count tarballs (tar.zst format).
        out_h5ad (str): Path to the output h5ad file.
        tmp_dir (str, optional): Temporary directory for extraction.
    '''
    assert out_h5ad, ('Output h5ad file must be specified.')
    kb_tool = kb.kb()
    kb_tool.merge_h5ads(
        in_count_tars=in_count_tars,
        out_h5ad=out_h5ad,
        tmp_dir_parent=tmp_dir
    )

__commands__.append(('kb_merge_h5ads', parser_kb_merge_h5ads))

def parser_krona_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Krona taxonomy database output directory.')
    parser.add_argument('--taxdump_tar_gz', help='NCBI taxdump.tar.gz file', default=None)
    parser.add_argument('--get_accessions', action='store_true', help='Fetch NCBI accession to taxid mappings. This is not required for processing kraken1/2/uniq hits, only for BLAST hits, and adds a significant amount of time and database space (default false).')
    cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, krona_build, split_args=True)
    return parser
def krona_build(db, taxdump_tar_gz=None, get_accessions=False):
    '''
    Builds a Krona taxonomy database
    '''
    krona.Krona().build_db(
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
    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kraken2_build, split_args=True)
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

    kraken_tool = kraken2.Kraken2()
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
    parser.add_argument('--index', help='kb output index file.')
    parser.add_argument('--workflow', choices=['standard', 'nac', 'kite', 'custom'],
                        default='standard', help='Type of index to create (default: %(default)s).')
    parser.add_argument('--kmer_len', type=int, help='k-mer length (default: 31).')
    parser.add_argument('--protein', action='store_true', help='True if sequence contains amino acids(default: False).')
    cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    cmd.attach_main(parser, kb_build, split_args=True)
    return parser
def kb_build(ref_fasta, index, workflow='standard', kmer_len=31, protein=False, threads=None):
    '''
    Builds a kb index from a reference fasta file.

    Args:
        ref_fasta (str): Path to the reference sequence fasta file.
        index (str): Path to the output kb index file.
        workflow (str): Type of index to create. Options are 'standard', 'nac', 'kite', 'custom'.
        kmer_len (int): k-mer length (default: 31).
        protein (bool): True if sequence contains amino acids (default: False).
        threads (int): Number of threads to use (default: None).
    '''
    kb_tool = kb.kb()
    kb_tool.build(ref_fasta,
                        index=index,
                        workflow=workflow,
                        kmer_len=kmer_len,
                        protein=protein,
                        num_threads=threads)
__commands__.append(('kb_build', parser_kb_build))


def full_parser():
    return cmd.make_parser(__commands__, __doc__)


def main():
    cmd.main_argparse(__commands__, __doc__)


if __name__ == '__main__':
    main()

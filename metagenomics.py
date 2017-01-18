#!/usr/bin/env python
''' This script contains a number of utilities for metagenomic analyses.
'''
from __future__ import print_function
from __future__ import division

__author__ = "yesimon@broadinstitute.org"

import argparse
import collections
import csv
import gzip
import itertools
import logging
import os.path
from os.path import join
import operator
import queue
import shutil
import sys

import pysam

import util.cmd
import util.file
import util.misc
import tools.bwa
import tools.diamond
import tools.kraken
import tools.krona
import tools.picard

__commands__ = []

log = logging.getLogger(__name__)


class TaxIdError(ValueError):
    '''Taxonomy ID couldn't be determined.'''


class TaxonomyDb(object):

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
            gis_paths = [join(tax_dir, 'gi_taxid_nucl.dmp'), join(tax_dir, 'gi_taxid_prot.dmp')]
            nodes_path = join(tax_dir, 'nodes.dmp')
            names_path = join(tax_dir, 'names.dmp')
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
        with open(dmp_path) as f:
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
        for line in open(names_db):
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
        with open(nodes_db) as f:
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


BlastRecord = collections.namedtuple(
    'BlastRecord', [
        'query_id', 'subject_id', 'percent_identity', 'aln_length', 'mismatch_count', 'gap_open_count', 'query_start',
        'query_end', 'subject_start', 'subject_end', 'e_val', 'bit_score'
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
        yield BlastRecord(*parts)


def paired_query_id(record):
    '''Replace paired suffixes in query ids.'''
    suffixes = ('/1', '/2')
    for suffix in suffixes:
        if record.query_id.endswith(suffix):
            rec_list = list(record)
            rec_list[0] = record.query_id[:-len(suffix)]
            return BlastRecord(*rec_list)
    return record


def translate_gi_to_tax_id(db, record):
    '''Replace gi headers in subject ids to int taxonomy ids.'''
    gi = int(record.subject_id.split('|')[1])
    tax_id = db.gis[gi]
    rec_list = list(record)
    rec_list[1] = tax_id
    return BlastRecord(*rec_list)


def extract_tax_id(sam1):
    '''Replace gi headers in subject ids to int taxonomy ids.'''
    parts = sam1.reference_name.split('|')
    if parts[0] == 'taxid':
        return int(parts[1])
    else:
        raise TaxIdError(parts)


def sam_lca(db, sam_file, output=None, top_percent=10, unique_only=True):
    ''' Calculate the LCA taxonomy id for multi-mapped reads in a samfile.

    Assumes the sam is sorted by query name. Writes tsv output: query_id \t tax_id.

    Args:
      db: (TaxonomyDb) Taxonomy db.
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
                tax_id = process_sam_hits(db, mapped_segs, top_percent)
                if tax_id is None:
                    log.warn('Query: {} has no valid taxonomy paths.'.format(query_name))
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


def blast_lca(db,
              m8_file,
              output,
              paired=False,
              min_bit_score=50,
              max_expected_value=0.01,
              top_percent=10,):
    '''Calculate the LCA taxonomy id for groups of blast hits.

    Writes tsv output: query_id \t tax_id

    Args:
      db: (TaxonomyDb) Taxonomy db.
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
        tax_id = process_blast_hits(db, blast_group, top_percent)
        query_id = blast_group[0].query_id
        if not tax_id:
            log.debug('Query: {} has no valid taxonomy paths.'.format(query_id))
        classified = 'C' if tax_id else 'U'
        output.write('{}\t{}\t{}\n'.format(classified, query_id, tax_id))


def process_sam_hits(db, sam_hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
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
    return coverage_lca(tax_ids, db.parents)


def process_blast_hits(db, blast_hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      blast_hits: []BlastRecord groups of hits.
      top_percent: (float) Only consider hits within this percent of top bit score.

    Return:
      (int) Tax id of LCA.
    '''
    hits = (translate_gi_to_tax_id(db, hit) for hit in blast_hits)
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
        tax_ids = [hit.subject_id for hit in valid_hits]
        return coverage_lca(tax_ids, db.parents)


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
                log.warn('Parent for query id: {} missing'.format(query_id))
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
    else:
        return "-"


def taxa_hits_from_tsv(f, tax_id_column=3):
    '''Return a counter of hits from tsv.'''
    c = collections.Counter()
    for row in csv.reader(f, delimiter='\t'):
        tax_id = int(row[tax_id_column - 1])
        c[tax_id] += 1
    return c


def kraken_dfs_report(db, taxa_hits):
    '''Return a kraken compatible DFS report of taxa hits.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      taxa_hits: (collections.Counter) # of hits per tax id.

    Return:
      []str lines of the report
    '''

    db.children = parents_to_children(db.parents)
    total_hits = sum(taxa_hits.values())
    lines = []
    kraken_dfs(db, lines, taxa_hits, total_hits, 1, 0)
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


def kraken_dfs(db, lines, taxa_hits, total_hits, taxid, level):
    '''Recursively do DFS for number of hits per taxa.'''
    cum_hits = num_hits = taxa_hits.get(taxid, 0)
    for child_taxid in db.children[taxid]:
        cum_hits += kraken_dfs(db, lines, taxa_hits, total_hits, child_taxid, level + 1)
    percent_covered = '%.2f' % (cum_hits / total_hits * 100)
    rank = rank_code(db.ranks[taxid])
    name = db.names[taxid]
    if cum_hits > 0:
        lines.append('\t'.join([percent_covered, str(cum_hits), str(num_hits), rank, str(taxid), '  ' * level + name]))
    return cum_hits


def kraken(inBam, db, outReport=None, outReads=None, filterThreshold=None, numThreads=1):
    '''
        Classify reads by taxon using Kraken
    '''

    assert outReads or outReport, ('Either --outReads or --outReport must be specified.')
    kraken_tool = tools.kraken.Kraken()

    # kraken classify
    tmp_reads = util.file.mkstempfname('.kraken')
    kraken_tool.classify(inBam, db, tmp_reads, numThreads=numThreads)

    # kraken filter
    if filterThreshold:
        tmp_filtered_reads = util.file.mkstempfname('.filtered-kraken')
        kraken_tool.filter(tmp_reads, db, tmp_filtered_reads, filterThreshold)
        os.unlink(tmp_reads)
    else:
        tmp_filtered_reads = tmp_reads

    # copy outReads
    if outReads:
        with open(tmp_filtered_reads, 'rb') as f_in:
            with util.file.open_or_gzopen(outReads, 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # kraken report
    if outReport:
        kraken_tool.report(tmp_filtered_reads, db, outReport)

    os.unlink(tmp_filtered_reads)


def parser_kraken(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('--outReport', help='Kraken report output file.')
    parser.add_argument('--outReads', help='Kraken per read output file.')
    parser.add_argument(
        '--filterThreshold', default=0.05, type=float, help='Kraken filter threshold (default %(default)s)'
    )
    parser.add_argument('--numThreads', type=int, default=1, help='Number of threads to run. (default %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken, split_args=True)
    return parser


def krona(inTsv, db, outHtml, queryColumn=None, taxidColumn=None, scoreColumn=None, noHits=None, noRank=None):
    '''
        Create an interactive HTML report from a tabular metagenomic report
    '''

    krona_tool = tools.krona.Krona()
    if inTsv.endswith('.gz'):
        tmp_tsv = util.file.mkstempfname('.tsv')
        with gzip.open(inTsv, 'rb') as f_in:
            with open(tmp_tsv, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                to_import = [tmp_tsv]
    else:
        to_import = [inTsv]

    krona_tool.import_taxonomy(
        db,
        to_import,
        outHtml,
        query_column=queryColumn,
        taxid_column=taxidColumn,
        score_column=scoreColumn,
        no_hits=noHits,
        no_rank=noRank
    )


def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inTsv', help='Input tab delimited file.')
    parser.add_argument('db', help='Krona taxonomy database directory.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--queryColumn', help='Column of query id. (default %(default)s)', type=int, default=2)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id. (default %(default)s)', type=int, default=3)
    parser.add_argument('--scoreColumn', help='Column of score. (default %(default)s)', type=int)
    parser.add_argument('--noHits', help='Include wedge for no hits.', action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.', action='store_true')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser


def diamond(inBam, db, taxDb, outReport, outM8=None, outLca=None, numThreads=1):
    '''
        Classify reads by the taxon of the Lowest Common Ancestor (LCA)
    '''
    tmp_fastq = util.file.mkstempfname('.fastq')
    tmp_fastq2 = util.file.mkstempfname('.fastq')
    # do not convert this to samtools bam2fq unless we can figure out how to replicate
    # the clipping functionality of Picard SamToFastq
    picard = tools.picard.SamToFastqTool()
    picard_opts = {
        'CLIPPING_ATTRIBUTE': tools.picard.SamToFastqTool.illumina_clipping_attribute,
        'CLIPPING_ACTION': 'X'
    }
    picard.execute(
        inBam,
        tmp_fastq,
        tmp_fastq2,
        picardOptions=tools.picard.PicardTools.dict_to_picard_opts(picard_opts),
        JVMmemory=picard.jvmMemDefault
    )

    diamond_tool = tools.diamond.Diamond()
    diamond_tool.install()
    tmp_alignment = util.file.mkstempfname('.daa')
    tmp_m8 = util.file.mkstempfname('.diamond.m8')
    diamond_tool.blastx(db, [tmp_fastq, tmp_fastq2], tmp_alignment, options={'--threads': numThreads})
    diamond_tool.view(tmp_alignment, tmp_m8, options={'--threads': numThreads})

    if outM8:
        with open(tmp_m8, 'rb') as f_in:
            with gzip.open(outM8, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    tax_db = TaxonomyDb(tax_dir=taxDb, load_names=True, load_nodes=True, load_gis=True)
    tmp_lca_tsv = util.file.mkstempfname('.tsv')
    with open(tmp_m8) as m8, open(tmp_lca_tsv, 'w') as lca:
        blast_lca(tax_db, m8, lca, paired=True, min_bit_score=50)

    if outLca:
        with open(tmp_lca_tsv, 'rb') as f_in:
            with gzip.open(outLca, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    with open(tmp_lca_tsv) as f:
        hits = taxa_hits_from_tsv(f)
    with open(outReport, 'w') as f:
        for line in kraken_dfs_report(tax_db, hits):
            print(line, file=f)


def parser_diamond(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Diamond database directory.')
    parser.add_argument('taxDb', help='Taxonomy database directory.')
    parser.add_argument('outReport', help='Output taxonomy report.')
    parser.add_argument('--outM8', help='Blast m8 formatted output file.')
    parser.add_argument('--outLca', help='Output LCA assignments for each read.')
    parser.add_argument('--numThreads', default=1, help='Number of threads (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, diamond, split_args=True)
    return parser


def align_rna_metagenomics(
    inBam,
    db,
    taxDb,
    outReport,
    dupeReport=None,
    outBam=None,
    dupeLca=None,
    outLca=None,
    sensitive=None,
    JVMmemory=None,
    numThreads=None,
    picardOptions=None,
    min_score_to_output=None,
):
    '''
        Align to metagenomics bwa index, mark duplicates, and generate LCA report
    '''
    picardOptions = picardOptions if picardOptions else []

    bwa = tools.bwa.Bwa()
    samtools = tools.samtools.SamtoolsTool()

    bwa_opts = ['-a']
    if sensitive:
        bwa_opts += '-k 12 -A 1 -B 1 -O 1 -E 1'.split()

    map_threshold = min_score_to_output or 30

    aln_bam = util.file.mkstempfname('.bam')
    bwa.mem(inBam, db, aln_bam, options=bwa_opts, min_qual=map_threshold)

    tax_db = TaxonomyDb(tax_dir=taxDb, load_names=True, load_nodes=True)

    if dupeReport:
        aln_bam_sorted = util.file.mkstempfname('.align_namesorted.bam')
        samtools.sort(aln_bam, aln_bam_sorted, args=['-n'], threads=numThreads)
        sam_lca_report(tax_db, aln_bam_sorted, outReport=dupeReport, outLca=dupeLca, unique_only=False)
        os.unlink(aln_bam_sorted)

    aln_bam_deduped = outBam if outBam else util.file.mkstempfname('.align_deduped.bam')
    opts = list(picardOptions)
    dupe_removal_out_metrics = util.file.mkstempfname('.metrics')
    pic = tools.picard.MarkDuplicatesTool()
    pic.execute([aln_bam], aln_bam_deduped, dupe_removal_out_metrics, picardOptions=opts, JVMmemory=JVMmemory)
    os.unlink(aln_bam)

    sam_lca_report(tax_db, aln_bam_deduped, outReport=outReport, outLca=outLca)

    if not outBam:
        os.unlink(aln_bam_deduped)


def sam_lca_report(tax_db, bam_aligned, outReport, outLca=None, unique_only=None):
    lca_tsv = outLca

    if outLca:
        lca_tsv = outLca
    else:
        lca_tsv = util.file.mkstempfname('.tsv')

    with util.file.open_or_gzopen(lca_tsv, 'wt') as lca:
        hits = sam_lca(tax_db, bam_aligned, lca, top_percent=10, unique_only=unique_only)

    with open(outReport, 'w') as f:

        for line in kraken_dfs_report(tax_db, hits):
            print(line, file=f)


def parser_align_rna_metagenomics(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Bwa index prefix.')
    parser.add_argument('taxDb', help='Taxonomy database directory.')
    parser.add_argument('outReport', help='Output taxonomy report.')
    parser.add_argument('--dupeReport', help='Generate report including duplicates.')
    parser.add_argument(
        '--sensitive',
        dest='sensitive',
        action="store_true",
        help='Use sensitive instead of default BWA mem options.'
    )
    parser.add_argument('--outBam', help='Output aligned, indexed BAM file. Default is to write to temp.')
    parser.add_argument('--outLca', help='Output LCA assignments for each read.')
    parser.add_argument('--dupeLca', help='Output LCA assignments for each read including duplicates.')
    parser.add_argument('--numThreads', default=1, help='Number of threads (default: %(default)s)')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.PicardTools.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, align_rna_metagenomics, split_args=True)
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

        kraken_tool = tools.kraken.Kraken()
        kraken_tool.report(tmp_metag_combined_txt, kraken_db.name, out_kraken_summary)


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


__commands__.append(('kraken', parser_kraken))
__commands__.append(('diamond', parser_diamond))
__commands__.append(('krona', parser_krona))
__commands__.append(('align_rna', parser_align_rna_metagenomics))
__commands__.append(('report_merge', parser_metagenomic_report_merge))


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

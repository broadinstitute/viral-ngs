"""
NCBI Taxonomy database utilities.

This module provides classes and functions for working with NCBI taxonomy data,
including taxonomy ID lookups, LCA calculations, and tree traversal.
"""

import collections
import itertools
import logging
import operator
import os

import pysam

from ..core import file as util_file

log = logging.getLogger(__name__)


class TaxIdError(ValueError):
    '''Taxonomy ID couldn't be determined.'''


def maybe_compressed(fn):
    """Check for file existence, returning .gz version if original not found."""
    fn_gz = fn + '.gz'
    if os.path.exists(fn):
        return fn
    elif os.path.exists(fn_gz):
        return fn_gz
    else:
        raise FileNotFoundError(fn)


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
    elif rank == "unclassified":
        return "U"
    else:
        return "-"


class TaxonomyDb(object):
    """
    NCBI Taxonomy database handler.

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
        with util_file.open_or_gzopen(dmp_path) as f:
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
        for line in util_file.open_or_gzopen(names_db):
            parts = line.strip().split('|')
            taxid = int(parts[0])
            name = parts[1].strip()
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
        with util_file.open_or_gzopen(nodes_db) as f:
            for line in f:
                parts = line.strip().split('|')
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
                rank = parts[2].strip()
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
            lca_tsv = util_file.mkstempfname('.tsv')
        with util_file.open_or_gzopen(lca_tsv, 'wt') as lca:
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

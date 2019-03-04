#!/usr/bin/env python
''' This script contains a number of utilities for metagenomic analyses.
'''
from __future__ import print_function
from __future__ import division

__author__ = "yesimon@broadinstitute.org"

import argparse
import codecs
import collections
import csv
import gzip
import io
import itertools
import logging
import os.path
from os.path import join
import operator
import queue
import re
import shutil
import sys
import subprocess
import tempfile
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import ncbitax
import ncbitax.subset
import pysam

import util.cmd
import util.file
import util.misc
import tools.bwa
import tools.diamond
import tools.kraken
import tools.krona
import tools.picard
import tools.samtools
from util.file import open_or_gzopen, maybe_compressed

__commands__ = []

log = logging.getLogger(__name__)


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


def translate_gi_to_tax_id(db, record):
    '''Replace gi headers in subject ids to int taxonomy ids.'''
    gi = int(record.subject_id.split('|')[1])
    tax_id = db.gis[gi]
    rec_list = list(record)
    rec_list[1] = tax_id
    return BlastRecord(*rec_list)


def blast_m8_taxids(record):
    return [int(record.subject_id)]


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
    return db.coverage_lca(tax_ids)


def process_blast_hits(db, hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      hits: []BlastRecord groups of hits.
      top_percent: (float) Only consider hits within this percent of top bit score.

    Return:
      (int) Tax id of LCA.
    '''
    hits = (translate_gi_to_tax_id(db, hit) for hit in hits)

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
        return db.coverage_lca(tax_ids)


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


def parser_kraken(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database directory.')
    parser.add_argument('inBams', nargs='+', help='Input unaligned reads, BAM format.')
    parser.add_argument('--outReports', nargs='+', help='Kraken summary report output file. Multiple filenames space separated.')
    parser.add_argument('--outReads', nargs='+', help='Kraken per read classification output file. Multiple filenames space separated.')
    parser.add_argument('--lockMemory', action='store_true', default=False, help='Lock kraken database in RAM. Requires high ulimit -l.')
    parser.add_argument(
        '--filterThreshold', default=0.05, type=float, help='Kraken filter threshold (default %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken, split_args=True)
    return parser
def kraken(db, inBams, outReports=None, outReads=None, lockMemory=False, filterThreshold=None, threads=None):
    '''
        Classify reads by taxon using Kraken
    '''

    assert outReads or outReports, ('Either --outReads or --outReport must be specified.')
    kraken_tool = tools.kraken.Kraken()
    kraken_tool.pipeline(db, inBams, outReports=outReports, outReads=outReads, lockMemory=lockMemory,
                         filterThreshold=filterThreshold, numThreads=threads)
__commands__.append(('kraken', parser_kraken))



def parser_krona(parser=argparse.ArgumentParser()):
    parser.add_argument('inTsv', help='Input tab delimited file.')
    parser.add_argument('db', help='Krona taxonomy database directory.')
    parser.add_argument('outHtml', help='Output html report.')
    parser.add_argument('--queryColumn', help='Column of query id. (default %(default)s)', type=int, default=2)
    parser.add_argument('--taxidColumn', help='Column of taxonomy id. (default %(default)s)', type=int, default=3)
    parser.add_argument('--scoreColumn', help='Column of score. (default %(default)s)', type=int, default=None)
    parser.add_argument('--magnitudeColumn', help='Column of magnitude. (default %(default)s)', type=int, default=None)
    parser.add_argument('--noHits', help='Include wedge for no hits.', action='store_true')
    parser.add_argument('--noRank', help='Include no rank assignments.', action='store_true')
    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, krona, split_args=True)
    return parser
def krona(inTsv, db, outHtml, queryColumn=None, taxidColumn=None, scoreColumn=None, magnitudeColumn=None, noHits=None, noRank=None):
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
    root_name = os.path.basename(inTsv)

    krona_tool.import_taxonomy(
        db,
        to_import,
        outHtml,
        query_column=queryColumn,
        taxid_column=taxidColumn,
        score_column=scoreColumn,
        magnitude_column=magnitudeColumn,
        root_name=root_name,
        no_hits=noHits,
        no_rank=noRank
    )

    if inTsv.endswith('.gz'):
        # Cleanup tmp .tsv files
        for tmp_tsv in to_import:
            os.unlink(tmp_tsv)
__commands__.append(('krona', parser_krona))


def parser_diamond(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam', help='Input unaligned reads, BAM format.')
    parser.add_argument('db', help='Diamond database directory.')
    parser.add_argument('taxDb', help='Taxonomy database directory.')
    parser.add_argument('outReport', help='Output taxonomy report.')
    parser.add_argument('--outReads', help='Output LCA assignments for each read.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, diamond, split_args=True)
    return parser
def diamond(inBam, db, taxDb, outReport, outReads=None, threads=None):
    '''
        Classify reads by the taxon of the Lowest Common Ancestor (LCA)
    '''
    # do not convert this to samtools bam2fq unless we can figure out how to replicate
    # the clipping functionality of Picard SamToFastq
    picard = tools.picard.SamToFastqTool()
    with util.file.fifo(2) as (fastq_pipe, diamond_pipe):
        s2fq = picard.execute(
            inBam,
            fastq_pipe,
            interleave=True,
            illuminaClipping=True,
            JVMmemory=picard.jvmMemDefault,
            background=True,
        )

        diamond_tool = tools.diamond.Diamond()
        taxonmap = join(taxDb, 'accession2taxid', 'prot.accession2taxid.gz')
        taxonnodes = join(taxDb, 'nodes.dmp')
        rutils = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'read_utils.py')
        cmd = '{read_utils} join_paired_fastq --outFormat fasta /dev/stdout {fastq}'.format(
            read_utils=rutils, fastq=fastq_pipe)

        cmd += ' | {} blastx --out {output} --outfmt 102 --sallseqid'.format(
            diamond_tool.install_and_get_path(), output=diamond_pipe)
        cmd += ' --threads {threads}'.format(threads=util.misc.sanitize_thread_count(threads))
        cmd += ' --db {db} --taxonmap {taxonmap} --taxonnodes {taxonnodes}'.format(
            db=db,
            taxonmap=taxonmap,
            taxonnodes=taxonnodes)

        if outReads is not None:
            # Interstitial save of stdout to output file
            cmd += ' | tee >(pigz --best > {out})'.format(out=outReads)

        diamond_ps = subprocess.Popen(cmd, shell=True, executable='/bin/bash')

        tax_db = ncbitax.TaxonomyDb(tax_dir=taxDb, load_names=True, load_nodes=True)

        with open(diamond_pipe) as lca_p:
            hits = ncbitax.taxa_hits_from_tsv(lca_p)
            with open(outReport, 'w') as f:
                for line in tax_db.kraken_dfs_report(hits):
                    print(line, file=f)

        s2fq.wait()
        assert s2fq.returncode == 0
        diamond_ps.wait()
        assert diamond_ps.returncode == 0
__commands__.append(('diamond', parser_diamond))


def parser_diamond_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta', help='Input sequences, FASTA format, optionally gzip compressed.')
    parser.add_argument('db', help='Diamond database file.')
    parser.add_argument('taxDb', help='Taxonomy database directory.')
    parser.add_argument('outFasta', help='Output sequences, same as inFasta, with taxid|###| prepended to each sequence identifier.')
    parser.add_argument('--memLimitGb', type=float, default=None, help='approximate memory usage in GB')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, diamond_fasta, split_args=True)
    return parser
def diamond_fasta(inFasta, db, taxDb, outFasta, threads=None, memLimitGb=None):
    '''
        Classify fasta sequences by the taxon of the Lowest Common Ancestor (LCA)
    '''

    with util.file.tmp_dir() as tmp_dir:
        # run diamond blastx on fasta sequences
        cmd = [tools.diamond.Diamond().install_and_get_path(), 'blastx',
            '-q', inFasta,
            '--db', db,
            '--outfmt', '102', # tsv: query name, taxid of LCA, e-value
            '--salltitles',# to recover the entire fasta sequence name
            '--sensitive', # this is necessary for longer reads or contigs
            '--algo', '1', # for small query files
            '--threads', str(util.misc.sanitize_thread_count(threads)),
            '--taxonmap', os.path.join(taxDb, 'accession2taxid', 'prot.accession2taxid.gz'),
            '--taxonnodes', os.path.join(taxDb, 'nodes.dmp'),
            '--tmpdir', tmp_dir,
            ]
        if memLimitGb:
            cmd.extend(['--block-size', str(round(memLimitGb / 5.0, 1))])
        log.debug(' '.join(cmd))
        diamond_p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        # read the output report and load into an in-memory map
        # of sequence ID -> tax ID (there shouldn't be that many sequences)
        seq_to_tax = {}
        for line in diamond_p.stdout:
            row = line.decode('UTF-8').rstrip('\n\r').split('\t')
            tax = row[1] if row[1] != '0' else '32644' # diamond returns 0 if unclassified, but the proper taxID for that is 32644
            seq_to_tax[row[0]] = tax
        if diamond_p.poll():
            raise subprocess.CalledProcessError(diamond_p.returncode, cmd)

    # copy inFasta to outFasta while prepending taxid|###| to all sequence names
    log.debug("transforming {} to {}".format(inFasta, outFasta))
    with util.file.open_or_gzopen(inFasta, 'rt') as inf:
        with util.file.open_or_gzopen(outFasta, 'wt') as outf:
            for seq in Bio.SeqIO.parse(inf, 'fasta'):
                taxid = seq_to_tax.get(seq.id, '32644') # default to "unclassified"
                for text_line in util.file.fastaMaker([(
                        '|'.join('taxid', taxid, seq.id),
                        str(seq.seq))]):
                    outf.write(text_line)

__commands__.append(('diamond_fasta', parser_diamond_fasta))


def parser_build_diamond_db(parser=argparse.ArgumentParser()):
    parser.add_argument('protein_fastas', nargs='+', help='Input protein fasta files')
    parser.add_argument('db', help='Output Diamond database file')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, build_diamond_db, split_args=True)
    return parser
def build_diamond_db(protein_fastas, db, threads=None):
    tool.diamond.Diamond().build(db, protein_fastas, options={'threads':str(util.misc.sanitize_thread_count(threads))})
__commands__.append(('build_diamond_db', parser_build_diamond_db))


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
    parser.add_argument('--outReads', help='Output LCA assignments for each read.')
    parser.add_argument('--dupeReads', help='Output LCA assignments for each read including duplicates.')
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.PicardTools.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, align_rna_metagenomics, split_args=True)
    return parser
def align_rna_metagenomics(
    inBam,
    db,
    taxDb,
    outReport,
    dupeReport=None,
    outBam=None,
    dupeReads=None,
    outReads=None,
    sensitive=None,
    JVMmemory=None,
    threads=None,
    picardOptions=None,
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

    # TODO: Use bwa.mem's min_score_to_filter argument to decrease false
    # positives in the output. Currently, it works by summing the alignment
    # score across all alignments output by bwa for each query (reads in a
    # pair, supplementary, and secondary alignments). This is not reasonable
    # for reads with secondary alignments because it will be easier for those
    # reads/queries to exceed the threshold given by the value of the argument.
    # In this context, bwa is called using '-a' as an option and its output
    # will likely include many secondary alignments. One option is to add
    # another argument to bwa.mem, similar to min_score_to_filter, that sets a
    # threshold on the alignment score of output alignments but only filters on
    # a per-alignment level (i.e., not by summing alignment scores across all
    # alignments for each query).

    aln_bam = util.file.mkstempfname('.bam')
    bwa.mem(inBam, db, aln_bam, options=bwa_opts)

    tax_db = ncbitax.TaxonomyDb(tax_dir=taxDb, load_names=True, load_nodes=True)

    if dupeReport:
        aln_bam_sorted = util.file.mkstempfname('.align_namesorted.bam')
        samtools.sort(aln_bam, aln_bam_sorted, args=['-n'], threads=threads)
        sam_lca_report(tax_db, aln_bam_sorted, outReport=dupeReport, outReads=dupeReads, unique_only=False)
        os.unlink(aln_bam_sorted)

    aln_bam_deduped = outBam if outBam else util.file.mkstempfname('.align_deduped.bam')
    opts = list(picardOptions)
    dupe_removal_out_metrics = util.file.mkstempfname('.metrics')
    pic = tools.picard.MarkDuplicatesTool()
    pic.execute([aln_bam], aln_bam_deduped, dupe_removal_out_metrics, picardOptions=opts, JVMmemory=JVMmemory)


    os.unlink(aln_bam)
    aln_bam_dd_sorted = util.file.mkstempfname('.bam')
    samtools.sort(aln_bam_deduped, aln_bam_dd_sorted, args=['-n'], threads=threads)
    sam_lca_report(tax_db, aln_bam_dd_sorted, outReport=outReport, outReads=outReads)

    if not outBam:
        os.unlink(aln_bam_deduped)
__commands__.append(('align_rna', parser_align_rna_metagenomics))


def sam_lca_report(tax_db, bam_aligned, outReport, outReads=None, unique_only=None):

    if outReads:
        lca_tsv = outReads
    else:
        lca_tsv = util.file.mkstempfname('.tsv')

    with util.file.open_or_gzopen(lca_tsv, 'wt') as lca:
        hits = sam_lca(tax_db, bam_aligned, lca, top_percent=10, unique_only=unique_only)

    with open(outReport, 'w') as f:

        for line in tax_db.kraken_dfs_report(hits):
            print(line, file=f)


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

        kraken_tool = tools.kraken.Kraken()
        kraken_tool.report(tmp_metag_combined_txt, kraken_db.name, out_kraken_summary)
__commands__.append(('report_merge', parser_metagenomic_report_merge))



def kraken_library_ids(library):
    '''Parse gi/accession from ids of fasta files in library directory. '''
    library_taxids = set()
    library_gis = set()
    library_accessions = set()
    for dirpath, dirnames, filenames in os.walk(library, followlinks=True):
        for filename in filenames:
            if not filename.endswith('.fna') and not filename.endswith('.fa') and not filename.endswith('.ffn'):
                continue
            filepath = os.path.join(dirpath, filename)
            for seqr in SeqIO.parse(filepath, 'fasta'):
                name = seqr.name
                # Search for encoded taxid
                mo = re.search('kraken:taxid\|(\d+)\|', name)
                if mo:
                    taxid = int(mo.group(1))
                    library_taxids.add(taxid)
                    continue
                # Search for gi
                mo = re.search('gi\|(\d+)\|', name)
                if mo:
                    gi = int(mo.group(1))
                    library_gis.add(gi)
                    continue
                # Search for accession
                mo = re.search('([A-Z]+\d+\.\d+)', name)
                if mo:
                    accession = mo.group(1)
                    library_accessions.add(accession)
    return library_taxids, library_gis, library_accessions


class KrakenBuildError(Exception):
    '''Error while building kraken database.'''

def parser_filter_bam_to_taxa(parser=argparse.ArgumentParser()):
    parser.add_argument('in_bam', help='Input bam file.')
    parser.add_argument('read_IDs_to_tax_IDs', help='TSV file mapping read IDs to taxIDs, Kraken-format by default. Assumes bijective mapping of read ID to tax ID.')
    parser.add_argument('out_bam', help='Output bam file, filtered to the taxa specified')
    parser.add_argument('nodes_dmp', help='nodes.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('names_dmp', help='names.dmp file from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/')
    parser.add_argument('--taxNames', nargs="+", dest="tax_names", help='The taxonomic names to include. More than one can be specified. Mapped to Tax IDs by lowercase exact match only. Ex. "Viruses" This is in addition to any taxonomic IDs provided.')
    parser.add_argument('--taxIDs', nargs="+", type=int, dest="tax_ids", help='The NCBI taxonomy IDs to include. More than one can be specified. This is in addition to any taxonomic names provided.')
    parser.add_argument('--without-children', action='store_true', dest="omit_children", help='Omit reads classified more specifically than each taxon specified (without this a taxon and its children are included).')
    parser.add_argument('--read_id_col', type=int, dest="read_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing read IDs. (default: %(default)s)', default=1)
    parser.add_argument('--tax_id_col', type=int, dest="tax_id_col", help='The (zero-indexed) number of the column in read_IDs_to_tax_IDs containing Taxonomy IDs. (default: %(default)s)', default=2)
    parser.add_argument(
        '--JVMmemory',
        default=tools.picard.FilterSamReadsTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)'
    )
    util.cmd.common_args(parser, (('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, filter_bam_to_taxa, split_args=True)
    return parser

def filter_bam_to_taxa( in_bam,
                        read_IDs_to_tax_IDs,
                        out_bam,
                        nodes_dmp,
                        names_dmp,
                        tax_names=None,
                        tax_ids=None,
                        omit_children=False,
                        read_id_col=1,
                        tax_id_col=2,
                        JVMmemory=None ):
    """
        Filter an (already classified) input bam file to only include reads that have been mapped to specified
        taxonomic IDs or scientific names. This requires a classification file, as produced
        by tools such as Kraken, as well as the NCBI taxonomy database.
    """
    tax_ids = set(tax_ids) if tax_ids else set()
    tax_names = tax_names or []
    db = ncbitax.TaxonomyDb(nodes_path=nodes_dmp, names_path=names_dmp, load_nodes=True, load_names=True)

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
            child_ids = ncbitax.collect_children(db.children, set([tax_id]))
            tax_ids_to_include |= set(child_ids)

    tax_ids_to_include = frozenset(tax_ids_to_include) # frozenset membership check slightly faster

    # perform the actual filtering to return a list of read IDs, writeen to a temp file
    with util.file.tempfname(".txt.gz") as temp_read_list:
        with open_or_gzopen(temp_read_list, "wt") as read_IDs_file:
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
                    log.debug("Found matching read ID: %s", read_id_match.group(1))
                    read_IDs_file.write(read_id_match.group(1)+"\n")
                    read_ids_written+=1

        # if we found reads matching the taxNames requested,
        if read_ids_written > 0:
            # filter the input bam to include only these
            tools.picard.FilterSamReadsTool().execute(in_bam,
                                                        False,
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
    parser.add_argument('--taxlevelFocus', dest="taxlevel_focus", help='The taxonomic heading to summarize (totals by Genus, etc.) (default: %(default)s).', default="species",
                        choices=["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"])
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

    Abundance = collections.namedtuple("Abundance", "percent,count")

    def indent_len(in_string):
        return len(in_string)-len(in_string.lstrip())

    for f in list(summary_files_in):
        sample_name, extension = os.path.splitext(f)
        sample_summary = {}
        sample_root_summary = {}
        tax_headings_copy = [s.lower() for s in tax_headings]

        with util.file.open_or_gzopen(f, 'rU') as inf:
            should_process = False
            indent_of_selection = -1
            currently_being_processed = ""
            for line in inf:
                if len(line.rstrip('\r\n').strip()) == 0:
                    continue
                csv.register_dialect('kraken_report', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
                fieldnames = ["pct_of_reads","num_reads","reads_exc_children","rank","NCBI_tax_ID","sci_name"]
                row = next(csv.DictReader([line.strip().rstrip('\n')], fieldnames=fieldnames, dialect="kraken_report"))

                indent_of_line = indent_len(row["sci_name"])
                # remove leading/trailing whitespace from each item
                row = { k:v.strip() for k, v in row.items()}

                # rows are formatted like so:
                # 0.00  16  0   D   10239     Viruses
                #
                # row["pct_of_reads"] Percentage of reads covered by the clade rooted at this taxon
                # row["num_reads"] Number of reads covered by the clade rooted at this taxon
                # row["reads_exc_children"] Number of reads assigned directly to this taxon
                # row["rank"] A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
                # row["NCBI_tax_ID"] NCBI taxonomy ID
                # row["sci_name"] indented scientific name

                # if the root-level bins (root, unclassified) should be included, do so, but bypass normal
                # stateful parsing logic since root does not have a distinct rank level
                if row["sci_name"].lower() in ["root","unclassified"] and include_root:
                    sample_root_summary[row["sci_name"]] = collections.OrderedDict()
                    sample_root_summary[row["sci_name"]][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]))
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
                        if row["rank"] == ncbitax.kraken_rank_code(taxlevel_focus):
                            same_level = True
                        if row["rank"] == "-":
                            log.warning("Non-taxonomic parent level selected")

                if should_process:
                    # skip "-" rank levels since they do not occur at the sample level
                    # otherwise include the taxon row if the rank matches the desired level of focus
                    if (row["rank"] != "-" and ncbitax.kraken_rank_code(taxlevel_focus) == row["rank"]):
                        if int(row["num_reads"])>=count_threshold:
                            sample_summary[currently_being_processed][row["sci_name"]] = Abundance(float(row["pct_of_reads"]), int(row["num_reads"]))


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


def parser_kraken_build(parser=argparse.ArgumentParser()):
    parser.add_argument('db', help='Kraken database output directory.')
    parser.add_argument('--library', help='Input library directory of fasta files. If not specified, it will be read from the "library" subdirectory of "db".')
    parser.add_argument('--taxonomy', help='Taxonomy db directory. If not specified, it will be read from the "taxonomy" subdirectory of "db".')
    parser.add_argument('--subsetTaxonomy', action='store_true', help='Subset taxonomy based on library fastas.')
    parser.add_argument('--minimizerLen', type=int, help='Minimizer length (kraken default: 15)')
    parser.add_argument('--kmerLen', type=int, help='k-mer length (kraken default: 31)')
    parser.add_argument('--maxDbSize', type=int, help='Maximum db size in GB (will shrink if too big)')
    parser.add_argument('--clean', action='store_true', help='Clean by deleting other database files after build')
    parser.add_argument('--workOnDisk', action='store_true', help='Work on disk instead of RAM. This is generally much slower unless the "db" directory lives on a RAM disk.')
    util.cmd.common_args(parser, (('threads', None), ('loglevel', None), ('version', None), ('tmp_dir', None)))
    util.cmd.attach_main(parser, kraken_build, split_args=True)
    return parser
def kraken_build(db, library, taxonomy=None, subsetTaxonomy=None,
                threads=None, workOnDisk=False,
                minimizerLen=None, kmerLen=None, maxDbSize=None, clean=False):
    '''
    Builds a kraken database from library directory of fastas and taxonomy db
    directory. The --subsetTaxonomy option allows shrinking the taxonomy to only
    include taxids associated with the library folders. For this to work, the
    library fastas must have the standard id names such as `>NC1234.1`
    accessions, `>gi|123456789|ref|XXXX||`, or custom kraken name
    `>kraken:taxid|1234|`.

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
            raise KrakenBuildError('Output db directory already contains taxonomy directory {}'.format(taxonomy_dir))
        if subsetTaxonomy:
            kraken_taxids, kraken_gis, kraken_accessions = kraken_library_ids(library)

            whitelist_taxid_f = util.file.mkstempfname()
            with open(whitelist_taxid_f, 'wt') as f:
                for taxid in kraken_taxids:
                    print(taxid, file=f)

            whitelist_gi_f = util.file.mkstempfname()
            with open(whitelist_gi_f, 'wt') as f:
                for gi in kraken_gis:
                    print(gi, file=f)

            whitelist_accession_f = util.file.mkstempfname()
            with open(whitelist_accession_f, 'wt') as f:
                for accession in kraken_accessions:
                    print(accession, file=f)

            # Context-managerize eventually
            taxonomy_tmp = tempfile.mkdtemp()
            ncbitax.subset.subset_taxonomy(
                taxonomy, taxonomy_tmp, whitelistTaxidFile=whitelist_taxid_f,
                whitelistGiFile=whitelist_gi_f, whitelistAccessionFile=whitelist_accession_f)
            shutil.move(taxonomy_tmp, taxonomy_dir)
        else:
            os.symlink(os.path.abspath(taxonomy), taxonomy_dir)
        for fn in os.listdir(os.path.join(taxonomy_dir, 'accession2taxid')):
            if not fn.endswith('accession2taxid'):
                continue

            os.symlink(os.path.join(taxonomy_dir, 'accession2taxid', fn),
                       os.path.join(taxonomy_dir, fn))
    else:
        if not taxonomy_exists:
            raise FileNotFoundError('Taxonomy directory {} not found'.format(taxonomy_dir))
        if subsetTaxonomy:
            raise KrakenBuildError('Cannot subset taxonomy if already in db folder')

    kraken_tool = tools.kraken.Kraken()
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
    kraken_tool.build(db, options=options)

    if clean:
        kraken_tool.execute('kraken-build', db, '', options={'--clean': None})
__commands__.append(('kraken_build', parser_kraken_build))



def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

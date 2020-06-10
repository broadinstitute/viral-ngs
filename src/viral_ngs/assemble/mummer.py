'''
    The MUMMER aligner
    http://mummer.sourceforge.net/
'''

import logging
import tools
import util.file
import util.misc
import os
import os.path
import random
import shutil
import subprocess
import Bio.SeqIO

TOOL_NAME = "mummer"

log = logging.getLogger(__name__)


class MummerTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(shutil.which(TOOL_NAME), require_executability=True)]
        super(MummerTool, self).__init__(install_methods=install_methods)

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        self.tool_version = subprocess.check_output([self.install_and_get_path(), '-version']).decode('UTF-8').strip()

    def executable_path(self):
        exec_path = tools.Tool.executable_path(self)
        if os.path.isdir(exec_path):
            return exec_path
        else:
            return os.path.dirname(exec_path)

    def execute(self, refFasta, qryFastas):
        toolCmd = [os.path.join(self.install_and_get_path(), 'mummer'),
            refFasta] + qryFastas
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def nucmer(self, refFasta, qryFasta, outDelta, extend=None, breaklen=None,
            maxgap=None, minmatch=None, mincluster=None):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'nucmer'),
            '--prefix={}'.format(outDelta)]
        if extend is not None:
            if not extend:
                # default behavior is --extend
                # mummer4 no longer recognizes --extend and we should only
                # specify an argument here if we want non-default behavior
                toolCmd.append('--noextend')
        if breaklen is not None:
            toolCmd.extend(['--breaklen', str(breaklen)])
        if maxgap is not None:
            toolCmd.extend(['--maxgap', str(maxgap)])
        if minmatch is not None:
            toolCmd.extend(['--minmatch', str(minmatch)])
        if mincluster is not None:
            toolCmd.extend(['--mincluster', str(mincluster)])
        toolCmd.extend([refFasta, qryFasta])
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def promer(self, refFasta, qryFasta, outDelta, extend=None, breaklen=None,
            maxgap=None, minmatch=None, mincluster=None):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'promer'),
            '--prefix={}'.format(outDelta)]
        if extend is not None:
            if not extend:
                # default behavior is --extend
                # mummer4 no longer recognizes --extend and we should only
                # specify an argument here if we want non-default behavior
                toolCmd.append('--noextend')
        if breaklen is not None:
            toolCmd.extend(['--breaklen', str(breaklen)])
        if maxgap is not None:
            toolCmd.extend(['--maxgap', str(maxgap)])
        if minmatch is not None:
            toolCmd.extend(['--minmatch', str(minmatch)])
        if mincluster is not None:
            toolCmd.extend(['--mincluster', str(mincluster)])
        toolCmd.extend([refFasta, qryFasta])
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def delta_filter(self, inDelta, outDelta):
        toolCmd = [os.path.join(self.install_and_get_path(), 'delta-filter'),
            '-q', inDelta]
        log.debug(' '.join(toolCmd))
        with open(outDelta, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)

    def show_tiling(self, inDelta, outTiling, outFasta=None,
            circular=False, min_pct_id=None, min_contig_len=None,
            min_pct_contig_aligned=None, tab_delim=False,
            min_contig_coverage_diff=None):
        opts = []
        if circular:
            opts.append('-c')
        if tab_delim:
            opts.append('-a')
        if min_pct_id is not None:
            opts.append('-i')
            opts.append(str(100.0 * min_pct_id))
        if min_contig_len is not None:
            opts.append('-l')
            opts.append(str(min_contig_len))
        if min_pct_contig_aligned is not None:
            opts.append('-v')
            opts.append(str(100.0 * min_pct_contig_aligned))
        if min_contig_coverage_diff is not None:
            opts.append('-V')
            opts.append(str(100.0 * min_contig_coverage_diff))
        toolCmd = [os.path.join(self.install_and_get_path(), 'show-tiling')]
        if outFasta:
            toolCmd = toolCmd + ['-p', outFasta]
        toolCmd = toolCmd + opts + [inDelta]
        log.debug(' '.join(toolCmd))
        with open(outTiling, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)

    def trim_contigs(self, refFasta, contigsFasta, outFasta,
            aligner='nucmer', circular=False, extend=False, breaklen=None,
            min_pct_id=0.6, min_pct_contig_aligned=0.5, min_contig_len=200):
        ''' Align contigs with MUMmer and trim off the unused portions.
        '''
        # run MUMmer to get best alignments
        if aligner=='nucmer':
            aligner = self.nucmer
        elif aligner=='promer':
            aligner = self.promer
        else:
            raise NameError()
        delta_1 = util.file.mkstempfname('.delta')
        delta_2 = util.file.mkstempfname('.delta')
        tiling = util.file.mkstempfname('.tiling')
        aligner(refFasta, contigsFasta, delta_1, extend=extend, breaklen=breaklen)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, circular=circular, tab_delim=True,
            min_pct_id=min_pct_id, min_contig_len=min_contig_len,
            min_pct_contig_aligned=min_pct_contig_aligned)
        os.unlink(delta_1)
        os.unlink(delta_2)

        # get trim boundaries per contig
        seq_bounds = {}
        with open(tiling, 'rt') as inf:
            for line in inf:
                row = line.rstrip().split('\t')
                seq_bounds.setdefault(row[-1], [float('inf'), 0])
                start = min(int(row[2]), int(row[3]))
                stop  = max(int(row[2]), int(row[3]))
                seq_bounds[row[-1]][0] = min(seq_bounds[row[-1]][0], start)
                seq_bounds[row[-1]][1] = max(seq_bounds[row[-1]][1], stop)
        os.unlink(tiling)

        # trim contigs
        out = []
        with open(contigsFasta, 'rt') as inf:
            for record in Bio.SeqIO.parse(inf, 'fasta'):
                if record.id in seq_bounds:
                    seq = record.seq[seq_bounds[record.id][0]-1:seq_bounds[record.id][1]]
                    record.seq = seq
                    out.append(record)
        with open(outFasta, 'wt') as outf:
            Bio.SeqIO.write(out, outf, 'fasta')

    def scaffold_contigs(self, refFasta, contigsFasta, outFasta,
            aligner='nucmer', circular=False, extend=None, breaklen=None,
            maxgap=None, minmatch=None, mincluster=None,
            min_contig_coverage_diff=0.0,
            min_pct_id=0.6, min_pct_contig_aligned=None, min_contig_len=200):
        ''' Use MUMmer's pseudomolecule feature to scaffold contigs
            onto a reference genome.
        '''
        if aligner=='nucmer':
            aligner = self.nucmer
        elif aligner=='promer':
            aligner = self.promer
        else:
            raise NameError()
        delta_1 = util.file.mkstempfname('.delta')
        delta_2 = util.file.mkstempfname('.delta')
        tiling = util.file.mkstempfname('.tiling')
        aligner(refFasta, contigsFasta, delta_1, extend=extend, breaklen=breaklen,
            maxgap=maxgap, minmatch=minmatch, mincluster=mincluster)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, outFasta=outFasta, circular=circular,
            min_pct_id=min_pct_id, min_contig_len=min_contig_len,
            min_contig_coverage_diff=min_contig_coverage_diff,
            min_pct_contig_aligned=min_pct_contig_aligned)
        os.unlink(delta_1)
        os.unlink(delta_2)
        os.unlink(tiling)

    def scaffold_contigs_custom(self, refFasta, contigsFasta, outFasta,
            outAlternateContigs=None,
            aligner='nucmer', extend=None, breaklen=None,
            maxgap=None, minmatch=None, mincluster=None,
            min_contig_coverage_diff=0.0,
            min_pct_id=0.6, min_pct_contig_aligned=None, min_contig_len=200,
            ambig_max_aligns=2, ambig_max_lens=1, ambig_max_frac=.01):
        ''' Re-implement a less buggy version of MUMmer's pseudomolecule
            feature to scaffold contigs onto a reference genome.
        '''

        # create tiling path with nucmer/promer and show-tiling
        if aligner=='nucmer':
            aligner = self.nucmer
        elif aligner=='promer':
            aligner = self.promer
            raise NotImplementedError('we have not implemented a show-aligns file reader that works for protein alignments')
        else:
            raise NameError()
        delta_1 = util.file.mkstempfname('.delta')
        delta_2 = util.file.mkstempfname('.delta')
        tiling = util.file.mkstempfname('.tiling')
        aligner(refFasta, contigsFasta, delta_1, extend=extend, breaklen=breaklen,
            maxgap=maxgap, minmatch=minmatch, mincluster=mincluster)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, tab_delim=True,
            min_pct_id=min_pct_id,
            min_contig_len=min_contig_len,
            min_contig_coverage_diff=min_contig_coverage_diff,
            min_pct_contig_aligned=min_pct_contig_aligned)
        os.unlink(delta_1)

        # load intervals into a FeatureSorter
        fs = util.misc.FeatureSorter()
        with util.file.open_or_gzopen(tiling, 'rU') as inf:
            for line in inf:
                row = line.rstrip('\n\r').split('\t')
                c = row[11]
                start, stop = (int(row[0]), int(row[1]))
                alt_seq = (row[12], int(row[2]), int(row[3]))
                if stop<start:
                    raise ValueError()
                if alt_seq[2]<alt_seq[1]:
                    s = '-'
                else:
                    s = '+'
                fs.add(c, start, stop, strand=s, other=alt_seq)
                log.info("mummer alignment %s:%d-%d - %s:%d-%d (%s)" % (
                    c, start, stop,
                    alt_seq[0], alt_seq[1], alt_seq[2],
                    s
                ))
        os.unlink(tiling)

        # load all contig-to-ref alignments into AlignsReaders
        alnReaders = {}
        aln_file = util.file.mkstempfname('.aligns')
        for c, start, stop, strand, other in fs.get_features():
            chr_pair = (c, other[0])
            if chr_pair not in alnReaders:
                toolCmd = [os.path.join(self.install_and_get_path(), 'show-aligns'),
                    '-r', delta_2, chr_pair[0], chr_pair[1]]
                #log.debug(' '.join(toolCmd))
                with open(aln_file, 'wt') as outf:
                    subprocess.check_call(toolCmd, stdout=outf)
                alnReaders[chr_pair] = AlignsReader(aln_file, refFasta)
        os.unlink(aln_file)
        os.unlink(delta_2)

        # for each chromosome, create the scaffolded sequence and write everything to fasta
        alternate_contigs = []
        with open(outFasta, 'wt') as outf:
            for c in [seqObj.id for seqObj in Bio.SeqIO.parse(refFasta, 'fasta')]:
                if c not in fs.get_seqids():
                    # segment c could not be assembled; skip it with the plan
                    # for assembly to fail later since (# reference segments) !=
                    # (# assembled segments)
                    continue

                def n_diff_vals(*vals): return len(set(vals))
                def n_diff_lens(seqs): return n_diff_vals(*map(len, seqs))
                def frac_unambig(seqs):
                    """Given a list of seqs of the same length, return the fraction of positions on which they all agree"""
                    util.misc.chk(n_diff_lens(alt_seqs_f) == 1, 'ambig_max_lens>1 not currently supported')
                    n_tot = len(seqs[0])
                    n_unambig = list(map(n_diff_vals, *seqs)).count(1)
                    return float(n_unambig) / float(n_tot or 1.0)

                # construct scaffolded sequence for this chromosome
                seq = []
                for _, left, right, n_features, features in fs.get_intervals(c):
                    # get all proposed sequences for this specific region
                    alt_seqs = []
                    for consider_ambig_aligns in (False, True):
                        for f in features:
                            alt_seqs_f = alnReaders[(c, f[-1][0])].retrieve_alts_by_ref(left, right, aln_start=f[1], aln_stop=f[2])
                            if len(alt_seqs_f) == 1:
                                alt_seqs.append(alt_seqs_f[0])
                            elif consider_ambig_aligns:
                                if len(alt_seqs_f) <= ambig_max_aligns and n_diff_lens(alt_seqs_f) <= ambig_max_lens and \
                                   frac_unambig(alt_seqs_f) > (1.0 - ambig_max_frac):
                                    alt_seqs.append(alt_seqs_f[0])
                                    log.info("using ambiguous alignment to ref seq {} at [{},{}]".format(c, f[1], f[2]))
                                else:
                                    log.warning("dropping ambiguous alignment to ref seq {} at [{},{}]".format(c, f[1], f[2]))
                        if alt_seqs:
                            # if have a non-unambiguous alignment, don't consider ambiguous ones
                            break

                    # pick the "right" one and glue together into a chromosome
                    ranked_unique_seqs = contig_chooser(alt_seqs, right-left+1, "%s:%d-%d" % (c, left, right))
                    seq.append(ranked_unique_seqs[0])
                    # emit the "alternates" in a separate file
                    if outAlternateContigs and len(ranked_unique_seqs) > 1:
                        alternate_contigs.append((c, left, right, ranked_unique_seqs))

                # write this chromosome to fasta file
                for line in util.file.fastaMaker([(str(c)+"_contigs_ordered_and_oriented", ''.join(seq))]):
                    outf.write(line)

        # if alternate scaffolds exist, emit to output fasta file (if specified)
        if outAlternateContigs and alternate_contigs:
            log.info("emitting alternative scaffold sequences to {}".format(outAlternateContigs))
            with open(outAlternateContigs, 'wt') as outf:
                for c, left, right, seqs in alternate_contigs:
                    for line in util.file.fastaMaker([
                        ("{}:{}-{}_option_{}".format(c,left,right,i), s)
                        for i,s in enumerate(seqs)]):
                        outf.write(line)


    def align_one_to_one(self, refFasta, otherFasta, outFasta):
        ''' Take two sequences, align with nucmer, and produce
            an aligned-fasta file of the two based on show-aligns
            output.
        '''
        # grab seq_ids (very inefficient, but whatever)
        ref_id = Bio.SeqIO.read(refFasta, 'fasta').id
        query_id = Bio.SeqIO.read(otherFasta, 'fasta').id

        # run nucmer
        delta_1 = util.file.mkstempfname('.delta')
        delta_2 = util.file.mkstempfname('.delta')
        self.nucmer(refFasta, otherFasta, delta_1)
        self.delta_filter(delta_1, delta_2)

        # run show-aligns
        aln_file = util.file.mkstempfname('.aligns')
        toolCmd = [os.path.join(self.install_and_get_path(), 'show-aligns'),
            '-r', delta_2, ref_id, query_id]
        log.debug(' '.join(toolCmd))
        with open(aln_file, 'wt') as outf:
            subprocess.check_call(toolCmd, stdout=outf)
        alns = AlignsReader(aln_file, refFasta)

        # write fasta
        seqs = [[alns.seq_ids[0], []], [alns.seq_ids[1], []]]
        for a in alns.get_intervals():
            seqs[0][1].append(a[6])
            seqs[1][1].append(a[7])
        seqs[0][1] = ''.join(seqs[0][1])
        seqs[1][1] = ''.join(seqs[1][1])
        util.file.makeFastaFile(seqs, outFasta)

        # cleanup
        for fn in (delta_1, delta_2, aln_file):
            os.unlink(fn)

def contig_chooser(alt_seqs, ref_len, coords_debug=""):
    ''' Our little heuristic to choose an alternative sequence from a pile
        of alignments to a reference. Takes a list of strings (one string per
        contig). This method will choose a single sequence from the input:
            1. if there are no alt_seqs, emit a stretch of Ns, same length as ref
            2. if there is only one alt_seq, emit that one
            3. if there are many alt_seqs, emit the most popular (if there is one)
            4. otherwise, if there is a most popular sequence length (including
                the ref_len as one vote), filter the alt_seqs to those of the
                most popular length and emit the most popular sequence (if there
                is one), otherwise, choose randomly amongst the remainder
            5. otherwise, choose randomly amongst the same-as-ref length sequences
            6. or just choose randomly if there are no same-as-ref length sequences
        The output will be a list of unique strings, where the first string
        is the "chosen" sequence, and the remaining strings are the alternative
        sequences (in no particular order, but guaranteed to be unique).
    '''
    if not alt_seqs:
        # no contigs align here, emit Ns of appropriate length
        new_seq = 'N' * ref_len
        other_seqs = []
    elif len(alt_seqs) == 1:
        # only one contig aligns here
        new_seq = alt_seqs[0]
        other_seqs = []
    else:
        # multiple contigs align here
        ranks = list(sorted(util.misc.histogram(alt_seqs).items(),
            key=lambda x:x[1], reverse=True))
        other_seqs = list(s for s,n in ranks)
        if len(ranks)==1:
            # only one unique sequence exists
            new_seq = ranks[0][0]
        elif ranks[0][1]>ranks[1][1]:
            # clear winner: a most popular sequence exists
            new_seq = ranks[0][0]
        else:
            # multiple possible replacement sequences
            len_ranks = list(sorted(util.misc.histogram(
                [len(s) for s in alt_seqs] + [ref_len] # let the reference have one vote
                ).items(), key=lambda x:x[1], reverse=True))
            if len(len_ranks)==1 or len_ranks[0][1]>len_ranks[1][1]:
                # a most popular replacement length exists
                # remove all alt_seqs that don't use that length
                alt_seqs = list(s for s in alt_seqs if len(s)==len_ranks[0][0])
                assert alt_seqs
                ranks = list(sorted(util.misc.histogram(alt_seqs).items(),
                    key=lambda x:x[1], reverse=True))
                if len(ranks)==1 or ranks[0][1]>ranks[1][1]:
                    # clear winner amongst remaining sequences of most popular length
                    new_seq = ranks[0][0]
                else:
                    # more complicated scenario. choose randomly.
                    # perhaps in future, vote based on aligned read count?
                    if len(alt_seqs)>1:
                        log.warning("choosing random contig from %d choices of most popular length in %s" % (len(alt_seqs), coords_debug))
                    new_seq = random.choice(alt_seqs)
            else:
                # no clear winner on replacement length
                alt_ref_len_seqs = list(s for s in alt_seqs if len(s)==ref_len)
                if alt_ref_len_seqs:
                    # choose randomly among same-as-ref-length sequences
                    alt_seqs = alt_ref_len_seqs
                    if len(alt_seqs)>1:
                        log.warning("choosing random contig from %d choices of reference length in %s" % (len(alt_seqs), coords_debug))
                    new_seq = random.choice(alt_seqs)
                else:
                    # no clear winner and all replacement lengths are different from reference length
                    # just choose randomly
                    if len(alt_seqs)>1:
                        log.warning("choosing random contig from %d choices in %s" % (len(alt_seqs), coords_debug))
                    new_seq = random.choice(alt_seqs)
        other_seqs = list(s for s in other_seqs if s!=new_seq)
    return [new_seq] + other_seqs


class AmbiguousAlignmentException(Exception):
    pass

class AlignsReader(object):
    ''' This class assists in the parsing and reading of show-aligns output.
    '''
    def __init__(self, aligns_file, ref_fasta=None):
        self.aligns_file = aligns_file
        self.ref_fasta = ref_fasta
        self.reference_seq = None
        self.alignments = []
        self.seq_ids = []
        self._load_align()
        self._load_fastas()

    def _load_align(self):
        with open(self.aligns_file, 'rt') as inf:
            # read ref_fasta, query_fasta from header
            header = inf.readline().strip().split()
            assert len(header) == 2
            if self.ref_fasta is None:
                self.ref_fasta = header[0]

            # iterate row by row
            mode = 'start'
            coords = None
            seqs = None
            for line in inf:
                line = line.rstrip()
                if not line:
                    pass # empty lines
                elif not line.strip('='):
                    pass # header and footer of file
                elif line[0] in (' ', '\t'):
                    pass # describes mismatches in alignment
                elif line.startswith('--'):
                    if line.startswith('-- Alignments between '):
                        assert mode == 'start'
                        mode = 'between'
                        self.seq_ids = line[22:].split(' and ')
                        assert len(self.seq_ids) == 2
                    elif line.startswith('-- BEGIN alignment [ '):
                        assert mode == 'between'
                        mode = 'align'
                        coords = list(x.split() for x in line[21:-2].split(' | '))
                        assert len(coords) == 2 and coords[0][2] == coords[1][2] == '-'
                        assert coords[0][0] == '+1', "error with line: %s" % line
                        seqs = [[], []]
                        align_lines = 0
                    elif line.startswith('--   END alignment [ '):
                        assert mode == 'align'
                        mode = 'between'
                        new_coords = list(x.split() for x in line[21:-2].split(' | '))
                        assert coords == new_coords, "error: %s != %s" % (new_coords, coords)
                        assert len(seqs[0]) == len(seqs[1])
                        seqs = list(''.join(x) for x in seqs)
                        assert len(seqs[0]) == len(seqs[1])
                        self.alignments.append([coords[0][0], int(coords[0][1]), int(coords[0][3]),
                               coords[1][0], int(coords[1][1]), int(coords[1][3]),
                               seqs[0], seqs[1]])
                        coords = None
                        seqs = None
                    else:
                        raise AssertionError("file format: line '%s'" % line)
                else:
                    # read in one line of an alignment
                    assert mode == 'align', "file format: line '%s' before alignment begins" % line
                    seq_str = line.split()[1].upper().replace('.','-')
                    seqs[align_lines % 2].append(seq_str)
                    align_lines += 1

    def _load_fastas(self):
        assert self.ref_fasta and self.seq_ids
        self.reference_seq = Bio.SeqIO.index(self.ref_fasta, 'fasta')[self.seq_ids[0]]

    def get_alignments(self):
        for a in self.alignments:
            yield a

    def get_intervals(self):
        prev = None
        for a in self.alignments:
            cur = a[1:3]
            assert cur[1] >= cur[0]
            if prev is None:
                if cur[0] > 1:
                    # emit leading reference sequence before first alignment
                    yield self._dummy_row(1, cur[0]-1, '-')
            else:
                assert cur[0] > prev[1], "overlaps not allowed"
                if cur[0] > prev[1] + 1:
                    # emit gap between alignments
                    yield self._dummy_row(prev[1]+1, cur[0]-1, 'N')
            # emit actual alignment
            yield a
            prev = cur
        if prev and prev[1] < len(self.reference_seq):
            # emit trailing reference sequence after last alignment
            yield self._dummy_row(prev[1]+1, len(self.reference_seq), '-')

    def _dummy_row(self, start, stop, filler='N'):
        return ['+1', start, stop, '+1', start, stop,
            self.get_ref_seq(start, stop),
            filler * (stop-start+1)]

    def get_ref_seq(self, start, stop):
        ''' Retrieve a sub-sequence from the reference (1st) sequence in the
            alignment using coordinates relative to the reference sequence.
            No gaps will be emitted.
        '''
        return str(self.reference_seq.seq[start-1:stop])

    def retrieve_alts_by_ref(self, start, stop, aln_start=None, aln_stop=None):
        ''' Retrieve sub-sequence(s) from the alternate (2nd) sequence in the
            alignment using coordinates relative to the reference sequence.
            No gaps will be emitted.
            Required: start-stop interval must be wholly contained within
            an alignment.
        '''

        # grab the one alignment that contains this window
        alns = list(a for a in self.alignments if a[1]<=start and a[2]>=stop)
        if aln_start is not None and aln_stop is not None:
            # if specified, restrict to a specific alignment that comes from show-tiling
            # (sometimes show-aligns is more promiscuous than show-tiling)
            new_alns = []
            for a in alns:
                if a[1] > aln_start or a[2] < aln_stop:
                    log.debug("dropping undesired alignment: %s(%s):%s-%s to %s(%s):%s-%s (%s:%s-%s requested)",
                        self.seq_ids[0], a[0], a[1], a[2],
                        self.seq_ids[1], a[3], a[4], a[5],
                        self.seq_ids[0], aln_start, aln_stop)
                else:
                    new_alns.append(a)
            alns = new_alns
        if len(alns) != 1:
            log.warning("invalid %s:%d-%d -> %s specified, %d alignments found that contain it",
                self.seq_ids[0], start, stop, self.seq_ids[1], len(alns))
            for aln in alns:
                log.debug("alignment: %s", str(aln[:6]))

        return [self._aln_to_alt_seq(aln, start, stop) for aln in alns]

    def _aln_to_alt_seq(self, aln, start, stop):
        """Given an alignment of a contig to ref, return the contig sequence aligned to a given stretch of ref"""
        ref_l, ref_r, ref_seq, alt_seq = (aln[1], aln[2], aln[-2], aln[-1])

        # convert desired start/stop relative to this reference window
        #  such that 0 <= start <= stop <= ref_r-ref_l+1
        aln_start = start - ref_l
        aln_stop = stop - ref_l

        # travel down alignment until we've reached the left edge
        #  (because of gaps, you must check each position one by one)
        #  end loop when ref_seq[:i_left] contains {aln_start} bases
        n_ref_bases = 0
        i_left = 0
        while n_ref_bases < aln_start:
            if ref_seq[i_left] != '-':
                n_ref_bases += 1
            i_left += 1

        # travel down alignment until we've reached the right edge
        #  (because of gaps, you must check each position one by one)
        #  end loop when ref_seq[:i_right] contains {aln_stop} bases
        i_right = i_left
        while n_ref_bases < aln_stop:
            if ref_seq[i_right] != '-':
                n_ref_bases += 1
            i_right += 1
        # consume and include any trailing gaps
        while i_right < len(ref_seq) and ref_seq[i_right] == '-':
            i_right += 1

        # grab the alternate sequence and strip gaps
        return alt_seq[i_left:i_right+1].replace('-','')

'''
    The MUMMER aligner
    http://mummer.sourceforge.net/
'''

import logging
import tools
import util.file
import os
import os.path
import subprocess
import Bio.SeqIO

tool_version = '3.23'
url = 'http://iweb.dl.sourceforge.net/project/mummer/mummer/{ver}/MUMmer{ver}.tar.gz'

log = logging.getLogger(__name__)


class MummerTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [
                tools.DownloadPackage(url.format(ver=tool_version),
                                      'MUMmer{}'.format(tool_version),
                                      post_download_command='cd MUMmer{}; make -s'.format(tool_version),
                                      verifycmd='{}/MUMmer{}/mummer -h > /dev/null 2>&1'.format(
                                          util.file.get_build_path(), tool_version))
                ]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return tool_version

    def execute(self, refFasta, qryFastas):
        toolCmd = [os.path.join(self.install_and_get_path(), 'mummer'),
            refFasta] + qryFastas
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def nucmer(self, refFasta, qryFasta, outDelta, extend=None, breaklen=None):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'nucmer'),
            '--prefix={}'.format(outDelta)]
        if extend is not None:
            if extend:
                toolCmd.append('--extend')
            else:
                toolCmd.append('--noextend')
        if breaklen is not None:
            toolCmd.extend(['--breaklen', str(breaklen)])
        toolCmd.extend([refFasta, qryFasta])
        log.debug(' '.join(toolCmd))
        subprocess.check_call(toolCmd)

    def promer(self, refFasta, qryFasta, outDelta, extend=None, breaklen=None):
        if not outDelta.endswith('.delta'):
            raise Exception()
        outDelta = outDelta[:-6]
        toolCmd = [os.path.join(self.install_and_get_path(), 'promer'),
            '--prefix={}'.format(outDelta)]
        if extend is not None:
            if extend:
                toolCmd.append('--extend')
            else:
                toolCmd.append('--noextend')
        if breaklen is not None:
            toolCmd.extend(['--breaklen', str(breaklen)])
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
            tab_delim=False):
        opts = []
        if circular:
            opts.append('-c')
        if tab_delim:
            opts.append('-a')
        if min_pct_id is not None:
            opts.append('-i')
            opts.append(str(min_pct_id))
        if min_contig_len is not None:
            opts.append('-l')
            opts.append(str(min_contig_len))
        toolCmd = [os.path.join(self.install_and_get_path(), 'show-tiling')]
        if outFasta:
            toolCmd = toolCmd + ['-p', outFasta]
        toolCmd = toolCmd + opts + [inDelta]
        log.debug(' '.join(toolCmd))
        with open(outTiling, 'w') as outf:
            subprocess.check_call(toolCmd, stdout=outf)
    
    def trim_contigs(self, refFasta, contigsFasta, outFasta,
            aligner='nucmer', circular=False, extend=False, breaklen=None,
            min_pct_id=0.6, min_contig_len=200):
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
        aligner(refFasta, contigsFasta, delta_1, extend=extend)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, circular=circular, tab_delim=True,
            min_pct_id=min_pct_id, min_contig_len=min_contig_len)
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
            min_pct_id=0.6, min_contig_len=200):
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
        aligner(refFasta, contigsFasta, delta_1, extend=extend)
        self.delta_filter(delta_1, delta_2)
        self.show_tiling(delta_2, tiling, outFasta=outFasta, circular=circular,
            min_pct_id=min_pct_id, min_contig_len=min_contig_len)
        os.unlink(delta_1)
        os.unlink(delta_2)
        os.unlink(tiling)
        
    def align_one_to_one(self, refFasta, otherFasta, outFasta):
        ''' Take a pseudomolecule (with N's and ragged ends) and produce
            an aligned-fasta file against a reference genome
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


class AlignsReader(object):
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
                        assert coords[0][0] == '+1'
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
                        self.alignments.append([
                            coords[0][0], int(coords[0][1]), int(coords[0][3]),
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
        assert self.ref_fasta
        self.reference_seq = Bio.SeqIO.read(self.ref_fasta, 'fasta')
        
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
        return str(self.reference_seq.seq[start-1:stop])

        
        


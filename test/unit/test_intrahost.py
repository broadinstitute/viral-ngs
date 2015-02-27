# Unit tests for intrahost.py

__author__ = "PLACEHOLDER"

import intrahost, util.file, test
import os, shutil, tempfile, itertools, argparse, unittest
import Bio, Bio.SeqRecord, Bio.Seq

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in intrahost.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

def makeTempFasta(seqs):
    fn = util.file.mkstempfname('.fasta')
    with open(fn, 'wt') as outf:
        for line in util.file.fastaMaker(seqs):
            outf.write(line)
    return fn

class MockVphaserOutput:
    ''' This creates test data that pretends to be the output from
        tools.vphaser2.Vphaser2Tool.iterate
    '''
    def __init__(self):
        self.isnvs = {}
        self.chroms = []
    def add_snp(self, chrom, pos, acounts):
        ''' Add an iSNP at this chrom,pos. acounts is a list of triples:
            (allele, fwd count, rev count).
        '''
        assert type(pos) == int and pos>0 and len(acounts)>1
        for a,f,r in acounts:
            assert a in ('A','C','G','T')
            assert f>=0 and r>=0 and f+r>0
        acounts = reversed(sorted((f+r,a,f,r) for a,f,r in acounts))
        acounts = list((a,f,r) for n,a,f,r in acounts)
        if chrom not in self.chroms:
            self.chroms.append(chrom)
        self.isnvs.setdefault(chrom, {})
        self.isnvs[chrom].setdefault(pos, {})
        self.isnvs[chrom][pos]['snp'] = acounts
    def add_indel(self, chrom, pos, acounts):
        ''' Add an iSNP at this chrom,pos. acounts is a list of triples:
            (allele, fwd count, rev count). allele is simply a small
            sequence with no "I" or "D" prefixing it (we'll add the I/D).
        '''
        assert type(pos) == int and pos>0 and len(acounts)>1
        for a,f,r in acounts:
            assert type(a) == str
            assert f>=0 and r>=0 and f+r>0
        assert '' in set(a for a,f,r in acounts)
        acounts = reversed(sorted((f+r,a,f,r) for a,f,r in acounts))
        acounts = list([a,f,r] for n,a,f,r in acounts)
        
        # vphaser funniness here
        if acounts[0][0] == '':
            # this is a set of insertions against the consensus
            acounts[0][0] = 'd'
            for i in range(1, len(acounts)):
                acounts[i][0] = 'I' + acounts[i][0]
        else:
            # this is a deletion against the consensus
            for i in range(1, len(acounts)):
                n_deleted = len(acounts[0][0]) - len(acounts[i][0])
                assert n_deleted > 0
                assert acounts[0][0][n_deleted:] == acounts[i][0]
                acounts[i][0] = 'D' + str(n_deleted)
            acounts[0][0] = 'i'
        
        if chrom not in self.chroms:
            self.chroms.append(chrom)
        self.isnvs.setdefault(chrom, {})
        self.isnvs[chrom].setdefault(pos, {})
        self.isnvs[chrom][pos]['lp'] = acounts
    def __iter__(self):
        for c in self.chroms:
            for p in sorted(self.isnvs[c].keys()):
                for model in self.isnvs[c][p].keys():
                    acounts = self.isnvs[c][p][model]
                    mac = sum(f+r for a,f,r in acounts[1:])
                    tot = sum(f+r for a,f,r in acounts)
                    yield [c, str(p), acounts[1][0], acounts[0][0],
                        '0.5', model, str(float(mac)/tot*100.0)] \
                        + ['{}:{}:{}'.format(a,f,r) for a,f,r in acounts]
    def dump_tmp_file(self):
        fn = util.file.mkstempfname('.txt')
        with open(fn, 'wt') as outf:
            for row in self:
                outf.write('\t'.join(map(str, row[:7] + ['', ''] + row[7:])) + '\n')
        return fn


class TestPerSample(test.TestCaseWithTmp):
    ''' This tests step 1 of the iSNV calling process
        (intrahost.vphaser_one_sample), which runs V-Phaser2 on
        a single sample, reformats the output slightly, and performs
        strand-bias filtering and adds library-bias statistics.
        These unit tests mock the vphaser tool output and just test
        the filtering/statistics/etc.
    '''
    def test_single_strand_bias_hard_filter(self):
        data = MockVphaserOutput()
        data.add_snp('c1', 100, [('A',10,20), ('T',5,2), ('C',30,500), ('G',60,40)])
        data.add_snp('c2', 100, [('C',10,2), ('T',2,8)])
        output = list(intrahost.filter_strand_bias(data))
        expected = ['c1', '100', 'A', 'G', None, 'snp', 23.076923076923078, 'G:60:40', 'A:10:20']
        self.assertEqual(len(output), 1)
        self.assertEqual(output[0][:4], expected[:4])
        self.assertEqual(output[0][5], expected[5])
        self.assertAlmostEqual(float(output[0][6]), expected[6], places=4)
        self.assertEqual(output[0][7:], expected[7:])
    
    @unittest.skip('not implemented')
    def test_single_lib(self):
        data = MockVphaserOutput()
        
    
    @unittest.skip('not implemented')
    def test_multi_lib(self):
        lib1 = MockVphaserOutput()
        lib2 = MockVphaserOutput()
        lib3 = MockVphaserOutput()
        pass



@unittest.skip('not implemented')
class TestVcfMerge(test.TestCaseWithTmp):
    ''' This tests step 2 of the iSNV calling process
        (intrahost.merge_to_vcf), which gets really nasty and tricky
        and has lots of edge cases. These unit tests mock the vphaser
        tool output and just test the merge and VCF stuff.
    '''
    def test_simple_snps(self):
        pass
    def test_sample_major_allele_not_ref_allele(self):
        pass
    def test_simple_insertions(self):
        # IA, ITCG, etc
        pass
    def test_simple_deletions(self):
        # D1, D2, etc...
        pass
    def test_deletion_spans_deletion(self):
        # sample assembly has deletion against reference and isnv deletes even more
        # POS is anchored right before the deletion
        # REF:  ATCGTTCA
        # S1:   ATCG--CA
        # isnv:       x  (position 5, D1)
        ref = makeTempFasta([('ref1', 'ATCGTTCA')])
        s1  = makeTempFasta([('s1_1', 'ATCGCA')])
        isnvs = MockVphaserOutput()
        isnvs.add_indel('s1_1', 5, [('C', 90, 90), ('', 10, 10)])
        outVcf = util.file.mkstempfname('.vcf')
        intrahost.merge_to_vcf(ref, outVcf, ['s1'], isnvs.dump_tmp_file(), [s1])
        # TO DO: test expected output
        pass
    def test_insertion_spans_deletion(self):
        # sample assembly has deletion against reference and isnv inserts back into it
        # POS is anchored right before the deletion
        # REF:  ATCGTTGA
        # S1:   ATCG--GA
        # isnv:     T     (position 4, IT)
        # isnv:     TT    (position 4, ITT)
        # isnv:     TTC   (position 4, ITTC)
        pass
    def test_deletion_within_insertion(self):
        # sample assembly has insertion against reference and isnv deletes from it
        # REF:  ATCG--GA
        # S1:   ATCGTTGA
        # isnv:      xx   (position 6, D2)
        # isnv:       x   (position 7, D1)
        # isnv:      x    (position 6, D1)
        # isnv:     x     (position 5, D1)
        # isnv:    x      (position 4, D1)
        # isnv:    xx     (position 4, D2)
        # isnv:    xxx    (position 4, D3)
        # isnv:    xxxx   (position 4, D4)
        pass
    def test_insertion_within_insertion(self):
        # sample assembly has insertion against reference and isnv puts even more in
        # REF:  ATCG--GA
        # S1:   ATCGTTGA
        # isnv:           (position 4, IA)
        # isnv:           (position 5, IA)
        # isnv:           (position 6, IA)
        pass
    def test_indel_collapse(self):
        # vphaser describes insertions and deletions separately
        # test appropriate collapse of coincident insertions and deletions into
        # a single output VCF row
        # isnv:           (position 4, IA, IAT)
        # isnv:           (position 5, D1, D2)
        pass




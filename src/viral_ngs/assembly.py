#!/usr/bin/env python
''' This script contains a number of utilities for viral sequence assembly
    from NGS reads.  Primarily used for Lassa and Ebola virus analysis in
    the Sabeti Lab / Broad Institute Viral Genomics.
'''

__author__ = "dpark@broadinstitute.org, rsealfon@broadinstitute.org"
__commands__ = []

import argparse, logging, random, os, os.path, shutil, subprocess, glob
import Bio.AlignIO, Bio.SeqIO, Bio.Data.IUPACData
import util.cmd, util.file, util.vcf
import read_utils, taxon_filter
import tools, tools.picard, tools.samtools, tools.gatk, tools.novoalign
import tools.trinity, tools.mosaik, tools.muscle

log = logging.getLogger(__name__)


def trim_rmdup_subsamp_reads(inBam, clipDb, outBam, n_reads=100000):
    ''' Take reads through Trimmomatic, Prinseq, and subsampling.
        This should probably move over to read_utils or taxon_filter.
    '''
    
    # BAM -> fastq
    infq = list(map(util.file.mkstempfname, ['.in.1.fastq', '.in.2.fastq']))
    tools.picard.SamToFastqTool().execute(inBam, infq[0], infq[1])
    
    # Trimmomatic
    trimfq = list(map(util.file.mkstempfname, ['.trim.1.fastq', '.trim.2.fastq']))
    taxon_filter.trimmomatic(infq[0], infq[1], trimfq[0], trimfq[1], clipDb)
    os.unlink(infq[0])
    os.unlink(infq[1])
    
    # Prinseq
    rmdupfq = list(map(util.file.mkstempfname, ['.rmdup.1.fastq', '.rmdup.2.fastq']))
    read_utils.rmdup_prinseq_fastq(trimfq[0], rmdupfq[0])
    read_utils.rmdup_prinseq_fastq(trimfq[1], rmdupfq[1])
    os.unlink(trimfq[0])
    os.unlink(trimfq[1])
    
    # Purge unmated
    purgefq = list(map(util.file.mkstempfname, ['.fix.1.fastq', '.fix.2.fastq']))
    read_utils.purge_unmated(rmdupfq[0], rmdupfq[1], purgefq[0], purgefq[1])
    os.unlink(rmdupfq[0])
    os.unlink(rmdupfq[1])

    # Log count
    with open(purgefq[0], 'rt') as inf:
        n = sum(1 for line in inf if line.strip()=='+')
        log.info("PRE-SUBSAMPLE COUNT: {} read pairs".format(n))
    
    # Subsample
    subsampfq = list(map(util.file.mkstempfname, ['.subsamp.1.fastq', '.subsamp.2.fastq']))
    cmd = [os.path.join(util.file.get_scripts_path(), 'subsampler.py'),
        '-n', str(n_reads),
        '-mode', 'p',
        '-in', purgefq[0], purgefq[1],
        '-out', subsampfq[0], subsampfq[1],
        ]
    subprocess.check_call(cmd)
    os.unlink(purgefq[0])
    os.unlink(purgefq[1])
    
    # Fastq -> BAM
    # Note: this destroys RG IDs! We should instead frun the BAM->fastq step in a way
    # breaks out the read groups and perform the above steps in a way that preserves
    # the RG IDs.
    tmp_bam = util.file.mkstempfname('.subsamp.bam')
    tmp_header = util.file.mkstempfname('.header.sam')
    tools.picard.FastqToSamTool().execute(
        subsampfq[0], subsampfq[1], 'Dummy', tmp_bam)
    tools.samtools.SamtoolsTool().dumpHeader(inBam, tmp_header)
    tools.samtools.SamtoolsTool().reheader(tmp_bam, tmp_header, outBam)
    os.unlink(tmp_bam)
    os.unlink(tmp_header)
    os.unlink(subsampfq[0])
    os.unlink(subsampfq[1])
def parser_trim_rmdup_subsamp(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam',
        help='Input reads, unaligned BAM format.')
    parser.add_argument('clipDb',
        help='Trimmomatic clip DB.')
    parser.add_argument('outBam',
        help='Output reads, unaligned BAM format (currently, read groups and other header information are destroyed in this process).')
    parser.add_argument('--n_reads', default=100000, type=int,
        help='Subsample reads to no more than this many pairs. (default %(default)s)')
    util.cmd.common_args(parser, (('loglevel',None), ('version',None), ('tmpDir',None)))
    util.cmd.attach_main(parser, trim_rmdup_subsamp_reads, split_args=True)
    return parser
__commands__.append(('trim_rmdup_subsamp', parser_trim_rmdup_subsamp))


def assemble_trinity(inBam, outFasta, clipDb, n_reads=100000, outReads=None):
    ''' This step runs the Trinity assembler.
        First trim reads with trimmomatic, rmdup with prinseq,
        and random subsample to no more than 100k reads.
    '''
    if outReads:
        subsamp_bam = outReads
    else:
        subsamp_bam = util.file.mkstempfname('.subsamp.bam')
    
    trim_rmdup_subsamp_reads(inBam, clipDb, subsamp_bam, n_reads=n_reads)
    subsampfq = list(map(util.file.mkstempfname, ['.subsamp.1.fastq', '.subsamp.2.fastq']))
    tools.picard.SamToFastqTool().execute(subsamp_bam, subsampfq[0], subsampfq[1])
    tools.trinity.TrinityTool().execute(subsampfq[0], subsampfq[1], outFasta)
    os.unlink(subsampfq[0])
    os.unlink(subsampfq[1])
    
    if not outReads:
        os.unlink(subsamp_bam)

def parser_assemble_trinity(parser=argparse.ArgumentParser()):
    parser.add_argument('inBam',
        help='Input reads, BAM format.')
    parser.add_argument('clipDb',
        help='Trimmomatic clip DB.')
    parser.add_argument('outFasta',
        help='Output assembly.')
    parser.add_argument('--n_reads', default=100000, type=int,
        help='Subsample reads to no more than this many pairs. (default %(default)s)')
    parser.add_argument('--outReads', default=None,
        help='Save the trimmomatic/prinseq/subsamp reads to a BAM file')
    util.cmd.common_args(parser, (('loglevel',None), ('version',None), ('tmpDir',None)))
    util.cmd.attach_main(parser, assemble_trinity, split_args=True)
    return parser
__commands__.append(('assemble_trinity', parser_assemble_trinity))


def order_and_orient(inFasta, inReference, outFasta, inReads=None):
    ''' This step cleans up the Trinity assembly with a known reference genome.
        Uses VFAT (switch to Bellini later):
        Take the Trinity contigs, align them to the known reference genome,
        switch it to the same strand as the reference, and produce
        chromosome-level assemblies (with runs of N's in between the Trinity
        contigs).
    '''
    # VFAT to order, orient, and merge contigs
    # TO DO: replace with Bellini
    musclepath = tools.muscle.MuscleTool().install_and_get_path()
    tmp_prefix = util.file.mkstempfname(prefix='VFAT-')
    cmd = [os.path.join(util.file.get_scripts_path(), 'vfat', 'orientContig.pl'),
        inFasta, inReference, tmp_prefix,
        '-musclepath', musclepath]
    subprocess.check_call(cmd)
    cmd = [os.path.join(util.file.get_scripts_path(), 'vfat', 'contigMerger.pl'),
        tmp_prefix+'_orientedContigs', inReference, tmp_prefix,
        '-musclepath', musclepath,
        '-samtoolspath', tools.samtools.SamtoolsTool().install_and_get_path()]
    if inReads:
        infq = list(map(util.file.mkstempfname, ['.in.1.fastq', '.in.2.fastq']))
        tools.picard.SamToFastqTool().execute(inReads, infq[0], infq[1])
        mosaik = tools.mosaik.MosaikTool()
        cmd = cmd + [
            '-readfq', infq[0], '-readfq2', infq[1],
            '-mosaikpath', os.path.dirname(mosaik.install_and_get_path()),
            '-mosaiknetworkpath', mosaik.get_networkFile(),
        ]
    subprocess.check_call(cmd)
    shutil.move(tmp_prefix+'_assembly.fa', outFasta)
    for fn in glob.glob(tmp_prefix+'*'):
        os.unlink(fn)
    with open(outFasta, 'rt') as inf:
        out_chr_count = len([1 for x in inf if x.startswith('>')])
    with open(inReference, 'rt') as inf:
        ref_chr_count = len([1 for x in inf if x.startswith('>')])
    if out_chr_count != ref_chr_count:
        raise Exception("error: expected {} chromosomes, only got {} chromosomes".format(ref_chr_count, out_chr_count))
    return 0

def parser_order_and_orient(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta',
        help='Input assembly/contigs, FASTA format.')
    parser.add_argument('inReference',
        help='Reference genome, FASTA format.')
    parser.add_argument('outFasta',
        help='Output assembly, FASTA format.')
    parser.add_argument('--inReads', default=None,
        help='Input reads in BAM format.')
    util.cmd.common_args(parser, (('loglevel',None), ('version',None), ('tmpDir',None)))
    util.cmd.attach_main(parser, order_and_orient, split_args=True)
    return parser
__commands__.append(('order_and_orient', parser_order_and_orient))


class PoorAssemblyError(Exception):
    def __init__(self, chr_idx, seq_len, non_n_count):
        super(PoorAssemblyError, self).__init__(
            'Error: poor assembly quality, chr {}: contig length {}, unambiguous bases {}'.format(
            chr_idx, seq_len, non_n_count))

def impute_from_reference(inFasta, inReference, outFasta,
        minLength, minUnambig, replaceLength, newName=None):
    '''
        This takes a de novo assembly, aligns against a reference genome, and
        imputes all missing positions (plus some of the chromosome ends)
        with the reference genome. This provides an assembly with the proper
        structure (but potentially wrong sequences in areas) from which
        we can perform further read-based refinement.
        Two steps:
        filter_short_seqs: We then toss out all assemblies that come out to
            < 15kb or < 95% unambiguous and fail otherwise.
        modify_contig: Finally, we trim off anything at the end that exceeds
            the length of the known reference assembly.  We also replace all
            Ns and everything within 55bp of the chromosome ends with the
            reference sequence.  This is clearly incorrect consensus sequence,
            but it allows downstream steps to map reads in parts of the genome
            that would otherwise be Ns, and we will correct all of the inferred
            positions with two steps of read-based refinement (below), and
            revert positions back to Ns where read support is lacking.
        FASTA indexing: output assembly is indexed for Picard, Samtools, Novoalign.
    '''
    # Halt if the assembly looks too poor at this point
    # TO DO: this can easily be generalized to multi-chr genomes by changing minLength to a fraction
    chr_idx = 0
    for seq in Bio.SeqIO.parse(inFasta, 'fasta'):
        chr_idx += 1
        non_n_count = unambig_count(seq.seq)
        seq_len = len(seq)
        if seq_len<minLength or non_n_count<seq_len*minUnambig:
            raise PoorAssemblyError(chr_idx, seq_len, non_n_count)
    
    # Align to known reference and impute missing sequences
    # TO DO: this can be iterated per chromosome
    concat_file = util.file.mkstempfname('.ref_and_actual.fasta')
    muscle_align = util.file.mkstempfname('.muscle.fasta')
    refName = None
    with open(concat_file, 'wt') as outf:
        with open(inReference, 'rt') as inf:
            for line in inf:
                if not refName and line.startswith('>'):
                    refName = line[1:].strip()
                outf.write(line)
        with open(inFasta, 'rt') as inf:
            for line in inf:
                outf.write(line)
    tools.muscle.MuscleTool().execute(concat_file, muscle_align, quiet=True)
    args = [muscle_align, outFasta, refName,
        '--call-reference-ns', '--trim-ends',
        '--replace-5ends', '--replace-3ends',
        '--replace-length', str(replaceLength),
        '--replace-end-gaps']
    if newName:
        args = args + ['--name', newName]
    args = parser_modify_contig().parse_args(args)
    args.func_main(args)
    os.unlink(concat_file)
    os.unlink(muscle_align)
    
    # Index final output FASTA for Picard/GATK, Samtools, and Novoalign
    tools.picard.CreateSequenceDictionaryTool().execute(outFasta, overwrite=True)
    tools.samtools.SamtoolsTool().faidx(outFasta, overwrite=True)
    tools.novoalign.NovoalignTool().index_fasta(outFasta)
    
    return 0

def parser_impute_from_reference(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta',
        help='Input assembly/contigs, FASTA format.')
    parser.add_argument('inReference',
        help='Reference genome, FASTA format.')
    parser.add_argument('outFasta',
        help='Output assembly, FASTA format.')
    parser.add_argument("--newName", default=None,
        help="rename output chromosome (default: do not rename)")
    parser.add_argument("--minLength", type=int, default=0,
        help="minimum length for contig (default: %(default)s)")
    parser.add_argument("--minUnambig", type=float, default=0.0,
        help="minimum percentage unambiguous bases for contig (default: %(default)s)")
    parser.add_argument("--replaceLength", type=int, default=0,
        help="length of ends to be replaced with reference (default: %(default)s)")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None), ('tmpDir',None)))
    util.cmd.attach_main(parser, impute_from_reference, split_args=True)
    return parser
__commands__.append(('impute_from_reference', parser_impute_from_reference))


def refine_assembly(inFasta, inBam, outFasta,
        outVcf=None, outBam=None, novo_params='', min_coverage=2,
        chr_names=[], keep_all_reads=False, JVMmemory=None):
    ''' This a refinement step where we take a crude assembly, align
        all reads back to it, and modify the assembly to the majority
        allele at each position based on read pileups.
        This step considers both SNPs as well as indels called by GATK
        and will correct the consensus based on GATK calls.
        Reads are aligned with Novoalign, then PCR duplicates are removed
        with Picard (in order to debias the allele counts in the pileups),
        and realigned with GATK's IndelRealigner (in order to call indels).
        Output FASTA file is indexed for Picard, Samtools, and Novoalign.
    '''
    # Get tools
    picard_index = tools.picard.CreateSequenceDictionaryTool()
    picard_mkdup = tools.picard.MarkDuplicatesTool()
    samtools = tools.samtools.SamtoolsTool()
    novoalign = tools.novoalign.NovoalignTool()
    gatk = tools.gatk.GATKTool()
    
    # Create deambiguated genome for GATK
    deambigFasta = util.file.mkstempfname('.deambig.fasta')
    deambig_fasta(inFasta, deambigFasta)
    picard_index.execute(deambigFasta, overwrite=True)
    samtools.faidx(deambigFasta, overwrite=True)
    
    # Novoalign reads to self
    novoBam = util.file.mkstempfname('.novoalign.bam')
    min_qual = 0 if keep_all_reads else 1
    novoalign.execute(inBam, inFasta, novoBam,
        options=novo_params.split(), min_qual=min_qual, JVMmemory=JVMmemory)
    rmdupBam = util.file.mkstempfname('.rmdup.bam')
    opts = ['CREATE_INDEX=true']
    if not keep_all_reads:
        opts.append('REMOVE_DUPLICATES=true')
    picard_mkdup.execute([novoBam], rmdupBam,
        picardOptions=opts, JVMmemory=JVMmemory)
    os.unlink(novoBam)
    realignBam = util.file.mkstempfname('.realign.bam')
    gatk.local_realign(rmdupBam, deambigFasta, realignBam, JVMmemory=JVMmemory)
    os.unlink(rmdupBam)
    if outBam:
        shutil.copyfile(realignBam, outBam)
    
    # Modify original assembly with VCF calls from GATK
    tmpVcf = util.file.mkstempfname('.vcf.gz')
    gatk.ug(realignBam, deambigFasta, tmpVcf, JVMmemory=JVMmemory)
    os.unlink(realignBam)
    os.unlink(deambigFasta)
    name_opts = []
    if chr_names:
        name_opts = ['--name'] + chr_names
    main_vcf_to_fasta(parser_vcf_to_fasta().parse_args([
        tmpVcf, outFasta, '--trim_ends', '--min_coverage', str(min_coverage),
        ] + name_opts))
    if outVcf:
        shutil.copyfile(tmpVcf, outVcf)
        if outVcf.endswith('.gz'):
            shutil.copyfile(tmpVcf+'.tbi', outVcf+'.tbi')
    os.unlink(tmpVcf)
    
    # Index final output FASTA for Picard/GATK, Samtools, and Novoalign
    picard_index.execute(outFasta, overwrite=True)
    samtools.faidx(outFasta, overwrite=True)
    novoalign.index_fasta(outFasta)
    return 0

def parser_refine_assembly(parser=argparse.ArgumentParser()):
    parser.add_argument('inFasta',
        help='Input assembly, FASTA format, pre-indexed for Picard, Samtools, and Novoalign.')
    parser.add_argument('inBam',
        help='Input reads, BAM format.')
    parser.add_argument('outFasta',
        help='Output refined assembly, FASTA format, indexed for Picard, Samtools, and Novoalign.')
    parser.add_argument('--outBam',
        default=None,
        help='Reads aligned to inFasta. Unaligned and duplicate reads have been removed. GATK indel realigned.')
    parser.add_argument('--outVcf',
        default=None,
        help='GATK genotype calls for genome in inFasta coordinate space.')
    parser.add_argument('--min_coverage',
        default=3, type=int,
        help='Minimum read coverage required to call a position unambiguous.')
    parser.add_argument('--novo_params',
        default='-r Random -l 40 -g 40 -x 20 -t 100',
        help='Alignment parameters for Novoalign.')
    parser.add_argument("--chr_names", dest="chr_names", nargs="*",
        help="Rename all output chromosomes (default: retain original chromosome names)",
        default=[])
    parser.add_argument("--keep_all_reads",
        help="""Retain all reads in BAM file? Default is to remove unaligned and duplicate reads.""",
        default=False, action="store_true", dest="keep_all_reads")
    parser.add_argument('--JVMmemory',
        default=tools.gatk.GATKTool.jvmMemDefault,
        help='JVM virtual memory size (default: %(default)s)')
    util.cmd.common_args(parser, (('loglevel',None), ('version',None), ('tmpDir',None)))
    util.cmd.attach_main(parser, refine_assembly, split_args=True)
    return parser
__commands__.append(('refine_assembly', parser_refine_assembly))




def unambig_count(seq):
    unambig = set(('A','T','C','G'))
    return sum(1 for s in seq if s.upper() in unambig)

def parser_filter_short_seqs(parser=argparse.ArgumentParser()):
    parser.add_argument("inFile", help="input sequence file")
    parser.add_argument("minLength", help="minimum length for contig", type=int)
    parser.add_argument("minUnambig", help="minimum percentage unambiguous bases for contig", type=float)
    parser.add_argument("outFile", help="output file")
    parser.add_argument("-f", "--format", help="Format for input sequence (default: %(default)s)", default="fasta")
    parser.add_argument("-of", "--output-format",
                        help="Format for output sequence (default: %(default)s)", default="fasta")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_filter_short_seqs)
    return parser
def main_filter_short_seqs(args):
    '''Check sequences in inFile, retaining only those that are at least minLength'''
    # orig by rsealfon, edited by dpark
    # TO DO: make this more generalized to accept multiple minLengths (for multiple chromosomes/segments)
    with util.file.open_or_gzopen(args.inFile) as inf:
        with util.file.open_or_gzopen(args.outFile, 'w') as outf:
            Bio.SeqIO.write(
                [s for s in Bio.SeqIO.parse(inf, args.format)
                    if len(s) >= args.minLength
                    and unambig_count(s.seq) >= len(s)*args.minUnambig],
                outf, args.output_format)
    return 0
__commands__.append(('filter_short_seqs', parser_filter_short_seqs))



def parser_modify_contig(parser=argparse.ArgumentParser()):
    parser.add_argument("input", help="input alignment of reference and contig (should contain exactly 2 sequences)")
    parser.add_argument("output", help="Destination file for modified contigs")
    parser.add_argument("ref", help="reference sequence name (exact match required)")
    parser.add_argument("-n", "--name",
        help="fasta header output name (default: existing header)",
        default=None)
    parser.add_argument("-cn", "--call-reference-ns",
        help="""should the reference sequence be called if there is an
        N in the contig and a more specific base in the reference (default: %(default)s)""",
        default=False, action="store_true", dest="call_reference_ns")
    parser.add_argument("-t", "--trim-ends",
        help="should ends of contig.fasta be trimmed to length of reference (default: %(default)s)",
        default=False, action="store_true", dest="trim_ends")
    parser.add_argument("-r5", "--replace-5ends",
        help="should the 5'-end of contig.fasta be replaced by reference (default: %(default)s)",
        default=False, action="store_true", dest="replace_5ends")
    parser.add_argument("-r3", "--replace-3ends",
        help="should the 3'-end of contig.fasta be replaced by reference (default: %(default)s)",
        default=False, action="store_true", dest="replace_3ends")
    parser.add_argument("-l", "--replace-length",
        help="length of ends to be replaced (if replace-ends is yes) (default: %(default)s)",
        default=10, type=int)
    parser.add_argument("-f", "--format",
        help="Format for input alignment (default: %(default)s)",
        default="fasta")
    parser.add_argument("-r", "--replace-end-gaps",
        help="Replace gaps at the beginning and end of the sequence with reference sequence (default: %(default)s)",
        default=False, action="store_true", dest="replace_end_gaps")
    parser.add_argument("-rn", "--remove-end-ns",
        help="Remove leading and trailing N's in the contig (default: %(default)s)",
        default=False, action="store_true", dest="remove_end_ns")
    parser.add_argument("-ca", "--call-reference-ambiguous",
        help="""should the reference sequence be called if the contig seq is ambiguous and
        the reference sequence is more informative & consistant with the ambiguous base
        (ie Y->C) (default: %(default)s)""",
        default=False, action="store_true", dest="call_reference_ambiguous")
    util.cmd.common_args(parser, (('tmpDir',None), ('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_modify_contig)
    return parser
def main_modify_contig(args):
    ''' Modifies an input contig. Depending on the options
        selected, can replace N calls with reference calls, replace ambiguous
        calls with reference calls, trim to the length of the reference, replace
        contig ends with reference calls, and trim leading and trailing Ns.
        Author: rsealfon.
    '''
    aln = Bio.AlignIO.read(args.input, args.format)
    
    if len(aln) != 2:
        raise Exception("alignment does not contain exactly 2 sequences")
    elif aln[0].name == args.ref:
        ref_idx = 0
        consensus_idx = 1
    elif aln[1].name == args.ref:
        ref_idx = 1
        consensus_idx = 0
    else:
        raise NameError("reference name '%s' not in alignment" % args.ref)
    
    mc = ContigModifier(str(aln[ref_idx].seq), str(aln[consensus_idx].seq))
    if args.remove_end_ns:
        mc.remove_end_ns()
    if args.call_reference_ns:
        mc.call_reference_ns()
    if args.call_reference_ambiguous:
        mc.call_reference_ambiguous()
    if args.trim_ends:
        mc.trim_ends()
    if args.replace_end_gaps:
        mc.replace_end_gaps()
    if args.replace_5ends:
        mc.replace_5ends(args.replace_length)
    if args.replace_3ends:
        mc.replace_3ends(args.replace_length)
    
    with open(args.output, "wt") as f:
        name = args.name and args.name or aln[consensus_idx].name
        for line in util.file.fastaMaker([(name, mc.get_stripped_consensus())]):
            f.write(line)
    return 0
__commands__.append(('modify_contig', parser_modify_contig))


class ContigModifier:
    ''' Initial modifications to Trinity+VFAT assembly output based on
        MUSCLE alignment to known reference genome
        author: rsealfon
    '''
    def __init__(self, ref, consensus):
        if len(ref) != len(consensus):
            raise Exception("improper alignment")
        self.ref       = list(ref)
        self.consensus = list(consensus)
        self.len = len(ref)
    
    def get_stripped_consensus(self):
        return ''.join(self.consensus).replace('-','')

    def call_reference_ns(self):
        log.debug("populating N's from reference...")
        for i in range(self.len):
            if self.consensus[i].upper() == "N":
                self.consensus[i] = self.ref[i]

    def call_reference_ambiguous(self):
        ''' This is not normally used by default in our pipeline '''
        log.debug("populating ambiguous bases from reference...")
        for i in range(self.len):
            if self.ref[i].upper() in Bio.Data.IUPACData.ambiguous_dna_values.get(self.consensus[i].upper(), []):
                self.consensus[i] = self.ref[i]

    def trim_ends(self):
        ''' This trims down the consensus so it cannot go beyond the given reference genome '''
        log.debug("trimming ends...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.ref[i] != "-":
                    break
                else:
                    self.consensus[i] = "-"

    def replace_end_gaps(self):
        ''' This fills out the ends of the consensus with reference sequence '''
        log.debug("populating leading and trailing gaps from reference...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if self.consensus[i] != "-":
                    break
                self.consensus[i] = self.ref[i]

    def replace_5ends(self, replace_length):
        ''' This replaces everything within <replace_length> of the ends of the
            reference genome with the reference genome.
        '''
        log.debug("replacing 5' ends...")
        ct = 0
        for i in range(self.len):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i+1):
                    self.consensus[j] = self.ref[j]
                break

    def replace_3ends(self, replace_length):
        log.debug("replacing 3' ends...")
        ct = 0
        for i in reversed(range(self.len)):
            if self.ref[i] != "-":
                ct = ct + 1
            if ct == replace_length:
                for j in range(i, self.len):
                    self.consensus[j] = self.ref[j]
                break

    def remove_end_ns(self):
        ''' This clips off any N's that begin or end the consensus.
            Not normally used in our pipeline
        '''
        log.debug("removing leading and trailing N's...")
        for end_iterator in (range(self.len), reversed(range(self.len))):
            for i in end_iterator:
                if (self.consensus[i].upper() == "N" or self.consensus[i] == "-"):
                    self.consensus[i] = "-"
                else:
                    break



class MutableSequence:
    def __init__(self, name, start, stop, init_seq=None):
        if not (stop>=start>=1):
            raise IndexError("coords out of bounds")
        if init_seq==None:
            self.seq = list('N' * (stop-start+1))
        else:
            self.seq = list(init_seq)
        if stop-start+1 != len(self.seq):
            raise Exception("wrong length")
        self.start = start
        self.stop = stop
        self.name = name
        self.deletions = []
    def modify(self, p, new_base):
        if not (self.start <= p <= self.stop):
            raise IndexError("position out of bounds")
        i = p-self.start
        self.seq[i] = new_base
    def replace(self, start, stop, new_seq):
        if stop>start:
            self.deletions.append((start, stop, new_seq))
        self.__change__(start, stop, new_seq)
    def __change__(self, start, stop, new_seq):
        if not (self.start <= start <= stop <= self.stop):
            raise IndexError("positions out of bounds")
        start -= self.start
        stop  -= self.start
        if start==stop:
            self.seq[start] = new_seq
        for i in range(max(stop-start+1, len(new_seq))):
            if start+i <= stop:
                if i<len(new_seq):
                    if start+i==stop:
                        # new allele is >= ref length, fill out the rest of the bases
                        self.seq[start+i] = new_seq[i:]
                    else:
                        self.seq[start+i] = new_seq[i]
                else:
                    # new allele is shorter than ref, so delete extra bases
                    self.seq[start+i] = ''
    def replay_deletions(self):
        for start, stop, new_seq in self.deletions:
            self.__change__(start, stop, new_seq)
    def emit(self):
        return (self.name, ''.join(self.seq))

def alleles_to_ambiguity(allelelist):
    ''' Convert a list of DNA bases to a single ambiguity base.
        All alleles must be one base long.  '''
    for a in allelelist:
        if len(a)!=1:
            raise Exception("all alleles must be one base long")
    if len(allelelist)==1:
        return allelelist[0]
    else:
        convert = dict([(tuple(sorted(v)),k) for k,v in Bio.Data.IUPACData.ambiguous_dna_values.items() if k!='X'])
        key = tuple(sorted(set(a.upper() for a in allelelist)))
        return convert[key]

def vcfrow_parse_and_call_snps(vcfrow, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
    ''' Parse a single row of a VCF file, emit an iterator over each sample,
        call SNP genotypes using custom viral method based on read counts.
    '''
    # process this row
    c = vcfrow[0]
    alleles = [vcfrow[3]] + [a for a in vcfrow[4].split(',') if a not in '.']
    start = int(vcfrow[1])
    stop = start + len(vcfrow[3]) - 1
    format = vcfrow[8].split(':')
    format = dict((format[i], i) for i in range(len(format)))
    assert 'GT' in format and format['GT']==0  # required by VCF spec
    assert len(vcfrow)==9+len(samples)
    info = [x.split('=') for x in vcfrow[7].split(';') if x != '.']
    info = dict(x for x in info if len(x)==2)
    info_dp = int(info.get('DP',0))

    # process each sample
    for i in range(len(samples)):
        sample = samples[i]
        rec = vcfrow[i+9].split(':')
        # require a minimum read coverage
        if len(alleles)==1:
            # simple invariant case
            dp = ('DP' in format and len(rec)>format['DP']) and int(rec[format['DP']]) or 0
            if dp < min_dp:
                continue
            geno = alleles
            if info_dp and float(dp)/info_dp < min_dp_ratio:
                log.warn("dropping invariant call at %s:%s %s (%s) due to low DP ratio (%s / %s = %s < %s)" % (
                    c,p,sample,geno,dp,info_dp,float(dp)/info_dp,min_dp_ratio))
                continue
        else:
            # variant: manually call the highest read count allele if it exceeds a threshold
            assert ('AD' in format and len(rec)>format['AD'])
            allele_depths = list(map(int, rec[format['AD']].split(',')))
            assert len(allele_depths)==len(alleles)
            allele_depths = [(allele_depths[i], alleles[i]) for i in range(len(alleles)) if allele_depths[i]>0]
            allele_depths = list(reversed(sorted((n,a) for n,a in allele_depths if n>=min_dp)))
            if not allele_depths:
                continue
            dp = sum(n for n,a in allele_depths)

            if allele_depths[0][0] > (dp*major_cutoff):
                # call a single allele at this position if it is a clear winner
                geno = [allele_depths[0][1]]
            else:
                # call multiple alleles at this position if there is no clear winner
                geno = [a for n,a in allele_depths]
        if geno:
            yield (c, start, stop, sample, geno)

def vcf_to_seqs(vcfIter, chrlens, samples, min_dp=0, major_cutoff=0.5, min_dp_ratio=0.0):
    ''' Take a VCF iterator and produce an iterator of chromosome x sample full sequences.'''
    seqs = {}
    cur_c = None
    for vcfrow in vcfIter:
        try:
            for c,start,stop,s,alleles in vcfrow_parse_and_call_snps(vcfrow, samples, min_dp=min_dp, major_cutoff=major_cutoff, min_dp_ratio=min_dp_ratio):
                # changing chromosome?
                if c != cur_c:
                    if cur_c!=None:
                        # dump the previous chromosome before starting a new one
                        for s in samples:
                            seqs[s].replay_deletions() # because of the order of VCF rows with indels
                            yield seqs[s].emit()

                    # prepare base sequences for this chromosome
                    cur_c = c
                    for s in samples:
                        name = len(samples)>1 and ("%s-%s" % (c,sample)) or c
                        seqs[s] = MutableSequence(name, 1, chrlens[c])

                # modify sequence for this chromosome/sample/position
                if len(alleles)==1:
                    # call a single allele
                    seqs[s].replace(start, stop, alleles[0])
                elif all(len(a)==1 for a in alleles):
                    # call an ambiguous SNP
                    seqs[s].replace(start, stop, alleles_to_ambiguity(alleles))
                else:
                    # mix of indels with no clear winner... force the most popular one
                    seqs[s].replace(start, stop, alleles[0])
        except:
            log.exception("Exception occurred while parsing VCF file.  Row: '%s'" % vcfrow)
            raise

    # at the end, dump the last chromosome
    if cur_c!=None:
        for s in samples:
            seqs[s].replay_deletions() # because of the order of VCF rows with indels
            yield seqs[s].emit()


def parser_vcf_to_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcf", help="Input VCF file")
    parser.add_argument("outFasta", help="Output FASTA file")
    parser.add_argument("--trim_ends",
        action="store_true", dest="trim_ends",
        default=False,
        help="""If specified, we will strip off continuous runs of N's from the beginning
        and end of the sequences before writing to output.  Interior N's will not be
        changed.""")
    parser.add_argument("--min_coverage", dest="min_dp", type=int,
        help="""Specify minimum read coverage (with full agreement) to make a call.
        [default: %(default)s]""",
        default=3)
    parser.add_argument("--major_cutoff", dest="major_cutoff", type=float,
        help="""If the major allele is present at a frequency higher than this cutoff,
        we will call an unambiguous base at that position.  If it is equal to or below
        this cutoff, we will call an ambiguous base representing all possible alleles at
        that position. [default: %(default)s]""",
        default=0.5)
    parser.add_argument("--min_dp_ratio", dest="min_dp_ratio", type=float,
        help="""The input VCF file often reports two read depth values (DP)--one for
        the position as a whole, and one for the sample in question.  We can optionally
        reject calls in which the sample read count is below a specified fraction of the
        total read count.  This filter will not apply to any sites unless both DP values
        are reported.  [default: %(default)s]""",
        default=0.0)
    parser.add_argument("--name", dest="name", nargs="*",
        help="output sequence names (default: reference names in VCF file)",
        default=[])
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, main_vcf_to_fasta)
    return parser
def main_vcf_to_fasta(args):
    ''' Take input genotypes (VCF) and construct a consensus sequence
        (fasta) by using majority-read-count alleles in the VCF.
        Genotypes in the VCF will be ignored--we will use the allele
        with majority read support (or an ambiguity base if there is no clear majority).
        Uncalled positions will be emitted as N's.
        Author: dpark.
    '''
    assert args.min_dp >= 0
    assert 0.0 <= args.major_cutoff < 1.0
    
    with util.vcf.VcfReader(args.inVcf) as vcf:
        chrlens = dict(vcf.chrlens())
        samples = vcf.samples()
    with open(args.outFasta, 'wt') as outf:
        chr_idx = 0
        for header, seq in vcf_to_seqs(util.file.read_tabfile(args.inVcf),
            chrlens, samples, min_dp=args.min_dp, major_cutoff=args.major_cutoff,
            min_dp_ratio=args.min_dp_ratio):
            if args.trim_ends:
                seq = seq.strip('Nn')
            if args.name:
                header = args.name[chr_idx % len(args.name)]
            for line in util.file.fastaMaker([(header, seq)]):
                outf.write(line)

    # done
    log.info("done")
    return 0
__commands__.append(('vcf_to_fasta', parser_vcf_to_fasta))



def parser_trim_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inFasta", help="Input fasta file")
    parser.add_argument("outFasta", help="Output (trimmed) fasta file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, trim_fasta, split_args=True)
    return parser
def trim_fasta(inFasta, outFasta):
    ''' Take input sequences (fasta) and trim any continuous sections of
        N's from the ends of them.  Write trimmed sequences to an output fasta file.
    '''
    with open(outFasta, 'wt') as outf:
        with open(inFasta, 'rt') as inf:
            for record in Bio.SeqIO.parse(inf, 'fasta'):
                for line in util.file.fastaMaker([(record.id, str(record.seq).strip('Nn'))]):
                    outf.write(line)
    log.info("done")
    return 0
__commands__.append(('trim_fasta', parser_trim_fasta))



def deambig_base(base):
    ''' Take a single base (possibly a IUPAC ambiguity code) and return a random
        non-ambiguous base from among the possibilities '''
    return random.choice(Bio.Data.IUPACData.ambiguous_dna_values[base.upper()])

def deambig_fasta(inFasta, outFasta):
    ''' Take input sequences (fasta) and replace any ambiguity bases with a
        random unambiguous base from among the possibilities described by the ambiguity
        code.  Write output to fasta file.
    '''
    with util.file.open_or_gzopen(outFasta, 'wt') as outf:
        with util.file.open_or_gzopen(inFasta, 'rt') as inf:
            for record in Bio.SeqIO.parse(inf, 'fasta'):
                for line in util.file.fastaMaker([(record.id, ''.join(map(deambig_base, str(record.seq))))]):
                    outf.write(line)
    return 0
def parser_deambig_fasta(parser=argparse.ArgumentParser()):
    parser.add_argument("inFasta", help="Input fasta file")
    parser.add_argument("outFasta", help="Output fasta file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, deambig_fasta, split_args=True)
    return parser
__commands__.append(('deambig_fasta', parser_deambig_fasta))


def vcf_dpdiff(vcfs):
    for vcf in vcfs:
        with util.vcf.VcfReader(vcf) as v:
            samples = v.samples()
        assert len(samples)==1
        for row in util.file.read_tabfile(vcf):
            dp1 = int(dict(x.split('=') for x in row[7].split(';') if x != '.').get('DP',0))
            dp2 = 0
            if 'DP' in row[8].split(':'):
                dpidx = row[8].split(':').index('DP')
                if len(row[9].split(':'))>dpidx:
                    dp2 = int(row[9].split(':')[dpidx])
            ratio = ''
            if dp1:
                ratio = float(dp2)/dp1
            yield (row[0],row[1],samples[0],dp1,dp2,dp1-dp2,ratio)


def parser_dpdiff(parser=argparse.ArgumentParser()):
    parser.add_argument("inVcfs", help="Input VCF file", nargs='+')
    parser.add_argument("outFile", help="Output flat file")
    util.cmd.common_args(parser, (('loglevel',None), ('version',None)))
    util.cmd.attach_main(parser, dpdiff, split_args=True)
    return parser
def dpdiff(inVcfs, outFile):
    ''' Take input VCF files (all with only one sample each) and report
        on the discrepancies between the two DP fields (one in INFO and one in the
        sample's genotype column).
    '''
    header = ['chr','pos','sample','dp_info','dp_sample','diff','ratio']
    with open(outFile, 'wt') as outf:
        outf.write('#'+'\t'.join(header)+'\n')
        for row in vcf_dpdiff(inVcfs):
            outf.write('\t'.join(map(str, row))+'\n')
    return 0
__commands__.append(('dpdiff', parser_dpdiff))


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)
if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

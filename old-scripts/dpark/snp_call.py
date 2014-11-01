#!/usr/bin/env python
#
# snp_call.py - scripts for aligning reads and calling snps
#
# requires python >= 2.6 and python < 3
#
# $Id: snp_call.py 7891 2014-07-02 19:18:32Z dpark $

__author__ = "dpark@broadinstitute.org"
__version__ = "$Revision: 7891 $".split(' ')[1]
__date__ = "$Date: 2014-07-02 15:18:32 -0400 (Wed, 02 Jul 2014) $"

import glob, os, tempfile, shutil, optparse, logging, urllib
import util_cmd, util_files, util_vcf

log = logging.getLogger(__name__)

global_tool_paths = {
    'bedtools': '/idi/sabeti-data/software/bedtools/bin',
    'bwa':      '/idi/sabeti-data/software/bwa/bwa-0.6.2',
    #'gatk_old_sen45_paper': '/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-1.3-25-g32cdef9',
    'gatk':     '/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-2.4-9-g532efad',
    'gatk-new':     '/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-2.8-1-g932cd3a',
    'picard':   '/seq/software/picard/current/bin',
    'samtools': '/idi/sabeti-data/software/samtools/samtools-0.1.19',
    'snpEff':   '/idi/sabeti-data/software/snpEff/snpEff_3.6-dev',
    'tabix':    '/idi/sabeti-data/software/tabix/tabix-0.2.6',
    'vcftools': '/idi/sabeti-data/software/vcftools/vcftools_0.1.11/bin',
    'bass':     '/prodinfo/prodapps/dmsClient',
    }

def main():
    version="$Id: snp_call.py 7891 2014-07-02 19:18:32Z dpark $"[5:-2]
    description = """
This script contains a number of tools to run a number of steps in the general
SNP-discovery pipeline for NGS reads.  From handling, aligning, and filtering
BAM files, to using GATK to call variants, filter, etc.  Outputs are always
properly indexed and compressed."""
    commands = [
        ('realign', main_realign, parser_realign, 3),
        ('realign_from_file', main_realign_from_file, parser_realign_from_file, 4),
        ('sort_index_bam', main_sort_index_bam, parser_sort_index_bam, 2),
        ('index_bam', main_index_bam, parser_index_bam, 1),
        ('revert_bam', main_revert_bam, parser_revert_bam, 2),
        ('rename_bam', main_rename_bam, parser_rename_bam, 3),
        ('subset_bam', main_subset_bam, parser_subset_bam, 2),
        ('post_process_bam', main_post_process_bam, parser_post_process_bam, 2),
        ('gatk_ug', main_gatk_ug, parser_gatk_ug, None),
        ('fastqs_to_bam', main_fastqs_to_bam, parser_fastqs_to_bam, 1),
        ('ena_maketable', main_ena_maketable, parser_ena_maketable, 5),
        ('vcf_merge', main_vcf_merge, parser_vcf_merge, None),
        ('vcf_rename', main_vcf_rename, parser_vcf_rename, 3),
        ('vcf_subset', main_vcf_subset, parser_vcf_subset, 3),
        ('post_process_vcf', main_post_process_vcf, parser_post_process_vcf, 5),
        ('vcf_bgzip_index', main_vcf_bgzip_index, parser_vcf_bgzip_index, 2),
        ('vcf_hard_filter', main_vcf_hard_filter, parser_vcf_hard_filter, 2),
        ('snpEff', main_snpEff, parser_snpEff, 3),
        ('vcf_to_biallelic', main_vcf_to_biallelic, parser_vcf_to_biallelic, 3),
        ('filter_sample_callrate', main_filter_sample_callrate, parser_filter_sample_callrate, 3),
        ('sen_unique_patients', main_sen_unique_patients, parser_sen_unique_patients, 2),
        ('broad_bams', main_broad_bams, parser_broad_bams, 3),
        ('merge_bams_from_file', main_merge_bams_from_file, parser_merge_bams_from_file, 2),
    ]
    return util_cmd.main(commands, version, global_tool_paths, description=description)


def rename_bam_sample(inBam, newSampleName, outBam, samtoolsPath=global_tool_paths['samtools']):
    ''' Rewrite BAM header so that all read groups in the header are assigned to a new
        sample name.
    '''
    assert inBam!=outBam
    header_file1 = util_files.mkstempfname(prefix='rename_bam-header-1', suffix='.sam')
    header_file2 = util_files.mkstempfname(prefix='rename_bam-header-2', suffix='.sam')
    cmdline = "%s/samtools view -H %s > %s" % (samtoolsPath, inBam, header_file1)
    assert not os.system(cmdline)
    with open(header_file1, 'rt') as inf:
        with open(header_file2, 'wt') as outf:
            for row in inf:
                row = row.rstrip('\n\r').split('\t')
                if row[0].startswith('@RG'):
                    for i in range(len(row)):
                        if row[i].startswith('SM:'):
                            #log.info("replacing existing sample name %s with %s" % (row[i][3:], newSampleName))
                            row[i] = 'SM:' + newSampleName
                outf.write('\t'.join(row)+'\n')
    #if inBam==outBam:
    cmdline = "%s/samtools reheader %s %s > %s" % (samtoolsPath, header_file2, inBam, outBam)
    #else:
    #   cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/ReplaceSamHeader.jar INPUT=%s HEADER=%s OUTPUT=%s" % (
    #       tempfile.tempdir, picardPath, inBam, header_file2, outBam)
    log.info("rename_bam_sample (%s): %s" % (newSampleName, cmdline))
    assert not os.system(cmdline)
    os.unlink(header_file1)
    os.unlink(header_file2)
    return outBam

def realign_bam(inBam, ref, outBam, newSampleName=None, paired=None, revertBam=True, num_threads=1,
    bwaPath=global_tool_paths['bwa'], samtoolsPath=global_tool_paths['samtools'],
    picardPath=global_tool_paths['picard']):
    ''' Use BWA to realign reads in a BAM file to a new reference.
        If paired = 2, we assume paired-end reads.  If paired = 1, we assume single-end.
        If paired = None, we auto-detect (might take a couple minutes).
        Input ref fasta must already be decompressed and indexed by *both* BWA and samtools.
    '''
    log.info("realign_bam starting on %s -> %s" % (inBam, outBam))
    # strip all info from original bam file
    if revertBam:
        inBam_strip = util_files.mkstempfname(prefix='realign-revert-', suffix='.bam')
        revert_bam(inBam, inBam_strip, picardPath=picardPath)
    else:
        inBam_strip = inBam

    # find read groups and split input bam by read group
    rg_file = util_files.mkstempfname(prefix='realign-read_groups-', suffix='.txt')
    cmdline = "%s/samtools view -H %s | grep ^@RG > %s" % (samtoolsPath, inBam_strip, rg_file)
    assert not os.system(cmdline)
    rg_ids = []
    rg_full = {}
    with open(rg_file, 'rt') as inf:
        for row in inf:
            row = row.rstrip('\n\r').split('\t')
            rg_id = None
            for i in range(len(row)):
                if newSampleName and row[i].startswith('SM:'):
                    if not rg_ids:
                        log.info("replacing existing sample name %s with %s" % (row[i][3:], newSampleName))
                    row[i] = 'SM:' + newSampleName
                if row[i].startswith('ID:'):
                    rg_id = row[i][3:]
            assert rg_id
            rg_ids.append(rg_id)
            rg_full[rg_id] = '\t'.join(row)
            log.info("loading read group: %s" % rg_full[rg_id])
    os.unlink(rg_file)
    bams_in = {}
    assert rg_ids
    for rg_id in rg_ids:
        log.info("splitting input bam file by read group: %s" % rg_id)
        bam = '%s/realign-split-bam_in-%s.bam' % (tempfile.tempdir, rg_id)
        cmdline = "%s/samtools view -h -b -r %s %s > %s" % (samtoolsPath, rg_id, inBam_strip, bam)
        log.info("samtools split by RG: %s" % cmdline)
        assert not os.system(cmdline)
        bams_in[rg_id] = bam
    if revertBam:
        os.unlink(inBam_strip)

    bams_realigned = []
    rg_idx = 1
    for rg_id in rg_ids:
        bam_in = bams_in[rg_id]
        log.info("processing read group %d out of %d" % (rg_idx, len(rg_ids)))
        rg_idx += 1

        # detect pairing
        if not paired:
            # if we detect even one paired read alignment, we assume the whole thing is paired, which might not be true.
            ptest = util_files.mkstempfname(prefix='realign-ptest-', suffix='.sam')
            cmdline = "%s/samtools view -f 1 %s | head -1 > %s" % (samtoolsPath, bam_in, ptest)
            assert not os.system(cmdline)
            pairing = 1
            with open(ptest, 'rt') as inf:
                for row in inf:
                    if row.rstrip():
                        pairing = 2
                        break
            os.unlink(ptest)
            log.info("input bam file autodetected as %s-end reads" % (pairing==1 and 'single' or 'paired'))
        else:
            pairing = paired

        # convert input file to sai files
        saiFiles = []
        pairs = pairing>1 and (1,2) or (0,)
        for pair in pairs:
            saiFile = util_files.mkstempfname(prefix='realign-%s-'%rg_id, suffix='-%d.sai'%pair)
            saiFiles.append(saiFile)
            cmdline = "%s/bwa aln %s -b%d -t %d %s > %s" % (bwaPath, ref, pair, num_threads, bam_in, saiFile)
            log.info("bwa: %s" % cmdline)
            assert not os.system(cmdline)

        # convert sai and input bam files to sam file
        samFile = util_files.mkstempfname(prefix='realign-%s-'%rg_id, suffix='.sam')
        read_group = rg_full[rg_id].replace('\t','\\t')
        if len(saiFiles)==1:
            cmdline = "%s/bwa samse -r '%s' %s %s %s > %s" % (
                bwaPath, read_group, ref, saiFiles[0], bam_in, samFile)
        else:
            cmdline = "%s/bwa sampe -r '%s' %s %s %s %s %s > %s" % (
                bwaPath, read_group, ref, saiFiles[0], saiFiles[1], bam_in, bam_in, samFile)
        log.info("bwa: %s" % cmdline)
        assert not os.system(cmdline)
        for t in saiFiles:
            os.unlink(t)
        os.unlink(bam_in)

        # convert sam file to bam file
        samFile2 = util_files.mkstempfname(prefix='realign-clean-%s-'%rg_id, suffix='.sam')
        clean_sam(samFile, samFile2, picardPath=picardPath)
        os.unlink(samFile)
        bam_tmp = util_files.mkstempfname(prefix='realign-presort-%s-'%rg_id, suffix='.bam')
        sam_to_bam(samFile2, bam_tmp, picardPath=picardPath)
        os.unlink(samFile2)
        bam_out = '%s/realign-split-bam_out-%s.bam' % (tempfile.tempdir, rg_id)
        sort_bam(bam_tmp, bam_out, picardPath=picardPath)
        index_bam(bam_out, picardPath=picardPath)
        os.unlink(bam_tmp)
        bams_realigned.append(bam_out)

    # merge, clean up, and exit
    log.info("done with per-lane alignments, merging results:")
    merge_bams(bams_realigned, outBam, picardPath=picardPath)
    for bam in bams_realigned:
        os.unlink(bam)
        os.unlink(bam[:-4]+'.bai')
    index_bam(outBam, picardPath=picardPath)
    return outBam

def clean_sam(inSam, outSam, picardPath=global_tool_paths['picard']):
    ''' Clean sam file.  Right now this only soft-clips an alignment that hangs off
        the end of the reference sequence and sets MAPQ to 0 for unmapped reads.
    '''
    assert inSam.endswith('.sam') and outSam.endswith('.sam')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/CleanSam.jar INPUT=%s OUTPUT=%s QUIET=true VALIDATION_STRINGENCY=SILENT" % (
        tempfile.tempdir, picardPath, inSam, outSam)
    log.info("clean_sam: %s" % cmdline)
    assert not os.system(cmdline)
    return

def revert_bam(inBam, outBam, picardPath=global_tool_paths['picard']):
    ''' Clean sam file.  Right now this only soft-clips an alignment that hangs off
        the end of the reference sequence and sets MAPQ to 0 for unmapped reads.
    '''
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/RevertSam.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=SILENT" % (
        tempfile.tempdir, picardPath, inBam, outBam)
    log.info("revert_bam: %s" % cmdline)
    assert not os.system(cmdline)
    return

def sam_to_bam(inSam, outBam, picardPath=global_tool_paths['picard']):
    ''' Index bam file. '''
    assert inSam.endswith('.sam') and outBam.endswith('.bam')
    #cmdline = "%s/samtools import %s.fai %s %s" % (samtoolsPath, refFasta, inSam, outBam)
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/SamFormatConverter.jar INPUT=%s OUTPUT=%s" % (
        tempfile.tempdir, picardPath, inSam, outBam)
    log.info("sam_to_bam: %s" % cmdline)
    assert not os.system(cmdline)
    return

def sort_bam(inBam, outBam=None, picardPath=global_tool_paths['picard']):
    ''' Take input bam file write a sorted bam file.
        If outBam=None or outBam=inBam, then replace input file with output file.
    '''
    replace=0
    if outBam==None or outBam==inBam:
        replace=1
        outBam = util_files.mkstempfname(prefix='sort-', suffix='.bam')
    assert outBam.endswith('.bam')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/SortSam.jar INPUT=%s OUTPUT=%s SORT_ORDER=coordinate" % (
        tempfile.tempdir, picardPath, inBam, outBam)
    log.info("sort_bam: %s" % cmdline)
    assert not os.system(cmdline)
    log.debug("bam file size (bytes): %d" % os.path.getsize(outBam))
    if replace:
        shutil.move(outBam, inBam)
    return

def index_bam(bam, picardPath=global_tool_paths['picard']):
    ''' Index bam file. '''
    assert bam.endswith('.bam')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/BuildBamIndex.jar INPUT=%s QUIET=true" % (
        tempfile.tempdir, picardPath, bam)
    log.info("index_bam: %s" % cmdline)
    assert not os.system(cmdline)
    return

def merge_bams(inBams, outBam, sort_order='coordinate', picardPath=global_tool_paths['picard']):
    ''' Merge multiple, sorted BAM files to one output, sorted BAM file. '''
    assert outBam.endswith('.bam')
    if len(inBams) == 1:
        log.info("only one bam file given as input, skipping merge and copying file instead")
        shutil.copy(inBams[0], outBam)
    else:
        cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/MergeSamFiles.jar %s OUTPUT=%s USE_THREADING=true MERGE_SEQUENCE_DICTIONARIES=true" % (
            tempfile.tempdir, picardPath,
            ' '.join(["INPUT="+b for b in inBams]),
            outBam)
        if sort_order != 'coordinate':
            assert sort_order in ('coordinate','queryname','unsorted')
            cmdline += ' SORT_ORDER=' + sort_order
        log.info("merge: %s" % cmdline)
        assert not os.system(cmdline)
        log.debug("bam file size (bytes): %d" % os.path.getsize(outBam))
    return

def markdup_bam(inBam, outBam=None, picardPath=global_tool_paths['picard']):
    ''' Mark duplicates in bam file '''
    if outBam==None or outBam==inBam:
        replace=1
        outBam = util_files.mkstempfname(prefix='markdup-', suffix='.bam')
        metrics = "%s.dedup.metrics" % inBam
    else:
        replace=0
        metrics = "%s.dedup.metrics" % outBam
    assert outBam.endswith('.bam')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/MarkDuplicates.jar INPUT=%s OUTPUT=%s METRICS_FILE=%s VALIDATION_STRINGENCY=SILENT AS=true" % (
        tempfile.tempdir, picardPath, inBam, outBam, metrics)
    assert not os.system(cmdline)
    log.debug("bam file size (bytes): %d" % os.path.getsize(outBam))
    if replace:
        shutil.move(outBam, inBam)
    return

def local_realign_bam(inBam, outBam, refFasta, knownSites=[], gatkPath=global_tool_paths['gatk']):
    ''' note: this process wants about 4GB of RAM! '''
    assert outBam.endswith('.bam')
    targetIntervals = util_files.mkstempfname(prefix='gatk_indel_realign-', suffix='.intervals')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -I %s -R %s -o %s %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inBam, refFasta, targetIntervals,
        ' '.join(["-known "+k for k in knownSites]))
    assert not os.system(cmdline)
    cmdline = "java -Xmx3g -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -I %s -R %s -targetIntervals %s -o %s %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inBam, refFasta, targetIntervals, outBam,
        ' '.join(["-known "+k for k in knownSites]))
    assert not os.system(cmdline)
    log.debug("bam file size (bytes): %d" % os.path.getsize(outBam))
    os.unlink(targetIntervals)
    return

def recal_bam(inBam, outBam, refFasta, knownSites=[], gatkPath=global_tool_paths['gatk']):
    ''' note: this process wants about 4GB of RAM! '''
    assert knownSites, "GATK docs say this is pointless without a knownSites list"
    assert outBam.endswith('.bam')
    recalReport = util_files.mkstempfname(prefix='gatk_bqsr-', suffix='.grp')
    cmdline = "java -Xmx3g -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -I %s -R %s -o %s %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inBam, refFasta, recalReport,
        ' '.join(["--knownSites "+k for k in knownSites]))
    assert not os.system(cmdline)
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T PrintReads -I %s -R %s -BQSR %s -o %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inBam, refFasta, recalReport, outBam)
    assert not os.system(cmdline)
    log.debug("bam file size (bytes): %d" % os.path.getsize(outBam))
    return

def create_seq_dict(fasta, picardPath=global_tool_paths['picard']):
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/CreateSequenceDictionary.jar REFERENCE=%s OUTPUT=%s" % (
        tempfile.tempdir, picardPath, fasta, "%s.dict" % fasta)
    assert not os.system(cmdline)
    return

def post_process_bam(inBam, outBam, refFasta, knownSites=[], sort=True,
    picardPath=global_tool_paths['picard'], gatkPath=global_tool_paths['gatk']):
    ''' Assume input is a newly created BAM file containing all reads for a given
        sample.  We will sort, index, mark duplicates, perform local realignment around
        indels, and possibly perform a base quality score recalibration (if given
        a knownSites file).
        inBam is allowed to equal outBam.  In that case, the original file will be
        overwritten.
    '''
    if sort:
        log.info("sorting and indexing BAM file with Picard")
        sorted_bam = util_files.mkstempfname(prefix='process_bam-sort-', suffix='.bam')
        sort_bam(inBam, sorted_bam, picardPath=picardPath)
        index_bam(sorted_bam, picardPath=picardPath)
    else:
        sorted_bam=inBam

    log.info("marking duplicate reads with Picard")
    mkdup_bam = util_files.mkstempfname(prefix='process_bam-markdup-', suffix='.bam')
    markdup_bam(sorted_bam, mkdup_bam, picardPath=picardPath)
    index_bam(mkdup_bam, picardPath=picardPath)
    if sort:
        os.unlink(sorted_bam)
        os.unlink(sorted_bam[:-4] + '.bai')

    log.info("performing local realignment around indels with GATK")
    if knownSites:
        localr_bam = util_files.mkstempfname(prefix='process_bam-local_realign-', suffix='.bam')
    else:
        localr_bam = outBam
    local_realign_bam(mkdup_bam, localr_bam, refFasta, knownSites=knownSites, gatkPath=gatkPath)
    index_bam(localr_bam, picardPath=picardPath)
    os.unlink(mkdup_bam)
    os.unlink(mkdup_bam[:-4] + '.bai')
    os.unlink(mkdup_bam+'.dedup.metrics')

    if knownSites:
        log.info("performing base quality recalibration with GATK")
        recal_bam(localr_bam, outBam, refFasta, knownSites=knownSites, gatkPath=gatkPath)
        index_bam(outBam, picardPath=picardPath)
        os.unlink(localr_bam)
        os.unlink(localr_bam[:-4] + '.bai')
    else:
        log.info("skipping base quality recalibration (no known sites file provided)")

    return outBam

def gatk_ug(inBams, inFasta, outVcf,
    glm='SNP', ploidy=2, num_threads=None, interval_idx=None, interval_total=None,
    gatkPath=global_tool_paths['gatk'], tabixPath=global_tool_paths['tabix'], vcftoolsPath=global_tool_paths['vcftools']):
    ''' Call SNPs with GATK '''
    assert outVcf.endswith('.vcf.gz') or outVcf.endswith('.vcf')

    # run GATK UG and write to temp file VCF
    if outVcf.endswith('.vcf'):
        tmpFile = outVcf
    else:
        tmpFile = util_files.mkstempfname(prefix='gatk_ug-', suffix='.vcf')
    # -A AlleleBalanceBySample -A DepthPerAlleleBySample
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T UnifiedGenotyper --output_mode EMIT_ALL_SITES -stand_emit_conf 0.0 -stand_call_conf 0.0 %s -R %s -o %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        ' '.join(['-I '+i for i in inBams]), inFasta, tmpFile)
    if ploidy != 2:
        cmdline += ' --sample_ploidy %d' % ploidy
    if glm != 'SNP':
        assert glm in ('SNP','INDEL','BOTH')
        cmdline += ' -glm %s' % glm
    if num_threads:
        cmdline += ' --num_threads %d' % num_threads
    if interval_idx!=None and interval_total:
        for c,start,stop in util_vcf.make_intervals(interval_idx, interval_total, inFasta):
            cmdline += ' -L %s:%d-%d' % (c,start,stop)
    assert not os.system(cmdline)

    # bgzip to final output file and index with tabix and vcftools
    if outVcf.endswith('.vcf.gz'):
        util_vcf.vcf_bgzip_index(tmpFile, outVcf, tabixPath=tabixPath, vcftoolsPath=vcftoolsPath)
        os.unlink(tmpFile)
    return

def bundle_fastqs_to_unaligned_bam(pairs_list, outBam, sample, picardPath=global_tool_paths['picard']):
    ''' Take a set of fastq files and lump them all into a single bam file that is
        unaligned to any reference.  This merges tons of files into one and also overlays
        read group information.
        The input pairs_list is a list of dicts.  Each dict corresponds to a single
        or paired fastq lane for this sample.  The following entries are mandatory:
            FASTQ               - input reads in fastq format
            QUALITY_FORMAT      - Solexa, Illumina, or Standard
                                    Solexa - pre-pipeline 1.3 scores (solexa scaling + 66)
                                    Illumina - pipeline >=1.3 scores (phred scaling + 64)
                                    Standard - phred scaling + 33
        The following entries are optional:
            FASTQ2              - second lane of reads.  If this does not exist, we treat
                                  this lane as single-end reads.  If it does, we treat it
                                  as paired-end reads.
            LIBRARY_NAME        - Library name.  Example: Pond-81791
            PLATFORM            - Sequencing platform.  Example: illumina
            PLATFORM_UNIT       - Platform unit.  Typically flowcell.lane.barcode.
                                  Example: C02B4ACXX110915.8.CAGCAAGG
            SEQUENCING_CENTER   - Sequencing center.  Example: BI
            RUN_DATE            - Date the run was produced in an ISO8601 format.
                                  Example: 2011-09-15T00:00:00-0400
        Fastq files may be gzip compressed.
    '''
    assert outBam.endswith('.bam') or outBam.endswith('.sam')

    # process each pair of fastqs
    first_bams = []
    for lane in pairs_list:
        assert 'FASTQ' in lane and 'QUALITY_FORMAT' in lane
        assert lane['QUALITY_FORMAT'] in ('Solexa','Illumina','Standard')
        rg_id = "%s.%d" % (urllib.quote_plus(sample), len(first_bams)+1)
        bam = util_files.mkstempfname(prefix='bundle-indiv_bam-%s-'%rg_id, suffix='.bam')
        cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s/FastqToSam.jar OUTPUT=%s READ_GROUP_NAME=%s SAMPLE_NAME=%s SORT_ORDER=unsorted" % (
            tempfile.tempdir, picardPath, bam, rg_id, sample)
        for k in ('QUALITY_FORMAT','FASTQ','FASTQ2','LIBRARY_NAME','PLATFORM','PLATFORM_UNIT','SEQUENCING_CENTER','RUN_DATE'):
            if k in lane:
                cmdline += " %s=%s" % (k,lane[k])
        log.info("fastq->bam: %s" % cmdline)
        assert not os.system(cmdline)
        first_bams.append(bam)

    # merge bams with picard
    merge_bams(first_bams, outBam, sort_order='unsorted', picardPath=picardPath)
    for bam in first_bams:
        os.unlink(bam)
    return outBam


def parser_fastqs_to_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inTable",
        version=version,
        description='''Take an input file describing a set of fastq files and associated
metadata.  Convert into a set of unaligned BAM files, one per sample, with metadata
loaded into the read group tags.  The input file should have a header line and the
following named fields: sample, fastq1, fastq2, bam, qual_format, platform, center.
Optionally also: library, run_date, platform_unit, and any number of other ignored
columns.  Reads will be read from fastq1 and fastq2 and written to bam.  If fastq2
is empty, we treat it as a single-end read lane, otherwise, it is paired-end.  Fastq
files may be gzip compressed.

Here are a description of the fields:
    qual_format         - Solexa, Illumina, or Standard
                            Solexa - pre-pipeline 1.3 scores (solexa scaling + 66)
                            Illumina - pipeline >=1.3 scores (phred scaling + 64)
                            Standard - phred scaling + 33
    library             - Library name.  Example: Pond-81791
    platform            - Sequencing platform.  Example: illumina
    platform_unit       - Platform unit.  Typically flowcell.lane.barcode.
                            Example: C02B4ACXX110915.8.CAGCAAGG
    center              - Sequencing center.  Example: BI
    run_date            - Date the run was produced in an ISO8601 format.
                            Example: 2011-09-15T00:00:00-0400
''')
    parser.add_option("--sample_name", dest="sample_name", type='string',
                      help="Only perform operations on the given sample name.",
                      default=None)
    parser.add_option("--sample_num", dest="sample_num", type='int',
                      help="Only perform operations on the ith sample.  1-based index.",
                      default=None)
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_fastqs_to_bam(args, options):
    inTable = args[0]
    assert not (options.sample_name and options.sample_num)
    assert options.sample_num==None or options.sample_num>0

    # parse input table
    samples = []
    pairs_list = []
    sample_to_bam = {}
    sample_idx = 0
    with open(inTable, 'rt') as inf:
        for row in util_files.FlatFileParser(inf):
            # basic checks
            for h in ('sample','fastq1','bam','qual_format','platform','center'):
                assert h in row and row[h]
            assert 'fastq2' in row
            assert row['qual_format'] in ('Solexa','Illumina','Standard')

            # if specified, restrict to a single sample
            sample_idx += 1
            if options.sample_num and options.sample_num != sample_idx:
                continue
            if options.sample_name and options.sample_name != row['sample']:
                continue

            # check for fastq files and fetch if missing
            if not os.access(row['fastq1'], os.R_OK):
                assert 'url1' in row
                log.info("fastq file %s does not exist, downloading from %s instead" % (row['fastq1'], row['url1']))
                cmdline = "wget %s -O %s -nc -nd --progress=dot:mega --ftp-password=dpark@broadinstitute.org" % (row['url1'], row['fastq1'])
                assert not os.system(cmdline)
            if row['fastq2'] and not os.access(row['fastq2'], os.R_OK):
                assert 'url2' in row
                log.info("fastq file %s does not exist, downloading from %s instead" % (row['fastq2'], row['url2']))
                cmdline = "wget %s -O %s -nc -nd --progress=dot:mega --ftp-password=dpark@broadinstitute.org" % (row['url2'], row['fastq2'])
                assert not os.system(cmdline)

            # check output file
            if row['sample'] in sample_to_bam:
                assert sample_to_bam[row['sample']] == row['bam']
            else:
                sample_to_bam[row['sample']] = row['bam']
                assert row['bam'].endswith('.bam')
                samples.append(row['sample'])

            # prep output row
            out = {'sample':row['sample'],
                'FASTQ':row['fastq1'], 'QUALITY_FORMAT':row['qual_format'],
                'PLATFORM':row['platform'], 'SEQUENCING_CENTER':row['center'],
                }
            if 'fastq2' in row and row['fastq2']:
                out['FASTQ2'] = row['fastq2']
            if 'library' in row and row['library']:
                out['LIBRARY_NAME'] = urllib.quote_plus(row['library'])
            if 'platform_unit' in row and row['platform_unit']:
                out['PLATFORM_UNIT'] = row['platform_unit']
            if 'run_date' in row and row['run_date']:
                out['RUN_DATE'] = row['RUN_DATE']
            pairs_list.append(out)

    assert pairs_list, "no samples found!"

    # perform bundling for each sample
    for s in samples:
        pairs = [x for x in pairs_list if x['sample']==s]
        log.info("bundling sample %s into a single bam from %d lanes to %s" % (s, len(pairs), sample_to_bam[s]))
        bundle_fastqs_to_unaligned_bam(pairs, sample_to_bam[s], s, picardPath=options.picardPath)
    return 0

def parser_merge_bams_from_file(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inTable outDir",
        version=version,
        description='''Take a two-column, headerless, tab-delimited input file describing
a set of bam files and sample names.  Convert into a set of aggregated BAM files, one per
sample, with metadata transferred over from the original read group tags, but with the
sample name re-written.  This provides a mechanism to aggregate bams that were provided
separately for some reason or another, but we didn't necessarily want sequencing to do
the aggregation work for us.  For example, we may want sequencing to treat culture-adapted
and hybrid-select runs of the same parasite separately, but we may want to merge the data
for certain analyses.  Output will go to a single bam file named by the sample name in
the outDir.
''')
    parser.add_option("--sample_name", dest="sample_name", type='string',
                      help="Only perform operations on the given sample name.",
                      default=None)
    parser.add_option("--sample_num", dest="sample_num", type='int',
                      help="Only perform operations on the ith sample.  1-based index.",
                      default=None)
    parser.add_option("--sort_order", dest="sort_order", type='string',
                      help="Sort order of input and output files.  Must be unsorted, queryname, or coordinate. [default: %default]",
                      default='coordinate')
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--samtoolsPath", dest="samtoolsPath", type='string',
                      help="Directory for samtools programs.  [default: %default]",
                      default=global_tool_paths['samtools'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_merge_bams_from_file(args, options):
    inTable, outDir = args
    assert not (options.sample_name and options.sample_num)

    # parse input table
    samples = []
    inFiles = {}
    sample_idx = 0
    with open(inTable, 'rt') as inf:
        for line in inf:
            sample, bam = line.rstrip('\r\n').split('\t')
            if sample not in inFiles:
                samples.append(sample)
                inFiles[sample] = []
            inFiles[sample].append(bam)

    # if sample_num is set to zero, just print out our list of samples instead
    if options.sample_num==0:
        for i in range(len(samples)):
            print("%d\t%s\t%d" % (i+1, samples[i], len(inFiles[samples[i]])))
        return 0

    # if specified, restrict to a single sample
    if options.sample_num:
        assert not options.sample_name
        assert 1 <= options.sample_num <= len(samples)
        log.info("restricted to sample %s (%d of %d)" % (samples[options.sample_num-1], options.sample_num, len(samples)))
        samples = [samples[options.sample_num-1]]
    if options.sample_name:
        assert not options.sample_num
        assert options.sample_name in inFiles
        log.info("restricted to sample %s (%d of %d)" % (options.sample_name, samples.index(options.sample_name), len(samples)))
        samples = [options.sample_name]

    # perform bundling for each sample
    tmpBam = util_files.mkstempfname(prefix='merge_bams-', suffix='.bam')
    for s in samples:
        outBam = "%s/%s.bam" % (outDir, s)
        log.info("bundling sample %s into a single bam from %d files, sort order is %s" % (s, len(inFiles[s]), options.sort_order))
        merge_bams(inFiles[s], tmpBam, sort_order=options.sort_order, picardPath=options.picardPath)
        rename_bam_sample(tmpBam, s, outBam, samtoolsPath=options.samtoolsPath)
        os.unlink(tmpBam)
        if options.sort_order=='coordinate':
            index_bam(outBam)
    return 0

def parser_realign(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam inRef outBam",
        version=version,
        description='''Re-align reads in a BAM file to a new reference.
Sort, index, mark duplicates and locally realign around indels.''')
    parser.add_option("--sample", dest="sample", type='string',
                      help="If provided, rename the sample name.",
                      default=None)
    parser.add_option("--paired", dest="paired", type='int',
                      help="Are reads paired (paired=2) or single (paired=1)?  Default is auto-detect.",
                      default=None)
    parser.add_option("--knownSites", dest="knownSites", type='string',
                      help="If a database of known variation is available, it will be used to provide base quality recalibration",
                      default=None)
    parser.add_option("--noRevert",
                      action="store_true", dest="noRevert",
                      default=False,
                      help="""By default, we will run Picard's RevertSam on the
input file before doing anything else to strip all of the recalibration and other
information.  Specify --noRevert in order to skip this step (for example, if the
input bam is already unaligned).""")
    parser.add_option("--samtoolsPath", dest="samtoolsPath", type='string',
                      help="Directory for samtools programs.  [default: %default]",
                      default=global_tool_paths['samtools'])
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--bwaPath", dest="bwaPath", type='string',
                      help="Directory for BWA.  [default: %default]",
                      default=global_tool_paths['bwa'])
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK programs.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_realign(args, options):
    inBam, inRef, outBam = args
    knownSites = []
    if options.knownSites:
        knownSites.append(options.knownSites)

    log.info("realigning reads to new reference with BWA")
    realignBam = util_files.mkstempfname(prefix='realign-bwa-', suffix='.bam')
    realign_bam(inBam, inRef, realignBam,
        newSampleName=options.sample, paired=options.paired, revertBam=not options.noRevert,
        bwaPath=options.bwaPath, samtoolsPath=options.samtoolsPath, picardPath=options.picardPath)
    post_process_bam(realignBam, outBam, inRef, knownSites=knownSites, sort=False,
        picardPath=options.picardPath, gatkPath=options.gatkPath)
    os.unlink(realignBam)
    os.unlink(realignBam[:-4]+'.bai')

    log.info("done")
    return 0

def parser_realign_from_file(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inTable inRef outBamDir index",
        version=version,
        description='''Re-align reads in a BAM file to a new reference.
Sort, index, mark duplicates and locally realign around indels.
We will read the index-th line from the inTable file, take the first column
as the sample name and the second column as the input bam file.  The output
bam file will be outBamDir/samplename.bam.  index is 1-based and there is no
header line in inTable.''')
    parser.add_option("--knownSites", dest="knownSites", type='string',
                      help="If a database of known variation is available, it will be used to provide base quality recalibration",
                      default=None)
    parser.add_option("--noRevert",
                      action="store_true", dest="noRevert",
                      default=False,
                      help="""By default, we will run Picard's RevertSam on the
input file before doing anything else to strip all of the recalibration and other
information.  Specify --noRevert in order to skip this step (for example, if the
input bam is already unaligned).""")
    parser.add_option("--samtoolsPath", dest="samtoolsPath", type='string',
                      help="Directory for samtools programs.  [default: %default]",
                      default=global_tool_paths['samtools'])
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--bwaPath", dest="bwaPath", type='string',
                      help="Directory for BWA.  [default: %default]",
                      default=global_tool_paths['bwa'])
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK programs.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_realign_from_file(args, options):
    inTable, inRef, outBamDir, index = args
    index = int(index)

    knownSites = []
    if options.knownSites:
        knownSites.append(options.knownSites)
        log.info("reading known sites from %s" % options.knownSites)

    log.info("reading %d-th line from %s" % (index, inTable))
    sample, inBam = (None, None)
    with open(inTable, 'rt') as inf:
        i=1
        assert 1<=index
        for line in inf:
            if i==index:
                sample, inBam = line.rstrip('\r\n').split('\t')
                break
            i += 1
    assert sample and inBam
    outBam = "%s/%s.bam" % (outBamDir, sample)
    log.info("sample %s, reading input from %s, writing output to %s" % (sample, inBam, outBam))

    realignBam = util_files.mkstempfname(prefix='realign-bwa-', suffix='.bam')
    realign_bam(inBam, inRef, realignBam, newSampleName=sample, revertBam=not options.noRevert,
        bwaPath=options.bwaPath, samtoolsPath=options.samtoolsPath, picardPath=options.picardPath)
    post_process_bam(realignBam, outBam, inRef, knownSites=knownSites, sort=False,
        picardPath=options.picardPath, gatkPath=options.gatkPath)
    os.unlink(realignBam)
    os.unlink(realignBam[:-4]+'.bai')

    log.info("done")
    return 0


def parser_subset_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam outBam",
        version=version,
        description='''Filter BAM file.''')
    parser.add_option("--region", dest="region", type='string',
                      help="Keep only reads aligned to the specified genomic region.  Example: chr2:100000-200000",
                      default='')
    parser.add_option("--read_group", dest="read_group", type='string',
                      help="Keep only reads from RG's that have the following field set.  Examples: CN=BI, DT=2012-02-07T00:00:00-0500, or SM=SenT179.08.h",
                      default=None)
    parser.add_option("--mapq", dest="mapq", type='int',
                      help="Minimum MAPQ.  [default: %default]",
                      default=0)
    parser.add_option("--flags_remove", dest="flags_remove", type='int',
                      help="Remove all reads that have any of the following flags set.  [default: %default]",
                      default=4)
    parser.add_option("--flags_keep", dest="flags_keep", type='int',
                      help="Keep only reads that have all of the following flags set.  [default: %default]",
                      default=0)
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--samtoolsPath", dest="samtoolsPath", type='string',
                      help="Directory for samtools programs.  [default: %default]",
                      default=global_tool_paths['samtools'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_subset_bam(args, options):
    inBam, outBam = args
    log.info("filtering %s to %s" % (inBam, outBam))

    cmdline = "%s/samtools view -bh" % options.samtoolsPath
    if options.read_group:
        key, val = options.read_group.split('=')
        rg_file = util_files.mkstempfname(prefix='read_groups-', suffix='.sam')
        cmdline2 = "%s/samtools view -H %s | grep ^@RG > %s" % (options.samtoolsPath, inBam, rg_file)
        assert not os.system(cmdline2)
        rg_id = None
        with open(rg_file, 'rt') as inf:
            for row in inf:
                for x in row.rstrip('\n\r').split('\t'):
                    if x.startswith('ID:'):
                        id = x[3:]
                    if x.startswith(key+':') and x[len(key)+1:]==val:
                        assert not rg_id, "multiple read groups matched %s=%s" % (key,val)
                        rg_id = id
        os.unlink(rg_file)
        assert rg_id, "failed to find any read groups that matched %s=%s" % (key,val)
        log.info("restricting to read group %s, which matches %s=%s" % (rg_id, key, val))
        cmdline += " -r %s" % rg_id
    if options.mapq:
        cmdline += " -q %d" % options.mapq
    if options.flags_remove:
        cmdline += " -F %d" % options.flags_remove
    if options.flags_keep:
        cmdline += " -f %d" % options.flags_keep
    cmdline += " -o %s %s %s" % (outBam, inBam, options.region)
    log.info("subset_bam: %s" % cmdline)
    assert not os.system(cmdline)

    log.info("indexing bam")
    index_bam(outBam, picardPath=options.picardPath)
    log.info("done")
    return 0

def parser_rename_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam sampleName outBam",
        version=version,
        description='''Replace sample name in read groups in BAM header.''')
    parser.add_option("--samtoolsPath", dest="samtoolsPath", type='string',
                      help="Directory for samtools programs.  [default: %default]",
                      default=global_tool_paths['samtools'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_rename_bam(args, options):
    inBam, sampleName, outBam = args
    rename_bam_sample(inBam, sampleName, outBam, samtoolsPath=options.samtoolsPath)
    return 0

def parser_sort_index_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam outBam",
        version=version,
        description='''Sort BAM file with samtools.  If inBam=outBam, replace
input file with output file.  Otherwise, write a new file.''')
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_sort_index_bam(args, options):
    inBam, outBam = args
    sort_bam(inBam, outBam, picardPath=options.picardPath)
    index_bam(outBam, picardPath=options.picardPath)
    return 0

def parser_index_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] bamFile",
        version=version,
        description='''Index BAM file with picard.''')
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    return parser
def main_index_bam(args, options):
    assert len(args)==1
    bamfile = args[0]
    index_bam(bamfile, picardPath=options.picardPath)
    return 0


def parser_revert_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam outBam",
        version=version,
        description='''Revert BAM file to unaligned state.''')
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    return parser
def main_revert_bam(args, options):
    inBam, outBam = args
    revert_bam(inBam, outBam, picardPath=options.picardPath)
    return 0


def parser_post_process_bam(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inBam outBam",
        version=version,
        description='''Mark read duplicates in BAM file with Picard.
If inBam=outBam, replace input file with output file.  Otherwise, write a new file.''')
    parser.add_option("--knownSites", dest="knownSites", type='string',
                      help="If a database of known variation is available, it will be used to provide base quality recalibration",
                      default=None)
    parser.add_option("--picardPath", dest="picardPath", type='string',
                      help="Directory for Picard.  [default: %default]",
                      default=global_tool_paths['picard'])
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK programs.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_post_process_bam(args, options):
    inBam, outBam = args
    knownSites = options.knownSites and [options.knownSites] or []
    post_process_bam(inBam, outBam, knownSites=knownSites, picardPath=options.picardPath, gatkPath=options.gatkPath)
    return 0


def parser_gatk_ug(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] bamfile1 [bamfile2 bamfile3 ..] ref outfile",
        version=version,
        description='''Take input BAM files of alignments and a reference genome
FASTA file.  Use GATK's Unified Genotyper to call consensus sequence.  BAM file must
be indexed by samtools/Picard (.bai) and FASTA file must be indexed by
Picard/GATK (.dict).''')
    parser.add_option("--ploidy", dest="ploidy", type='int',
                      help="# of chromosome copies per individual  [default: %default]",
                      default=1)
    parser.add_option("--interval_idx", dest="interval_idx", type='int',
                      help="""Break the genome into interval_total pieces and compute
only on the interval_idx'th piece.  interval_idx is one-based.  Default is to compute
on the whole genome (unbroken).""",
                      default=None)
    parser.add_option("--interval_total", dest="interval_total", type='int',
                      help="""Break the genome into interval_total pieces and compute
only on the interval_idx'th piece.  interval_idx is one-based.  Default is to compute
on the whole genome (unbroken).""",
                      default=None)
    parser.add_option("--num_threads", dest="num_threads", type='int',
                      help="# of threads to execute in [default unspecified]",
                      default=None)
    parser.add_option("--glm", dest="glm", type='choice',
                      help="Calculates SNP, INDEL, or BOTH.  [default: %default]",
                      choices=['SNP','INDEL','BOTH'],
                      default='BOTH')
    parser.add_option("--minQUAL", dest="minQUAL", type='float',
                      help="Min QUAL value.  [default: %default]",
                      default=60.0)
    parser.add_option("--minGQ", dest="minGQ", type='int',
                      help="Min GQ value.  [default: %default]",
                      default=30)
    parser.add_option("--keepHets",
                      action="store_true", dest="keepHets",
                      default=False,
                      help="Default is to filter out all het genotypes.  Specify --keepHets to disable.")
    parser.add_option("--hardFilter",
                      action="store_true", dest="hardFilter",
                      default=False,
                      help="Default is to soft-filter by marking low quality and het bases with the FILTER and FT flags.  Specify this to physically remove these positions and genotypes from the output VCF altogether.")
    parser.add_option("--oneSample", dest="oneSample", type='int',
                      help="""Default is to run all samples in all bam files into a single
output VCF file (outfile).  If this option is specified (a number between 1 and the number
of input bam files), only run one of the bam files.  Outfile is interpreted as a directory
and the output will go to a .vcf.gz file in the output directory with the same base
file name as the input bam file.""",
                      default=None)
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK programs.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser

def main_gatk_ug(args, options):
    assert len(args)>=3
    assert options.interval_idx==None and options.interval_total==None or options.interval_idx<=options.interval_total and options.interval_total>1
    ref, outfile = args[-2:]
    bamfiles = args[:-2]
    if options.oneSample:
        assert os.path.isdir(outfile) and 1 <= options.oneSample <= len(bamfiles)
        bamfiles = bamfiles[options.oneSample-1:options.oneSample]
        assert bamfiles[0].endswith('.bam')
        outfile = outfile.endswith('/') and outfile[:-1] or outfile
        outfile = "%s/%s.vcf.gz" % (outfile, bamfiles[0][:-4].split('/')[-1])
        log.info("--oneSample=%d specified: reading from only %s and writing to %s" % (options.oneSample, bamfiles[0], outfile))
    assert outfile.endswith('.vcf.gz')
    tmpVcf1 = util_files.mkstempfname(prefix='gatk_ug-', suffix='.vcf')
    tmpVcf2 = util_files.mkstempfname(prefix='gatk_ug-', suffix='.vcf')

    log.info("starting gatk_ug on %d bam files" % len(bamfiles))
    gatk_ug(bamfiles, ref, tmpVcf1,
        glm=options.glm, ploidy=options.ploidy, num_threads=options.num_threads,
        interval_idx=options.interval_idx, interval_total=options.interval_total,
        gatkPath=options.gatkPath, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)

    log.info("starting VCF soft filter")
    #vcf_filter(tmpVcf1, ref, tmpVcf2, gatkPath=options.gatkPath, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    vcf_filter_soft(tmpVcf1, tmpVcf2, minQUAL=options.minQUAL, minGQ=options.minGQ, removeHets=not options.keepHets, resetFilters=True)
    os.unlink(tmpVcf1)

    if options.hardFilter:
        log.info("starting VCF hard filter")
        vcf_filter_hard(tmpVcf2, tmpVcf1)
        os.unlink(tmpVcf2)
    else:
        tmpVcf1 = tmpVcf2

    log.info("compressing and indexing")
    util_vcf.vcf_bgzip_index(tmpVcf1, outfile, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf1)
    return 0


def parser_vcf_merge(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf1 inVcf2 [inVcf3 inVcf4 ...] refFasta outFile",
        version=version,
        description='''Take input VCF files and merge together to a single VCF file
using GATK CombineVariants.''')
    parser.add_option("--concat",
                      action="store_true", dest="concat",
                      default=False,
                      help="""If true, we pass --assumeIdenticalSamples to
CombineVariants.  This assumes that all input files describe the same set of samples
across disjoint sets of markers, and that the files simply need to be concatenated
together in the right order.  If false (default), we assume input files describe
the same set of markers across disjoint sets of samples.""")
    parser.add_option("--interval_idx", dest="interval_idx", type='int',
                      help="""Break the genome into interval_total pieces and compute
only on the interval_idx'th piece.  interval_idx is one-based.  Default is to compute
on the whole genome (unbroken).""",
                      default=None)
    parser.add_option("--interval_total", dest="interval_total", type='int',
                      help="""Break the genome into interval_total pieces and compute
only on the interval_idx'th piece.  interval_idx is one-based.  Default is to compute
on the whole genome (unbroken).""",
                      default=None)
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools programs.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix and bgzip.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser

def main_vcf_merge(args, options):
    assert len(args)>=4
    refFasta, outFile = args[-2:]
    inFiles = args[:-2]
    for f in inFiles:
        assert os.access(f, os.R_OK), "input file '%s' is not readable" % f
        assert f.endswith('.vcf.gz')
    assert not os.access(outFile, os.F_OK) or os.access(outFile, os.W_OK), "outFile '%s' is not a writeable destination" % outFile
    tmpVcf = util_files.mkstempfname(prefix='gatk_cv-', suffix='.vcf')

    # vcf-merge
    log.info("running CombineVariants on %d input files" % len(inFiles))
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T CombineVariants -R %s %s -o %s --genotypemergeoption REQUIRE_UNIQUE" % (
        tempfile.tempdir, options.gatkPath+'/GenomeAnalysisTK.jar',
        refFasta,
        ' '.join(['--variant '+i for i in inFiles]),
        tmpVcf)
    if options.concat:
        cmdline += " --assumeIdenticalSamples"
    if options.interval_idx!=None and options.interval_total:
        for c,start,stop in util_vcf.make_intervals(options.interval_idx, options.interval_total, refFasta):
            cmdline += ' -L %s:%d-%d' % (c,start,stop)
    assert not os.system(cmdline)

    # bgzip to outFile and index with tabix and vcftools
    util_vcf.vcf_bgzip_index(tmpVcf, outFile, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf)
    log.info("done")
    return 0


def parser_vcf_rename(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf conversionString outVcf",
        version=version,
        description='''Rename samples in a VCF file by changing its header row.
Do this based on a conversion string which is a comma-separated list of
equal-sign-separated pairs which contain the old sample name and the new sample name.
Example: oldsample1=newsample1,oldsample2=newsample2.''')
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix and bgzip.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser

def main_vcf_rename(args, options):
    inVcf, conversionString, outVcf = args
    tmpFile1 = util_files.mkstempfname(suffix='.vcf', prefix='header_old-')
    tmpFile2 = util_files.mkstempfname(suffix='.vcf', prefix='header_new-')

    sampleMap = dict([x.split('=') for x in conversionString.split(',')])
    for k,v in sampleMap.items():
        log.info("mapping sample name %s to %s" % (k,v))

    cmdline = "%s/tabix -h %s null:1-1 > %s" % (options.tabixPath, inVcf, tmpFile1)
    assert not os.system(cmdline)

    with open(tmpFile1, 'rt') as inf:
        with open(tmpFile2, 'wt') as outf:
            for line in inf:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        row = line.rstrip('\r\n').split('\t')
                        head = row[:9]
                        samples = row[9:]
                        samples = map(sampleMap.get, samples)
                        line = '\t'.join(head+samples)+'\n'
                    outf.write(line)

    log.info("rewriting header with tabix to %s" % outVcf)
    cmdline = "%s/tabix -r %s %s > %s" % (options.tabixPath, tmpFile2, inVcf, outVcf)
    assert not os.system(cmdline)

    log.info("indexing with tabix")
    cmdline = "%s/tabix %s -f -p vcf" % (options.tabixPath, outVcf)
    assert not os.system(cmdline)

    os.unlink(tmpFile1)
    os.unlink(tmpFile2)
    log.info("done")
    return 0


def parser_ena_maketable(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] samples_file erp_file fastq_dir bam_dir outFile",
        version=version,
        description='''Take samples listed in samples_file, cross reference to the
erp_file, find fastqs in fastq_dir, and produce a table suitable for fastqs_to_bam.''')
    parser.add_option("--qual_format", dest="qual_format", type='string',
                      help="Fastq quality score format (Solexa, Illumina, Standard) [default %default]",
                      default='Standard')
    parser.add_option("--center", dest="center", type='string',
                      help="Sequencing center name.  [default: %default]",
                      default='BI')

    return parser
def main_ena_maketable(args, options):
    samples_file, erp_file, fastq_dir, bam_dir, outFile = args
    assert options.qual_format in ('Solexa','Illumina','Standard')

    # sample info from samples_file and erp_file
    samples = []
    with open(samples_file, 'rt') as inf:
        for line in inf:
            line = line.rstrip('\r\n').split('\t')
            assert len(line) == 2
            samples.append((line[0], line[1].split(',')))
    with open(erp_file, 'rt') as inf:
        erp_table = [row for row in util_files.FlatFileParser(inf)]

    out = []
    for sample, ers_ids in samples:
        outBam = "%s/%s.bam" % (bam_dir, sample)

        # subset main table and verify that all ers_ids are present
        erp_rows = [row for row in erp_table if row['Sample'] in ers_ids]
        for ers in ers_ids:
            assert len([1 for row in erp_rows if row['Sample']==ers]), "sample %s, accession %s missing from %s" % (sample, ers, erp_file)

        # group into runs (lanes) and build info for each
        runs = unique([r['Run'] for r in erp_rows])
        for run in runs:
            erp_run_rows = [row for row in erp_rows if row['Run']==run]
            assert len(erp_run_rows) in (1,2)
            fastqs = ["%s/%s" % (fastq_dir, row['Ftp'].split('/')[-1]) for row in erp_run_rows]
            run = {
                'sample':sample, 'bam':outBam,
                'fastq1':fastqs[0], 'url1':erp_run_rows[0]['Ftp'],
                'center':options.center, 'qual_format':options.qual_format,
                'library':erp_run_rows[0]['Library Name'],
                'platform':erp_run_rows[0]['Instrument Platform'],
                'platform_unit':erp_run_rows[0]['Run'],
            }
            if len(erp_run_rows)==2:
                run['fastq2'] = fastqs[1]
                run['url2'] = erp_run_rows[1]['Ftp']
                for k in ('Library Name', 'Instrument Platform'):
                    assert erp_run_rows[0][k] == erp_run_rows[1][k]
            else:
                run['fastq2'] = ''
                run['url2'] = ''
            out.append(run)

    with open(outFile, 'wt') as outf:
        header = ('sample','center','library','platform','platform_unit','qual_format','bam','fastq1','fastq2','url1','url2')
        outf.write('\t'.join(header)+'\n')
        for run in out:
            outf.write('\t'.join([run[h] for h in header])+'\n')
    log.info("done")
    return 0


def vcf_filter_soft(inVcf, outVcf,
    minQUAL=60, minGQ=30, removeHets=True, resetFilters=True):
    ''' Filter a VCF file based on our own criteria.  Unfortunately, GATK's
        VariantFiltration seems unable to perform these filters properly
        even though it ought to be able to.
        This handles VCF files with any ploidy, unlike vcftools.  Handles any phasing.
    '''
    assert (inVcf.endswith('.vcf') or inVcf.endswith('.vcf.gz')) and outVcf.endswith('.vcf')
    assert 0 <= minQUAL <= 255
    assert 0 <= minGQ <= 99
    log.info("soft filtering VCF file from %s to %s" % (inVcf, outVcf))
    with util_files.open_or_gzopen(inVcf, 'rt') as inf:
        with open(outVcf, 'wt') as outf:
            nrow = 0
            for row in inf:
                if not row.startswith('#'):
                    row = vcfRow_soft_filter(row, minQUAL, minGQ, removeHets, resetFilters)
                    nrow += 1
                    if not (nrow % 1000000):
                        log.info("processed %dM sites" % (nrow/1000000))
                else:
                    if row.startswith('##fileformat'):
                        filt_head = [
                            ('LowQUAL', 'QUAL<%s'%minQUAL),
                            ('LowGQ', 'GQ<%s'%minGQ),
                            ('MissingGT', 'GT is missing'),
                            ('NoPassingGenotypes', 'Non-passing FT values for all samples'),
                            ]
                        if removeHets:
                            filt_head.append(('het', 'isHet==1'))
                        for id,desc in filt_head:
                            row += '##FILTER=<ID=%s,Description="%s">\n' % (id,desc)
                        row += '##VcfFilter="minQUAL=%s minGQ=%s removeHets=%s resetFilters=%s"\n' % (minQUAL, minGQ, removeHets, resetFilters)
                    elif row.startswith('##FILTER'):
                        if resetFilters or row[13:row.find(',')] in ('LowQUAL','LowGQ','MissingGT','NoPassingGenotypes','het'):
                            continue
                outf.write(row)
    return outVcf

def vcf_filter_hard(inVcf, outVcf):
    ''' Hard-filter a VCF file, removing genotype calls with FT flags and rows
        with FILTER flags and rows where all genotypes are either missing or FT.
        If dropped genotype calls result in less alleles, adjust the ALT column
        accordingly.  Handles genotypes of any ploidy or phasing.
        Surprising that GATK and vcftools all lack a way to do this.
    '''
    assert (inVcf.endswith('.vcf') or inVcf.endswith('.vcf.gz')) and outVcf.endswith('.vcf')
    log.info("hard filtering VCF file from %s to %s" % (inVcf, outVcf))
    with util_files.open_or_gzopen(inVcf, 'rt') as inf:
        with open(outVcf, 'wt') as outf:
            n_pos_in = 0
            n_pos_out = 0
            n_genos_in = 0
            n_genos_out = 0
            counters = {}
            for row in inf:
                if not row.startswith('#'):
                    row, row_stats = vcfRow_hard_filter(row)
                    n_pos_in += 1
                    n_pos_out += (row and 1 or 0)
                    n_genos_in += row_stats['genos_in']
                    n_genos_out += row_stats['genos_out']
                    counters.setdefault(row_stats['drop_reason'], 0)
                    counters[row_stats['drop_reason']] += 1
                    if not (n_pos_in % 1000000):
                        log.info("processed %dM sites" % (n_pos_in/1000000))
                if row:
                    outf.write(row)
    log.info(".. # genotypes in:%d\t# genotypes out: %d" % (n_genos_in, n_genos_out))
    log.info(".. # positions in:%d\t# positions out: %d" % (n_pos_in, n_pos_out))
    for k,v in counters.items():
        log.info(".... # positions with %s: %d" % (k,v))
    return outVcf

def vcfRow_soft_filter(row, minQUAL, minGQ, removeHets, resetFilters):
    row = row.rstrip('\r\n').split('\t')
    genos = row[9:]
    row = row[:9]
    c, p, id, ref, alt, qual, filter, info, format = row
    n_genos_out = 0

    # set up FILTER as list
    if resetFilters or filter in ('PASS','','.'):
        filter = []
    else:
        filter = filter.split(';')

    # filter if QUAL too low
    if qual=='.' or minQUAL > float(qual):
        filter.append('LowQUAL')

    # get GT,GQ,FT fields on genotypes
    format = format.split(':')
    assert format[0]=='GT'
    if 'FT' not in format:
        format.append('FT')
        row[8] = ':'.join(format)
    ft_idx = format.index('FT')
    if 'GQ' in format:
        gq_idx = format.index('GQ')
    else:
        gq_idx = -1

    # loop through each sample
    for i in range(len(genos)):
        g = genos[i].split(':')
        assert len(g)<=len(format), "bad VCF file at %s:%s - FORMAT (%d): %s - geno for sample %d (%d): %s" % (c,p,len(format),':'.join(format), i+1, len(g), ':'.join(g))
        if len(g)<len(format):
            # apparently, this is allowed in the VCF spec.. just pad it out with missing dots
            g.extend(list('.' * (len(format) - len(g))))

        # set up FT as list:
        if resetFilters or g[ft_idx] in ('PASS','.'):
            ft = []
        else:
            ft = g[ft_idx].split(';')

        # parse GQ
        if gq_idx>=0 and g[gq_idx]!='.' and minGQ>int(g[gq_idx]):
            ft.append('LowGQ')

        # parse GT
        alleles = set(g[0].replace('|','/').split('/'))
        if '.' in alleles:
            ft.append('MissingGT')
        elif removeHets and len(alleles)>1:
            ft.append('het')

        # finalize changes to sample
        if ft:
            g[ft_idx] = ';'.join(ft)
        else:
            g[ft_idx] = 'PASS'
            n_genos_out += 1
        genos[i] = ':'.join(g)

    # filter if all genotypes were filtered out
    if n_genos_out==0:
        filter.append('NoPassingGenotypes')

    # finalize changes to marker and return new row
    row[6] = filter and ';'.join(filter) or 'PASS'
    return '\t'.join(row + genos)+'\n'

def vcfRow_hard_filter(row):
    row = row.rstrip('\r\n').split('\t')
    genos = row[9:]
    row = row[:9]
    c, p, id, ref, alt, qual, filter, info, format = row
    row_stats = {'genos_in':0,'genos_out':0,'drop_reason':'PASS'}

    # drop if FILTER
    if filter and filter not in ('PASS','.'):
        row_stats['genos_in'] = len(genos)
        row_stats['drop_reason'] = filter
        return (None, row_stats)

    # check for FT fields on genotypes
    format = format.split(':')
    assert format[0]=='GT'
    if 'FT' in format:
        ft_idx = format.index('FT')
    else:
        ft_idx = -1

    # count ALT alleles seen
    alt_alleles = set()

    # loop through each sample
    for i in range(len(genos)):
        g = genos[i].split(':')
        assert len(g)<=len(format), "bad VCF file at %s:%s - FORMAT (%d): %s - geno for sample %d (%d): %s" % (c,p,len(format),':'.join(format), i+1, len(g), ':'.join(g))
        if len(g)!=len(format):
            # apparently, this is allowed in the VCF spec.. just pad it out with missing dots
            g.extend(list('.' * (len(format) - len(g))))
        if '.' not in g[0]:
            # genotype is non-missing
            row_stats['genos_in'] += 1
            alleles = [int(x) for x in g[0].replace('|','/').split('/')]
            if ft_idx<0 or g[ft_idx] in ('PASS','','.'):
                # genotype should be kept
                row_stats['genos_out'] += 1
                for a in alleles:
                    if a:
                        alt_alleles.add(a)
            else:
                # genotype should be converted to missing (and unphased)
                ploidy = len(alleles)
                g[0] = '/'.join('.' * ploidy)
                genos[i] = ':'.join(g)

    # if we see a reduction in ALT alleles, adjust everything accordingly
    alt = alt!='.' and alt.split(',') or []
    if len(alt) > len(alt_alleles):
        # we lost ALT alleles in the filtering
        if not alt_alleles:
            # simple case: we lost all ALT alleles and became ref-only
            alt = '.'
        else:
            # we lost some ALT alleles
            alt_alleles = sorted(alt_alleles)
            if alt_alleles==range(1,1+len(alt_alleles)):
                # we lost only the last few ALT alleles and require no remapping
                alt = alt[:len(alt_alleles)]
            else:
                # we lost some non-terminal ALT alleles and require remapping
                alt = [alt[i-1] for i in alt_alleles]
                a_map = dict([(alt_alleles[i],i+1) for i in range(len(alt_alleles))])
                a_map[0] = 0
                for i in range(len(genos)):
                    g = genos[i].split(':')

                    # remap GT values
                    try:
                        alleles = [x!='.' and str(a_map[int(x)]) or '.'
                            for x in g[0].replace('|','/').split('/')]
                        g[0] = '/'.join(alleles)
                    except ValueError as e:
                        log.error("error with GT at sample %d at %s:%s - %s" % (i+1,c,p, g[0]))
                        raise

                    # remap PL values (TODO)
                    # remap GL values (TODO)
                    # remap GP values (TODO)

                    genos[i] = ':'.join(g)
            alt = ','.join(alt)
        #log.debug("reduction in number of ALT alleles at %s:%s from %s to %s" % (c,p,row[4],alt))
        row[4] = alt

    # drop row if all genotypes were filtered out
    if row_stats['genos_out']==0:
        row_stats['drop_reason'] = 'NoPassingGenotypes'
        row = None
    else:
        row = row + genos
        row = '\t'.join(row)+'\n'
    return (row, row_stats)

def vcf_filter_invariants(inVcf, outVcf):
    ''' Remove all invariant markers from VCF file, leaving only variants.
        Only looks at the ALT column of the VCF, does not actually check to see
        if all ALT genotypes were actually referred to in the genotype data.
    '''
    assert (inVcf.endswith('.vcf') or inVcf.endswith('.vcf.gz')) and outVcf.endswith('.vcf')
    log.info("filtering VCF file for variants only from %s to %s" % (inVcf, outVcf))
    with util_files.open_or_gzopen(inVcf, 'rt') as inf:
        with open(outVcf, 'wt') as outf:
            n_pos_in = 0
            n_pos_out = 0
            for row in inf:
                if not row.startswith('#'):
                    if row.split('\t',5)[4]=='.':
                        row = None
                    else:
                        n_pos_out += 1
                    n_pos_in += 1
                    if not (n_pos_in % 1000000):
                        log.info("processed %d M sites" % (n_pos_in/1000000))
                if row:
                    outf.write(row)
    log.info(".. # positions in:%d\t# positions out: %d" % (n_pos_in, n_pos_out))
    return outVcf


def parser_post_process_vcf(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf genome filteredVcf variantVcf annotatedVcf",
        version=version,
        description='''Take VCF and produce 1. a soft-filtered VCF, 2. a hard-filtered, variant-only VCF, and 3. a snpEff annotated VCF''')
    parser.add_option("--minQUAL", dest="minQUAL", type='float',
                      help="Min QUAL value.  [default: %default]",
                      default=60.0)
    parser.add_option("--minGQ", dest="minGQ", type='int',
                      help="Min GQ value.  [default: %default]",
                      default=30)
    parser.add_option("--keepHets",
                      action="store_true", dest="keepHets",
                      default=False,
                      help="Default is to filter out all het genotypes.  Specify --keepHets to disable.")
    parser.add_option("--keepOldFilters",
                      action="store_true", dest="keepOldFilters",
                      default=False,
                      help="Default is to ignore old filters and start fresh.  Specify --keepOldFilters to honor old FILTER and FT values.")
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--snpEffPath", dest="snpEffPath", type='string',
                      help="Directory for snpEff.  [default: %default]",
                      default=global_tool_paths['snpEff'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_post_process_vcf(args, options):
    inVcf, genome, filteredVcf, variantVcf, annotatedVcf = args
    tmpVcf1 = util_files.mkstempfname(prefix='vcf_filter-', suffix='.vcf')
    tmpVcf2 = util_files.mkstempfname(prefix='vcf_filter-', suffix='.vcf')

    vcf_filter_soft(inVcf, tmpVcf1,
        minQUAL=options.minQUAL, minGQ=options.minGQ,
        removeHets=not options.keepHets, resetFilters=not options.keepOldFilters)
    util_vcf.vcf_bgzip_index(tmpVcf1, filteredVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)

    vcf_filter_hard(filteredVcf, tmpVcf1)
    vcf_filter_invariants(tmpVcf1, tmpVcf2)
    os.unlink(tmpVcf1)
    util_vcf.vcf_bgzip_index(tmpVcf2, variantVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf2)

    vcf_snpEff(variantVcf, genome, annotatedVcf,
        snpEffPath=options.snpEffPath, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    return 0

def parser_vcf_bgzip_index(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf outVcf",
        version=version,
        description='''Take uncompressed VCF, produce bgzip-compressed vcf.gz, and index with tabix and vcftools.''')
    parser.add_option("--deleteInput",
                      action="store_true", dest="deleteInput",
                      default=False,
                      help="Delete input file when done.")
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_vcf_bgzip_index(args, options):
    inVcf, outVcf = args
    util_vcf.vcf_bgzip_index(inVcf, outVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    if options.deleteInput:
        os.unlink(inVcf)
    return 0

def parser_vcf_hard_filter(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf outVcf",
        version=version,
        description='''Take a soft-filtered VCF file and strip away all filtered
positions and genotypes.  Optionally reduce to just variant sites.''')
    parser.add_option("--removeInvariants",
                      action="store_true", dest="removeInvariants",
                      default=False,
                      help="Remove all invariant sites.")
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_vcf_hard_filter(args, options):
    inVcf, outVcf = args
    tmpVcf1 = util_files.mkstempfname(prefix='vcf_filter-', suffix='.vcf')
    tmpVcf2 = util_files.mkstempfname(prefix='vcf_filter-', suffix='.vcf')

    vcf_filter_hard(inVcf, tmpVcf1)
    if options.removeInvariants:
        vcf_filter_invariants(tmpVcf1, tmpVcf2)
        os.unlink(tmpVcf1)
    else:
        os.unlink(tmpVcf2)
        tmpVcf2 = tmpVcf1

    util_vcf.vcf_bgzip_index(tmpVcf2, outVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf2)
    return 0


def parser_snpEff(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf genome outVcf",
        version=version,
        description='''Run snpEff.''')
    parser.add_option("--asText",
                      action="store_true", dest="asText",
                      default=False,
                      help="If set, output will be in snpEff's plaintext output format.  Default is VCF format.")
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=None)
    parser.add_option("--snpEffPath", dest="snpEffPath", type='string',
                      help="Directory for snpEff.  [default: %default]",
                      default=global_tool_paths['snpEff'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser

def main_snpEff(args, options):
    inVcf, genome, outVcf = args
    vcf_snpEff(inVcf, genome, outVcf,
        vcf_out = not options.asText,
        snpEffPath=options.snpEffPath, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    return 0


def parser_vcf_subset(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf refFasta outVcf",
        version=version,
        description='''Filter VCF using SelectVariants''')
    parser.add_option("--selectVariantsOptions", dest="selectVariantsOptions", type='string',
                      help="Options for SelectVariants.  Often things like -selectType SNP --sample_file include.txt --exclude_sample_file exclude.txt -ef  -restrictAllelesTo BIALLELIC -env, etc.",
                      default='')
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_vcf_subset(args, options):
    inVcf, refFasta, outVcf = args
    tmpVcf = util_files.mkstempfname(prefix='vcf_subset-', suffix='.vcf')
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T SelectVariants -U LENIENT_VCF_PROCESSING --variant %s -R %s -o %s %s" % (
        tempfile.tempdir, options.gatkPath+'/GenomeAnalysisTK.jar',
        inVcf, refFasta, tmpVcf, options.selectVariantsOptions)
    log.info("vcf_subset: %s" % cmdline)
    assert not os.system(cmdline)
    util_vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf)
    log.info("done")
    return 0


def parser_vcf_to_biallelic(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf refFasta outVcf",
        version=version,
        description='''Filter VCF down to markers we may want to use for GWAS:
Markers must be biallelic (SNPs and indels allowed) and must be non-singletons.
Also filter down to markers that have at most a 20% missing rate (--geno 0.8).
At some point, we should add an option to filter to SNPs only or indels only
(GATK has a method for this though).''')
    parser.add_option("--selectVariantsOptions", dest="selectVariantsOptions", type='string',
                      help="Extra options for SelectVariants first pass.  Often things like -selectType SNP --sample_file include.txt --exclude_sample_file exclude.txt -ef, etc.",
                      default='')
    parser.add_option("--vcftoolsOptions", dest="vcftoolsOptions", type='string',
                      help="Extra options for vcftools first pass.  Often things like --mind, --thin, etc.",
                      default='')
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_vcf_to_biallelic(args, options):
    inVcf, refFasta, outVcf = args
    tmpVcf1 = util_files.mkstempfname(prefix='vcf_to_biallelic-1-')
    tmpVcf2 = util_files.mkstempfname(prefix='vcf_to_biallelic-2-')

    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T SelectVariants -U LENIENT_VCF_PROCESSING --variant %s -R %s -o %s -restrictAllelesTo BIALLELIC -env %s" % (
        tempfile.tempdir, options.gatkPath+'/GenomeAnalysisTK.jar',
        inVcf, refFasta, tmpVcf1+'.vcf', options.selectVariantsOptions)
    log.info("vcf_to_biallelic: %s" % cmdline)
    assert not os.system(cmdline)

    # notes: the --mac 1 is necessary because the -env in SelectVariants above
    # doesn't actually test whether or not there's variation in the population
    # (it only looks a the ALT column) so this allows us to eliminate true
    # novars.  We don't bother with --mac 2 (and instead do this two-step thing
    # to identify and exclude singletons) because these VCFs are encoded as
    # diploids and --mac looks at the allele count (not the individual count)
    # so you'd need --mac 4 to remove singletons, but unfortunately that would
    # allow for minor alleles spread as hets throughout the population.
    # Luckily, the vcftools --singletons command also identifies "private
    # doubletons" which is exactly what we want to eliminate, so we do the
    # two-step process.
    cmdline = "%s/vcftools --gzvcf %s --out %s --geno 0.8 --recode --missing --mac 1 --singletons %s" % (
        options.vcftoolsPath, tmpVcf1+'.vcf', tmpVcf2, options.vcftoolsOptions)
    log.info("vcf_to_biallelic: %s" % cmdline)
    assert not os.system(cmdline)
    os.unlink(tmpVcf1+'.vcf')

    cmdline = "%s/vcftools --vcf %s --out %s --recode --exclude-positions %s.singletons" % (
        options.vcftoolsPath, tmpVcf2+'.recode.vcf', tmpVcf1, tmpVcf2)
    log.info("vcf_to_biallelic: %s" % cmdline)
    assert not os.system(cmdline)
    os.unlink(tmpVcf2+'.recode.vcf')

    util_vcf.vcf_bgzip_index(tmpVcf1+'.recode.vcf', outVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf1+'.recode.vcf')
    log.info("done")
    return 0


def parser_filter_sample_callrate(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inVcf refFasta outVcf",
        version=version,
        description='''Take a VCF representing consensus calls (whole genome)
and filter the VCF to those samples that have sufficient genome coverage.
The cutoff is an absolute number, not a percentage (like --mind in vcftools).''')
    parser.add_option("--minCall", dest="minCall", type='int',
                      help="Minimum number of bases required per sample.  [default: %default]",
                      default=23332831/2)
    parser.add_option("--statsOut", dest="statsOut", type='string',
                      help="Output file for genome coverage stats.",
                      default=None)
    parser.add_option("--gatkPath", dest="gatkPath", type='string',
                      help="Directory for GATK.  [default: %default]",
                      default=global_tool_paths['gatk'])
    parser.add_option("--vcftoolsPath", dest="vcftoolsPath", type='string',
                      help="Directory for vcftools.  [default: %default]",
                      default=global_tool_paths['vcftools'])
    parser.add_option("--tabixPath", dest="tabixPath", type='string',
                      help="Directory for tabix.  [default: %default]",
                      default=global_tool_paths['tabix'])
    parser.add_option("--tmpDir", dest="tmpDir", type='string',
                      help="Directory for temp files.  [default: %default]",
                      default=global_tool_paths['tmpDir'])
    return parser
def main_filter_sample_callrate(args, options):
    inVcf, refFasta, outVcf = args
    tmpStats = util_files.mkstempfname(prefix='stats-',suffix='.stats')
    #tmpBadSamps = util_files.mkstempfname(prefix='bad_samples-',suffix='.txt')
    tmpVcf = util_files.mkstempfname(prefix='filtered-',suffix='.vcf')

    cmdline = "%s/vcftools --gzvcf %s --out %s --missing" % (
        options.vcftoolsPath, inVcf, tmpStats)
    assert not os.system(cmdline)
    coverage = []
    bad_samps = []
    with open(tmpStats+'.imiss', 'rt') as inf:
        #with open(tmpBadSamps, 'wt') as outf:
            n=0
            for line in inf:
                row = line.rstrip('\n').split('\t')
                if row[0]!='INDV':
                    s, cvg = (row[0], int(row[1])-int(row[3]))
                    coverage.append((s,cvg))
                    if cvg<options.minCall:
                        bad_samps.append(s)
                        #outf.write(s+'\n')
                        log.info("dropping sample %s (only %d bp covered)" % (s,cvg))
                        n += 1
            log.info("dropped %d samples" % n)
    os.unlink(tmpStats+'.imiss')
    os.unlink(tmpStats+'.lmiss')
    os.unlink(tmpStats+'.log')

    if options.statsOut:
        with open(options.statsOut, 'wt') as outf:
            for s,cvg in coverage:
                outf.write("%s\t%d\n" % (s,cvg))

    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T SelectVariants --variant %s -R %s -o %s %s" % (
        tempfile.tempdir, options.gatkPath+'/GenomeAnalysisTK.jar',
        inVcf, refFasta, tmpVcf,
        ''.join([" --exclude_sample_name "+s for s in bad_samps]))
    #cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T SelectVariants --variant %s -R %s -o %s --exclude_sample_file %s" % (
    #   tempfile.tempdir, options.gatkPath+'/GenomeAnalysisTK.jar',
    #   inVcf, refFasta, tmpVcf, tmpBadSamps)
    assert not os.system(cmdline)
    #os.unlink(tmpBadSamps)

    util_vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=options.tabixPath, vcftoolsPath=options.vcftoolsPath)
    os.unlink(tmpVcf)
    log.info("done")
    return 0


def parser_sen_unique_patients(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inImiss outExcludeSamples",
        version=version,
        description='''Read through an .imiss file (produced by vcftools --missing).
For Senegal Broad samples, find samples that were sequenced multiple ways (c/d/h/hh)
and find the best performing run for each patient sample.  Produce an output list
of sample names to exclude that are *not* the best performing run for each sample.''')
    return parser
def main_sen_unique_patients(args, options):
    inImiss, outExcludeSamples = args
    patients = {}
    call_rates = {}
    with open(inImiss, 'rt') as inf:
        for line in inf:
            row = line.rstrip('\n\r').split('\t')
            if row[0]!='INDV':
                r = row[0]
                cr = 1.0 - float(row[4])
                call_rates[r] = cr
                if r.startswith('Sen'):
                    samp_parts = r.split('.')
                    assert len(samp_parts)==3 and samp_parts[2] in ('c','d','h','hh')
                    s = '.'.join(samp_parts[:2])
                    patients.setdefault(s, [])
                    patients[s].append(r)
    toss = []
    for p,runs in patients.items():
        if len(runs)>1:
            runs = sorted([(call_rates[r],r) for r in runs], reverse=True)
            log.info("preferring %s:%0.3f over %s" % (runs[0][1], runs[0][0], ', '.join(['%s:%0.3f'%(r,cr) for cr,r in runs[1:]])))
            toss = toss + [r for cr,r in runs][1:]
    with open(outExcludeSamples, 'wt') as outf:
        for r in toss:
            outf.write(r+'\n')
    return 0


def vcf_filter(inVcf, refFasta, outVcf,
    gatkPath=global_tool_paths['gatk'], tabixPath=global_tool_paths['tabix'], vcftoolsPath=global_tool_paths['vcftools']):
    ''' This uses GATK VariantFiltration to soft-filter low quality and het bases
        using the FILTER and FT fields.
    '''
    log.info("filtering VCF file for quality and homozygous genotypes")
    if outVcf.endswith('.gz'):
        tmpVcf = util_files.mkstempfname(prefix='vcf_filter-', suffix='.vcf')
    else:
        tmpVcf = outVcf
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T VariantFiltration --variant %s -R %s -o %s" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inVcf, refFasta, tmpVcf)
    cmdline += " -G_filterName 'het'"
    cmdline += " -G_filter 'isHet==1'"
    cmdline += " -G_filterName 'LowQual'"
    cmdline += " -G_filter 'GQ<30'"
    cmdline += " --filterName 'LowQual'"
    cmdline += " --filterExpression 'QUAL<60'"
    log.info("vcf_filter: %s" % cmdline)
    assert not os.system(cmdline)
    if outVcf.endswith('.gz'):
        util_vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=tabixPath, vcftoolsPath=vcftoolsPath)
        os.unlink(tmpVcf)
    return outVcf

def vcf_variants_only(inVcf, refFasta, outVcf,
    gatkPath=global_tool_paths['gatk'], tabixPath=global_tool_paths['tabix'], vcftoolsPath=global_tool_paths['vcftools']):
    ''' This uses GATK SelectVariants to hard-remove FILTERed loci
        (but not FT'ed genotypes) and removes all invariant loci as well
        (leaving SNPs, MNPs, indels, etc).
    '''
    log.info("filtering VCF file to variant sites")
    if outVcf.endswith('.gz'):
        tmpVcf = util_files.mkstempfname(prefix='vcf_variants_only-', suffix='.vcf')
    else:
        tmpVcf = outVcf
    cmdline = "java -Xmx2g -Djava.io.tmpdir=%s -jar %s -T SelectVariants --variant %s -R %s -o %s -ef -env" % (
        tempfile.tempdir, gatkPath+'/GenomeAnalysisTK.jar',
        inVcf, refFasta, tmpVcf)
    log.info("vcf_variants_only: %s" % cmdline)
    assert not os.system(cmdline)
    if outVcf.endswith('.gz'):
        util_vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=tabixPath, vcftoolsPath=vcftoolsPath)
        os.unlink(tmpVcf)
    return outVcf

def vcf_snpEff(inVcf, genome, outVcf, vcf_out=True,
    snpEffPath=global_tool_paths['snpEff'], tabixPath=global_tool_paths['tabix'], vcftoolsPath=global_tool_paths['vcftools']):
    log.info("adding snpEff annotation")
    if outVcf.endswith('.gz'):
        tmpVcf = util_files.mkstempfname(prefix='vcf_snpEff-', suffix='.vcf')
    else:
        tmpVcf = outVcf
    cmdline = "cat %s | java -Xmx2g -Djava.io.tmpdir=%s -jar %s/snpEff.jar eff -c %s/snpEff.config -o %s %s -treatAllAsProteinCoding false -noLog -ud 0 -noStats > %s" % (
        inVcf, tempfile.tempdir, snpEffPath, snpEffPath, (vcf_out and 'vcf' or 'txt'), genome, tmpVcf)
    if inVcf.endswith('.gz'):
        cmdline = 'z'+cmdline
    log.info("vcf_snpEff: %s" % cmdline)
    assert not os.system(cmdline)
    if outVcf.endswith('.gz'):
        if vcf_out:
            util_vcf.vcf_bgzip_index(tmpVcf, outVcf, tabixPath=tabixPath, vcftoolsPath=vcftoolsPath)
        else:
            cmdline = "%s/bgzip -c %s > %s" % (tabixPath, tmpVcf, outVcf)
            assert not os.system(cmdline)
            cmdline = "%s/tabix %s -f -p vcf" % (tabixPath, outVcf)
            assert not os.system(cmdline)
        os.unlink(tmpVcf)
    return outVcf


def parser_broad_bams(cmd, version):
    parser = optparse.OptionParser(
        usage="usage: python %prog "+cmd+" [options] inReport outReport linkDir",
        version=version,
        description='''Take Broad sequencing report and find current locations of
bam files.  Send output to a text file and add symlinks to a directory of symlinks.''')
    parser.add_option("--sampleCol", dest="sampleCol", type='string',
                      help="""Column name in the inReport that corresponds to the real
sample name we want to use from here on out.  [default: %default]""",
                      default='Sample')
    parser.add_option("--queryKeys", dest="queryKeys", type='string',
                      help="""Which query keys do we use to extract BAM files from BASS
using information in the inReport file?  We assume the inReport file is a tab text file
with a single row.  The queryKeys parameter is a semicolon delimited list of key pairs
that takes the form of key1:col1;key2:col2;key3:col3;etc.  The "key" is the parameter
we pass to BASS (see dmsClient -help -attribute_list for a list of parameters).  The
"col" is the column header (1st line) of the inReport.  All other columns of the
inReport are ignored.  This script will error if any row of the inReport fails to
return exactly one BAM file.  No outputs will be written until we know that all
rows succeed.  [default: %default]""",
                      default='initiative:Initiative;project:Project;sample:External ID')
    parser.add_option("--bassPath", dest="bassPath", type='string',
                      help="Directory for BASS command line clients.  [default: %default]",
                      default=global_tool_paths['bass'])
    return parser
def main_broad_bams(args, options):
    inReport, outReport, linkDir = args
    sample_to_bam = {}

    log.info("querying for bam files")
    with open(inReport, 'rt') as inf:
        tmpFile = util_files.mkstempfname(suffix='.txt', prefix='broad_bams-')
        queryKeys = [pair.split(':') for pair in options.queryKeys.split(';')]
        for row in util_files.FlatFileParser(inf):
            sample = row[options.sampleCol]
            assert not os.access('%s/%s.bam' % (linkDir, sample), os.R_OK), "error: %s symlink already exists" % sample
            assert sample not in sample_to_bam, "%s seen twice in %s!" % (sample, inReport)
            keys = ' '.join(["-%s '%s'" % (k,row[c]) for k,c in queryKeys])
            cmdline = "%s/dmsList -type bam -report files %s > %s" % (options.bassPath, keys, tmpFile)
            assert not os.system(cmdline)
            with open(tmpFile, 'rt') as intmp:
                bams = [l.rstrip('\r\n') for l in intmp]
            assert len(bams)==1, "error: %d bam files returned on query for sample %s!  query: %s  bam files: %s" % (len(bams), sample, keys, str(bams))
            assert bams[0].endswith('.bam'), "error: %s bam file %s does not end with bam!" % (sample, bams[0])
            assert os.access(bams[0], os.R_OK), "error: %s bam file %s is not readable" % (sample, bams[0])
            sample_to_bam[sample] = bams[0]
    log.info("found %d bam files" % len(sample_to_bam))

    with open(outReport, 'wt') as outf:
        for s in sorted(sample_to_bam.keys()):
            b = sample_to_bam[s]
            os.symlink(b, '%s/%s.bam' % (linkDir, s))
            outf.write('%s\t%s\n' % (s, b))
    log.info("done")
    return 0


######## misc utility stuff below

def unique(items):
    ''' Return unique items in the same order as seen in the input. '''
    seen = set()
    out = []
    for i in items:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out

if __name__ == '__main__':
    main()


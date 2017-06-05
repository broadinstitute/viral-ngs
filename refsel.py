#!/usr/bin/env python
''' Selecting, for each sample, the reference closest to that sample.  "refsel" refers to "reference selection"
throughout.
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import logging
import random
import os
import os.path
import re
import shutil
import subprocess

try:
    from itertools import zip_longest    # pylint: disable=E0611
except ImportError:
    from itertools import izip_longest as zip_longest    # pylint: disable=E0611

# intra-module
import util.cmd
import util.file
import util.genbank
import util.misc
import util.vcf
import read_utils
import tools
import tools.picard
import tools.samtools
import tools.gatk
import tools.novoalign
import tools.trimmomatic
import tools.trinity
import tools.mafft
import tools.mummer
import tools.muscle
import tools.kmc
import tools.clark

# third-party
import Bio.AlignIO
import Bio.SeqIO, Bio.Seq, Bio.SeqRecord, Bio.Alphabet
import Bio.Data.IUPACData

log = logging.getLogger(__name__)

def build_refsel_db_with_kmc(refselDir, refFasta, kmcKmerSize=35, clarkKmerSize=35,
                             JVMmemory=None, threads=1,novoalign_license_path=None):
    '''Precompute info that will help us identify, for each sample, a reference genome closest to that sample.

    Args:
       refselDir: directory where refsel data will be stored.  It must contain one file, named refsel.fasta,
           containing a list of possible references to select from.  If the reference genome consists of 
           N segments, the first N entries of refsel.fasta give one genome, the next N entries give another
           genome, and so on.
       refFasta: a fasta file giving the standard reference genome.
    '''

    refselFN=os.path.join(refselDir,'refsel.fasta')

    refsel_ids_seen=set()
    for r in Bio.SeqIO.parse(refselFN, 'fasta'):
        assert r.id not in refsel_ids_seen, 'Duplicate entry {} in {}'.format(r.id, refselFN)
        refsel_ids_seen.add(r.id)

    doneFile=refselFN+'.indexing.done'
    if os.path.isfile(doneFile): os.unlink(doneFile)

    nsegs=util.file.fasta_length(refFasta)

    shutil.copyfile(refFasta, os.path.join(refselDir, 'reference.fasta'))


    read_utils.index_fasta_all(refselFN, JVMmemory=JVMmemory, threads=threads, 
                               novoalign_license_path=novoalign_license_path)

    clarkDir=os.path.join(refselDir,'CLARK', 'k'+str(clarkKmerSize))
    clarkDbDir=os.path.join(clarkDir,'clarkDb')
    util.file.mkdir_p(clarkDbDir)

    refsDir=os.path.join(refselDir,'refs')
    util.file.mkdir_p(refsDir)
    refsSplitDir=os.path.join(clarkDir,'refs')
    util.file.mkdir_p(refsSplitDir)
    kmcDir=os.path.join(refselDir,'kmc','k'+str(kmcKmerSize))
    util.file.mkdir_p(kmcDir)

    kmc = tools.kmc.KmcTool()

    refselKmcDb=os.path.join(kmcDir, 'refsel.kmc')

    kmc.build_kmer_db(kmerDb=refselKmcDb,
                      inFiles=refselFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
    kmc.histogram(refselKmcDb, refselKmcDb+'.histogram.txt')
    kmc.dump(refselKmcDb, refselKmcDb+'.kmers.txt')
    for cx in (1,2,3):
        kmc.reduce(kmerDbIn=refselKmcDb, kmerDbOut=os.path.join(kmcDir, 'refsel.cx{}.kmc'.format(cx)),
                   kmerDbIn_opts='-ci1 -cx{}'.format(cx))

    clarkTargsFN=os.path.join(clarkDir,'targets_addresses.txt')

    segsRecs = [ [] for seg_idx in range(nsegs) ]

    refListFN=os.path.join(refselDir, 'reflist.txt')

    with open(clarkTargsFN, 'wt') as outTargAddr, open(refListFN, 'wt') as refListOut:
        recs = list(Bio.SeqIO.parse(refselFN, 'fasta'))
        assert (len(recs) % nsegs)==0, "Number of ref seqs not a multiple of the number of ref segments"
        nrefs = int(len(recs) / nsegs)

        log.info('Processing {} scaffolding references, with {} segments each'.format(nrefs,nsegs))

        for ref_num in range(nrefs):
            ref_idx_in_list = ref_num*nsegs
            # we'll denote a ref by the id of its first segment
            ref_id = recs[ref_idx_in_list].id

            log.info('Processing ref {} ({}/{})'.format(ref_id,ref_num,nrefs))
            
            refFNbase=util.file.string_to_file_name(ref_id)
            refFN=os.path.join(refsDir, refFNbase+'.fa')
            Bio.SeqIO.write(recs[ref_idx_in_list:ref_idx_in_list+nsegs], refFN, 'fasta')
            kmc.build_kmer_db(kmerDb=os.path.join(kmcDir, refFNbase+'.kmc'),
                              inFiles=refFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
            
            segFNs=[]
            for seg_idx in range(nsegs):
                seg_rec = recs[ref_idx_in_list+seg_idx]
                segsRecs[seg_idx].append(seg_rec)
                seg_id = seg_rec.id


                segFN=os.path.join(refsSplitDir,util.file.string_to_file_name(seg_id) + '.fa')
                segFNs.append(segFN)
                log.info('Procesing SEG {} ref_num={} ref_idx_in_list={} seg_id={} segFN={}'.format(seg_idx, ref_num, ref_idx_in_list, seg_id, segFN))
                Bio.SeqIO.write([seg_rec], segFN, 'fasta')
                targId='REF{:04d}_{}'.format(ref_num,refFNbase)
                targId=targId[:25]
                outTargAddr.write( segFN + ' ' + targId + '\n' )

            refListOut.write('\t'.join([refFN] + segFNs) + '\n')


    for seg_idx in range(nsegs):
        segFasta=os.path.join(refselDir,'segs{}.fasta'.format(seg_idx))
        Bio.SeqIO.write(segsRecs[seg_idx], segFasta, 'fasta')
        read_utils.index_fasta_all(segFasta, JVMmemory=JVMmemory, threads=threads, 
                                   novoalign_license_path=novoalign_license_path)

    clark = tools.clark.ClarkTool()

    clark.execute(['-T', clarkTargsFN,
                   '-D', clarkDbDir+'/',
                   '-O', os.path.join(refselDir, 'refsel.fasta'),
                   '-m', '0', '--extended',
                   '-R', os.path.join(clarkDir, 'refsel.test')])

    shutil.copyfile(refselFN, doneFile)
    
    # kmc = tools.kmc.KmcTool()
    
    # kmc.builds_kmer_db( inFiles = [refsFasta], 
    #                     kmerDb = kmerDb, kmerSize=kmerSize, kmerOpts=kmerOpts, threads=threads )


def build_refsel_db(refselDir, refFasta, kmcKmerSize=35, clarkKmerSize=35,
                    JVMmemory=None, threads=1,novoalign_license_path=None):
    '''Precompute info that will help us identify, for each sample, a reference genome closest to that sample.

    Args:
       refselDir: directory where refsel data will be stored.  It must contain one file, named refsel.fasta,
           containing a list of possible references to select from.  If the reference genome consists of 
           N segments, the first N entries of refsel.fasta give one genome, the next N entries give another
           genome, and so on.
       refFasta: a fasta file giving the standard reference genome.
    '''

    refselFN=os.path.join(refselDir,'refsel.fasta')

    refsel_ids_seen=set()
    for r in Bio.SeqIO.parse(refselFN, 'fasta'):
        assert r.id not in refsel_ids_seen, 'Duplicate entry {} in {}'.format(r.id, refselFN)
        refsel_ids_seen.add(r.id)

    doneFile=refselFN+'.indexing.done'
    if os.path.isfile(doneFile): os.unlink(doneFile)

    nsegs=util.file.fasta_length(refFasta)

    shutil.copyfile(refFasta, os.path.join(refselDir, 'reference.fasta'))


    read_utils.index_fasta_all(refselFN, JVMmemory=JVMmemory, threads=threads, 
                               novoalign_license_path=novoalign_license_path)

    clarkDir=os.path.join(refselDir,'CLARK', 'k'+str(clarkKmerSize))
    clarkDbDir=os.path.join(clarkDir,'clarkDb')
    util.file.mkdir_p(clarkDbDir)

    refsDir=os.path.join(refselDir,'refs')
    util.file.mkdir_p(refsDir)
    refsSplitDir=os.path.join(clarkDir,'refs')
    util.file.mkdir_p(refsSplitDir)
    kmcDir=os.path.join(refselDir,'kmc','k'+str(kmcKmerSize))
    util.file.mkdir_p(kmcDir)

    kmc = tools.kmc.KmcTool()

    refselKmcDb=os.path.join(kmcDir, 'refsel.kmc')

    kmc.build_kmer_db(kmerDb=refselKmcDb,
                      inFiles=refselFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
    kmc.histogram(refselKmcDb, refselKmcDb+'.histogram.txt')
    kmc.dump(refselKmcDb, refselKmcDb+'.kmers.txt')
    for cx in (1,2,3):
        kmc.reduce(kmerDbIn=refselKmcDb, kmerDbOut=os.path.join(kmcDir, 'refsel.cx{}.kmc'.format(cx)),
                   kmerDbIn_opts='-ci1 -cx{}'.format(cx))

    clarkTargsFN=os.path.join(clarkDir,'targets_addresses.txt')

    segsRecs = [ [] for seg_idx in range(nsegs) ]

    refListFN=os.path.join(refselDir, 'reflist.txt')

    with open(clarkTargsFN, 'wt') as outTargAddr, open(refListFN, 'wt') as refListOut:
        recs = list(Bio.SeqIO.parse(refselFN, 'fasta'))
        assert (len(recs) % nsegs)==0, "Number of ref seqs not a multiple of the number of ref segments"
        nrefs = int(len(recs) / nsegs)

        log.info('Processing {} scaffolding references, with {} segments each'.format(nrefs,nsegs))

        for ref_num in range(nrefs):
            ref_idx_in_list = ref_num*nsegs
            # we'll denote a ref by the id of its first segment
            ref_id = recs[ref_idx_in_list].id

            log.info('Processing ref {} ({}/{})'.format(ref_id,ref_num,nrefs))
            
            refFNbase=util.file.string_to_file_name(ref_id)
            refFN=os.path.join(refsDir, refFNbase+'.fa')
            Bio.SeqIO.write(recs[ref_idx_in_list:ref_idx_in_list+nsegs], refFN, 'fasta')
            if False:
                kmc.build_kmer_db(kmerDb=os.path.join(kmcDir, refFNbase+'.kmc'),
                                  inFiles=refFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
            
            segFNs=[]
            for seg_idx in range(nsegs):
                seg_rec = recs[ref_idx_in_list+seg_idx]
                segsRecs[seg_idx].append(seg_rec)
                seg_id = seg_rec.id


                segFN=os.path.join(refsSplitDir,util.file.string_to_file_name(seg_id) + '.fa')
                segFNs.append(segFN)
                log.info('Procesing SEG {} ref_num={} ref_idx_in_list={} seg_id={} segFN={}'.format(seg_idx, ref_num, ref_idx_in_list, seg_id, segFN))
                Bio.SeqIO.write([seg_rec], segFN, 'fasta')
                targId='REF{:04d}_{}'.format(ref_num,refFNbase)
                targId=targId[:25]
                outTargAddr.write( segFN + ' ' + targId + '\n' )

            refListOut.write('\t'.join([refFN] + segFNs) + '\n')


    for seg_idx in range(nsegs):
        segFasta=os.path.join(refselDir,'segs{}.fasta'.format(seg_idx))
        Bio.SeqIO.write(segsRecs[seg_idx], segFasta, 'fasta')
        read_utils.index_fasta_all(segFasta, JVMmemory=JVMmemory, threads=threads, 
                                   novoalign_license_path=novoalign_license_path)

    clark = tools.clark.ClarkTool()

    clark.execute(['-T', clarkTargsFN,
                   '-D', clarkDbDir+'/',
                   '-O', os.path.join(refselDir, 'refsel.fasta'),
                   '-m', '0', '--extended',
                   '-R', os.path.join(clarkDir, 'refsel.test')])

    shutil.copyfile(refselFN, doneFile)
    
    # kmc = tools.kmc.KmcTool()
    
    # kmc.builds_kmer_db( inFiles = [refsFasta], 
    #                     kmerDb = kmerDb, kmerSize=kmerSize, kmerOpts=kmerOpts, threads=threads )

    
def parser_build_refsel_db(parser=argparse.ArgumentParser()):
    parser.add_argument('refselDir', help='Directory where refsel info is stored.')
    parser.add_argument('refFasta', help='Main fasta reference.')
    util.cmd.attach_main(parser, build_refsel_db, split_args=True)
    return parser


__commands__.append(('build_refsel_db', parser_build_refsel_db))


def select_ref(refselDir, contigs, out_ref, clarkKmerSize=35):
    '''Choose a reference closest to the sample.
    '''
    
    clark = tools.clark.ClarkTool()
    clarkDir=os.path.join(refselDir,'CLARK', 'k'+str(clarkKmerSize))
    clarkDbDir=os.path.join(clarkDir,'clarkDb')
    clarkTargsFN=os.path.join(clarkDir,'targets_addresses.txt')

    with util.file.tempfnames(('refsel_contigs_joined', 'refsel_clark_results')) as (contigs_joined, clark_results):
        Bio.SeqIO.write([Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(('N'*clarkKmerSize).join([str(r.seq) for r in Bio.SeqIO.parse(contigs, 'fasta')]), Bio.Alphabet.generic_dna), id='joined_contigs')], contigs_joined, 'fasta')

        log.info('running clark...')
        clark.execute(['-T', clarkTargsFN,
                       '-D', clarkDbDir+'/',
                       '-O', contigs_joined,
                       '-m', '0', '--extended',
                       '-R', clark_results])
        log.info('clark done')
        shutil.copyfile(clark_results+'.csv', out_ref+'.results.csv')
        closest_ref_contents=util.file.slurp_file(clark_results+'.csv')
        closest_ref_lines=closest_ref_contents.split('\n')
        closest_ref_data_line=closest_ref_lines[1]
        closest_ref_name_list=closest_ref_data_line.split(',')
        closest_ref_name=closest_ref_name_list[-5]
        found=False
        for target_addr in util.file.slurp_file(clarkTargsFN).split('\n'):
            refFasta, refName = target_addr.strip().split()
            if refName==closest_ref_name:
                shutil.copyfile(refFasta, out_ref)
                log.info('Copying best ref {} to {}'.format(refFasta, out_ref))
                found=True
                break
        assert found
    
def parser_select_ref(parser=argparse.ArgumentParser(description='Select reference closest to sample')):
    parser.add_argument('refselDir', help='Directory where refsel info was stored by the build_refsel_db command.')
    parser.add_argument('contigs', help='Contigs assembled from sample.')
    parser.add_argument('out_ref', help='Output the selected reference to this file.')
    util.cmd.attach_main(parser, select_ref, split_args=True)
    return parser

__commands__.append(('select_ref', parser_select_ref))


def kmc_build_db(inFile, kmerDb, kmerSize, kmerOpts, histogram='', dump=''):
    '''Build a kmc kmer database'''

    kmc = tools.kmc.KmcTool()
    kmc.build_kmer_db(inFiles=[inFile], kmerDb=kmerDb, kmerSize=kmerSize, kmerOpts=kmerOpts)
    if histogram: kmc.histogram(kmerDb, histogram)
    if dump: kmc.dump(kmerDb, dump)
    
    
def parser_kmc_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inFile', help='Bam or fasta file to go into the database.')
    parser.add_argument('kmerDb', help='kmer db name.')
    parser.add_argument('kmerSize', help='kmer size.')
    parser.add_argument('--kmerOpts', help='kmer opts.')
    parser.add_argument('--histogram', help='also output histogram.')
    parser.add_argument('--dump', help='also output dump.')
    util.cmd.attach_main(parser, kmc_build_db, split_args=True)
    return parser


__commands__.append(('kmc_build_db', parser_kmc_build_db))


#def subtract_reads_from_refs(scafrefDir, readsKmcDb):
#    '''For each ref, see which kmers in that ref are _not_ in the db'''

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)

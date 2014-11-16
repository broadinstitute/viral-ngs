## note: this is tailored to Lassa, see ebov_assembly.sh for Ebola-specific tweaks

#-------- VIRAL GENOME ASSEMBLY @ BROAD --------#
use BLAST+
use Perl-5.10
use Samtools
use Python-2.7
use BWA
use BamTools
use Java-1.7

# MAKE REQUIRED SUB-DIRECTORIES
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o ~/log.txt -P sabeti_align -J Directories "mkdir $directory/_logs $directory/_temp $directory/_reads $directory/_bams $directory/_reports $directory/_pileup $directory/_meta $directory/_refs"
done

# CONVERT BAM FILES TO FASTQ
for sample in
do
for directory in
do
for suffix in sorted.bam
do
for in_url in /idi/sabeti-data/kandersen/analysis/140416_automation-samples/final_files
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.fq "java -Xmx2g -jar /seq/software/picard/current/bin/SamToFastq.jar INPUT=$in_url/$sample.$suffix FASTQ=$directory/_reads/$sample.reads1.fastq SECOND_END_FASTQ=$directory/_reads/$sample.reads2.fastq VALIDATION_STRINGENCY=SILENT"
done
done
done
done

#-------- PERFORM DE NOVO ASSEMBLY OF LASSA READS USING TRINITY --------#
# TRIM THE READS WITH TRIMMOMATIC
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=4]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.tr "java -Xmx2g -classpath /idi/sabeti-scratch/kandersen/bin/trimmomatic/trimmomatic-0.32.jar org.usadellab.trimmomatic.TrimmomaticPE $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq $directory/_reads/$sample.trimmed.1.fastq $directory/_temp/$sample.reads1.trimmed_unpaired.fastq $directory/_reads/$sample.trimmed.2.fastq $directory/_temp/$sample.reads2.trimmed_unpaired.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:30 ILLUMINACLIP:/idi/sabeti-scratch/kandersen/references/contaminants/contaminants.fasta:2:30:12 && rm $directory/_temp/$sample.reads?.trimmed_unpaired.fastq"
done
done

# LASTAL ALIGNMENT TO FIND RELEVANT READS
for sample in
do
for directory in
do
for database in arena
do
bsub -n 1 -R "span[hosts=1]" -q week -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.l1 "/idi/sabeti-scratch/kandersen/bin/last/lastal -Q1 /idi/sabeti-scratch/kandersen/references/lastal/$database $directory/_reads/$sample.trimmed.1.fastq | /idi/sabeti-scratch/kandersen/bin/last/scripts/maf-sort.sh -n2 | /idi/sabeti-scratch/kandersen/bin/last/scripts/maf-convert.py tab /dev/stdin > $directory/_temp/$sample.reads1.lastal.txt && python /idi/sabeti-scratch/kandersen/bin/scripts/noBlastLikeHits.py -b $directory/_temp/$sample.reads1.lastal.txt -r $directory/_reads/$sample.trimmed.1.fastq -m hit | perl /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -ns_max_n 1 -derep 1 -fastq stdin -out_bad null -line_width 0 -out_good $directory/_temp/$sample.lastal.1 && rm $directory/_temp/$sample.reads1.lastal.txt"
bsub -n 1 -R "span[hosts=1]" -q week -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.l2 "/idi/sabeti-scratch/kandersen/bin/last/lastal -Q1 /idi/sabeti-scratch/kandersen/references/lastal/$database $directory/_reads/$sample.trimmed.2.fastq | /idi/sabeti-scratch/kandersen/bin/last/scripts/maf-sort.sh -n2 | /idi/sabeti-scratch/kandersen/bin/last/scripts/maf-convert.py tab /dev/stdin > $directory/_temp/$sample.reads2.lastal.txt && python /idi/sabeti-scratch/kandersen/bin/scripts/noBlastLikeHits.py -b $directory/_temp/$sample.reads2.lastal.txt -r $directory/_reads/$sample.trimmed.2.fastq -m hit | perl /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -ns_max_n 1 -derep 1 -fastq stdin -out_bad null -line_width 0 -out_good $directory/_temp/$sample.lastal.2 && rm $directory/_temp/$sample.reads2.lastal.txt"
done
done
done

# DE NOVO ASSEMBLY USING TRINITY
use Java-1.6
for sample in
do
for directory in
do
bsub -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -q hour -W 4:00 -o $directory/_logs/$sample.log.denovo.txt -P sabeti_meta -J $sample.dn "/idi/sabeti-scratch/kandersen/bin/scripts/mergeShuffledFastqSeqs.pl -t -r '^@(\S+)/[1|2]$' -f1 $directory/_temp/$sample.lastal.1.fastq -f2 $directory/_temp/$sample.lastal.2.fastq -o $directory/_temp/$sample.clean && /idi/sabeti-scratch/kandersen/bin/scripts/subsampler.py -n 100000 -mode p -in $directory/_temp/$sample.clean.1.fastq $directory/_temp/$sample.clean.2.fastq -out $directory/_temp/$sample.reads1.sub.fastq $directory/_temp/$sample.reads2.sub.fastq && wc -l $directory/_temp/$sample.clean.?.fastq > $directory/_logs/$sample.log.lastal.txt && perl /idi/sabeti-scratch/kandersen/bin/trinity_old/Trinity.pl --CPU 1 --min_contig_length 300 --seqType fq --left $directory/_temp/$sample.reads1.sub.fastq --right $directory/_temp/$sample.reads2.sub.fastq --output $directory/_temp/$sample.trinity && mv $directory/_temp/$sample.trinity/Trinity.fasta $directory/_pileup/$sample.contigs.fasta && rm $directory/_temp/$sample.clean.nomatch.fastq && rm $directory/_temp/$sample.reads?.sub.fastq"
done
done

# ALIGN AND ORIENT CONTIGS TO REFERENCES
for sample in
do
for directory in
do
for country in NG # SL
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c1 "touch $directory/_temp/$sample.s_segment1_merger_assembly.fa && perl /seq/viral/analysis/xyang/scripts/others/VfatSoftwarePackage_201401/orientContig.pl $directory/_pileup/$sample.contigs.fasta /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.s.fasta $directory/_temp/$sample.s_segment1 && perl /seq/viral/analysis/xyang/scripts/others/VfatSoftwarePackage_201401/contigMerger.pl $directory/_temp/$sample.s_segment1_orientedContigs /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.s.fasta -readfq $directory/_temp/$sample.clean.1.fastq -readfq2 $directory/_temp/$sample.clean.2.fastq -fakequals 30 $directory/_temp/$sample.s_segment1 && cat $directory/_temp/$sample.s_segment1*assembly.fa > $directory/_temp/$sample.s_segment1.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c2 "touch $directory/_temp/$sample.l_segment1_merger_assembly.fa && perl /seq/viral/analysis/xyang/scripts/others/VfatSoftwarePackage_201401/orientContig.pl $directory/_pileup/$sample.contigs.fasta /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.l.fasta $directory/_temp/$sample.l_segment1 && perl /seq/viral/analysis/xyang/scripts/others/VfatSoftwarePackage_201401/contigMerger.pl $directory/_temp/$sample.l_segment1_orientedContigs /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.l.fasta -readfq $directory/_temp/$sample.clean.1.fastq -readfq2 $directory/_temp/$sample.clean.2.fastq -fakequals 30 $directory/_temp/$sample.l_segment1 && cat $directory/_temp/$sample.l_segment1*assembly.fa > $directory/_temp/$sample.l_segment1.fasta"
done
done
done

# SORT CONTIGS BASED ON LENGTH
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.s1 "python /idi/sabeti-data/rsealfon/assembly_scripts/filter_short_seqs.py $directory/_temp/$sample.s_segment1.fasta 2800 $directory/_temp/$sample.s_segment1.bad.fasta $directory/_temp/$sample.s_segment1.good.fasta && rm $directory/_temp/$sample.*.fa $directory/_temp/$sample.*.txt $directory/_temp/$sample.*.qlx $directory/_temp/$sample.*.pdf $directory/_temp/$sample.*.R $directory/_temp/$sample.*.afa $directory/_temp/$sample.*.mfa"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.s2 "python /idi/sabeti-data/rsealfon/assembly_scripts/filter_short_seqs.py $directory/_temp/$sample.l_segment1.fasta 6000 $directory/_temp/$sample.l_segment1.bad.fasta $directory/_temp/$sample.l_segment1.good.fasta"
done
done

#-------- MAP LASSA READS TO MODIFIED CONTIGS USING NOVOALIGN --------#
# ALIGN REFERENCES AND CONTIGS THAT PASSED LENGTH FILTER
for sample in
do
for directory in
do
for country in NG # SL
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.m1 "cat $directory/_temp/$sample.s_segment1.good.fasta /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.s.fasta | /idi/sabeti-scratch/kandersen/bin/muscle/muscle -out $directory/_temp/$sample.s_alignment1.fasta -quiet"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.m2 "cat $directory/_temp/$sample.l_segment1.good.fasta /idi/sabeti-scratch/kandersen/references/annotations/lasv/$country.l.fasta | /idi/sabeti-scratch/kandersen/bin/muscle/muscle -out $directory/_temp/$sample.l_alignment1.fasta -quiet"
done
done
done

# CLEANUP CONTIGS
for sample in
do
for directory in
do
for country in NG # SL
do
for year_species in 2013H
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c1 "python /idi/sabeti-data/rsealfon/assembly_scripts/modified_contig/modify_contig.py -n LASV-$sample-S-$country-$year_species --call-reference-ns yes --trim-ends yes --replace-5ends yes --replace-3ends yes --replace-length 20 --replace-end-gaps yes --remove-end-ns no --call-reference-ambiguous no $directory/_temp/$sample.s_alignment1.fasta $directory/_temp/$sample.s_segment2.fasta $country.s && python ~dpark/dev/sabetilab-dpark/scripts/viral.py deambig_fasta $directory/_temp/$sample.s_segment2.fasta $directory/_temp/$sample.s_segment2a.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.c2 "python /idi/sabeti-data/rsealfon/assembly_scripts/modified_contig/modify_contig.py -n LASV-$sample-L-$country-$year_species --call-reference-ns yes --trim-ends yes --replace-5ends yes --replace-3ends yes --replace-length 20 --replace-end-gaps yes --remove-end-ns no --call-reference-ambiguous no $directory/_temp/$sample.l_alignment1.fasta $directory/_temp/$sample.l_segment2.fasta $country.l && python ~dpark/dev/sabetilab-dpark/scripts/viral.py deambig_fasta $directory/_temp/$sample.l_segment2.fasta $directory/_temp/$sample.l_segment2a.fasta"
done
done
done
done

# INDEX MODIFIED CONTIGS
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_kga -J $sample.i1 "java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.s_segment2.fasta O=$directory/_temp/$sample.s_segment2.dict && /idi/sabeti-scratch/kandersen/bin/novocraft/novoindex $directory/_temp/$sample.s_segment2.nix $directory/_temp/$sample.s_segment2.fasta && samtools faidx $directory/_temp/$sample.s_segment2.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.s_segment2a.fasta O=$directory/_temp/$sample.s_segment2a.dict && samtools faidx $directory/_temp/$sample.s_segment2a.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_kga -J $sample.i2 "java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.l_segment2.fasta O=$directory/_temp/$sample.l_segment2.dict && /idi/sabeti-scratch/kandersen/bin/novocraft/novoindex $directory/_temp/$sample.l_segment2.nix $directory/_temp/$sample.l_segment2.fasta && samtools faidx $directory/_temp/$sample.l_segment2.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.l_segment2a.fasta O=$directory/_temp/$sample.l_segment2a.dict && samtools faidx $directory/_temp/$sample.l_segment2a.fasta"
done
done

# NOVOALIGN LASSA READS TO MODIFIED CONTIGS
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a1 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -f $directory/_temp/$sample.clean.1.fastq $directory/_temp/$sample.clean.2.fastq -r Random -l 30 -g 40 -x 20 -t 502 -F STDFQ -d $directory/_temp/$sample.s_segment2.nix -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' | samtools view -buS -q 1 - | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.ref_mapped_s1.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $directory/_temp/$sample.s_segment2a.fasta -I $directory/_temp/$sample.ref_mapped_s1.bam -o $directory/_temp/$sample.gatk_s1.vcf --baq OFF --useOriginalQualities -out_mode EMIT_ALL_SITES -dt NONE --num_threads 1 --min_base_quality_score 15 -ploidy 4 -stand_call_conf 0 -stand_emit_conf 0 -A AlleleBalance && python ~dpark/dev/sabetilab-dpark/scripts/viral.py vcf_to_fasta --trim_ends --min_coverage 2 $directory/_temp/$sample.gatk_s1.vcf $directory/_temp/$sample.s_segment3.fasta && python ~dpark/dev/sabetilab-dpark/scripts/viral.py deambig_fasta $directory/_temp/$sample.s_segment3.fasta $directory/_temp/$sample.s_segment3a.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a2 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -f $directory/_temp/$sample.clean.1.fastq $directory/_temp/$sample.clean.2.fastq -r Random -l 30 -g 40 -x 20 -t 502 -F STDFQ -d $directory/_temp/$sample.l_segment2.nix -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' | samtools view -buS -q 1 - | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.ref_mapped_l1.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $directory/_temp/$sample.l_segment2a.fasta -I $directory/_temp/$sample.ref_mapped_l1.bam -o $directory/_temp/$sample.gatk_l1.vcf --baq OFF --useOriginalQualities -out_mode EMIT_ALL_SITES -dt NONE --num_threads 1 --min_base_quality_score 15 -ploidy 4 -stand_call_conf 0 -stand_emit_conf 0 -A AlleleBalance && python ~dpark/dev/sabetilab-dpark/scripts/viral.py vcf_to_fasta --trim_ends --min_coverage 2 $directory/_temp/$sample.gatk_l1.vcf $directory/_temp/$sample.l_segment3.fasta && python ~dpark/dev/sabetilab-dpark/scripts/viral.py deambig_fasta $directory/_temp/$sample.l_segment3.fasta $directory/_temp/$sample.l_segment3a.fasta"
done
done

# INDEX CONSENSUS SEQUENCES
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_kga -J $sample.i1 "java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.s_segment3.fasta O=$directory/_temp/$sample.s_segment3.dict && /idi/sabeti-scratch/kandersen/bin/novocraft/novoindex $directory/_temp/$sample.s_segment3.nix $directory/_temp/$sample.s_segment3.fasta && samtools faidx $directory/_temp/$sample.s_segment3.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.s_segment3a.fasta O=$directory/_temp/$sample.s_segment3a.dict && samtools faidx $directory/_temp/$sample.s_segment3a.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_kga -J $sample.i2 "java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.l_segment3.fasta O=$directory/_temp/$sample.l_segment3.dict && /idi/sabeti-scratch/kandersen/bin/novocraft/novoindex $directory/_temp/$sample.l_segment3.nix $directory/_temp/$sample.l_segment3.fasta && samtools faidx $directory/_temp/$sample.l_segment3.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_temp/$sample.l_segment3a.fasta O=$directory/_temp/$sample.l_segment3a.dict && samtools faidx $directory/_temp/$sample.l_segment3a.fasta"
done
done

# NOVOALIGN ALL READS TO CONSENSUS SEQUENCES
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a1 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -f $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -r Random -l 40 -g 40 -x 20 -t 100 -F STDFQ -d $directory/_temp/$sample.s_segment3.nix -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' | samtools view -buS -q 1 - | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.ref_mapped_s2.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $directory/_temp/$sample.s_segment3a.fasta -I $directory/_temp/$sample.ref_mapped_s2.bam -o $directory/_temp/$sample.gatk_s2.vcf --baq OFF --useOriginalQualities -out_mode EMIT_ALL_SITES -dt NONE --num_threads 1 --min_base_quality_score 15 -ploidy 4 -stand_call_conf 0 -stand_emit_conf 0 -A AlleleBalance && python ~dpark/dev/sabetilab-dpark/scripts/viral.py vcf_to_fasta --trim_ends --min_coverage 2 $directory/_temp/$sample.gatk_s2.vcf $directory/_temp/$sample.s_segment4.fasta"
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.a2 "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -f $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -r Random -l 40 -g 40 -x 20 -t 100 -F STDFQ -d $directory/_temp/$sample.l_segment3.nix -o SAM $'@RG\tID:$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' | samtools view -buS -q 1 - | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.ref_mapped_l2.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $directory/_temp/$sample.l_segment3a.fasta -I $directory/_temp/$sample.ref_mapped_l2.bam -o $directory/_temp/$sample.gatk_l2.vcf --baq OFF --useOriginalQualities -out_mode EMIT_ALL_SITES -dt NONE --num_threads 1 --min_base_quality_score 15 -ploidy 4 -stand_call_conf 0 -stand_emit_conf 0 -A AlleleBalance && python ~dpark/dev/sabetilab-dpark/scripts/viral.py vcf_to_fasta --trim_ends --min_coverage 2 $directory/_temp/$sample.gatk_l2.vcf $directory/_temp/$sample.l_segment4.fasta"
done
done

#-------- PREPARE FINAL MAPPED ASSEMBLY USING NOVOALIGN --------#
# INDEX FINAL CONSENSUS SEQUENCES
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_kga -J $sample.ix "cat $directory/_temp/$sample.l_segment4.fasta $directory/_temp/$sample.s_segment4.fasta > $directory/_refs/$sample.fasta && java -Xmx2g -jar /seq/software/picard/current/bin/CreateSequenceDictionary.jar R=$directory/_refs/$sample.fasta O=$directory/_refs/$sample.dict && /idi/sabeti-scratch/kandersen/bin/novocraft/novoindex $directory/_refs/$sample.nix $directory/_refs/$sample.fasta && samtools faidx $directory/_refs/$sample.fasta"
done
done

# ALIGN ALL READS TO ITS OWN LASSA CONSENSUS
for sample in
do
for directory in
do
for date in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -n 4 -R "span[hosts=1]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.al "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -k -c 3 -f $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -r Random -l 40 -g 40 -x 20 -t 100 -F STDFQ -d $directory/_refs/$sample.nix -o SAM $'@RG\tID:$date.$sample\tSM:$sample\tPL:Illumina\tPU:HiSeq\tLB:BroadPE\tCN:Broad' 2> $directory/_logs/$sample.log.novoalign.txt | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_bams/$sample.sorted.bam CREATE_INDEX=true && samtools view -b -q 1 -u $directory/_bams/$sample.sorted.bam | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.mapped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && java -Xmx2g -jar /seq/software/picard/current/bin/MarkDuplicates.jar I=$directory/_temp/$sample.mapped.bam O=$directory/_temp/$sample.mappedNoDub.bam METRICS_FILE=$directory/_temp/$sample.log.markdups.txt CREATE_INDEX=true REMOVE_DUPLICATES=true"
done
done
done

#-------- PREPARE ALIGNMENT AND CALCULATE METRICS --------#
# RUN FASTQC REPORT
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.dx "/idi/sabeti-scratch/kandersen/bin/fastqc/fastqc -f bam $directory/_bams/$sample.sorted.bam -o $directory/_reports/"
done
done

# PERFORM LOCAL REALIGNMENT
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.lr "java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $directory/_refs/$sample.fasta -o $directory/_temp/$sample.log.realigner.intervals -I $directory/_temp/$sample.mappedNoDub.bam && java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T IndelRealigner -R $directory/_refs/$sample.fasta -targetIntervals $directory/_temp/$sample.log.realigner.intervals -I $directory/_temp/$sample.mappedNoDub.bam -o $directory/_bams/$sample.realigned.bam"
done
done

# CALCULATE COVERAGE & SNPS
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.st "/idi/sabeti-scratch/kandersen/bin/bedtools/genomeCoverageBed -d -ibam $directory/_bams/$sample.realigned.bam -g $directory/_refs/$sample.fasta > $directory/_pileup/$sample.coverage.txt && bamtools stats -insert -in $directory/_bams/$sample.realigned.bam > $directory/_logs/$sample.log.bamstats_realigned.txt && java -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $directory/_refs/$sample.fasta -I $directory/_bams/$sample.realigned.bam -o $directory/_temp/$sample.snps.vcf --baq OFF --useOriginalQualities -out_mode EMIT_VARIANTS_ONLY -dt NONE --num_threads 1 --min_base_quality_score 20 -ploidy 10 -stand_call_conf 80.0 -stand_emit_conf 30.0 -A AlleleBalance"
done
done

# FILTER SNPS
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.gk2 java -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T VariantFiltration -R $directory/_refs/$sample.fasta -V $directory/_temp/$sample.snps.vcf -o $directory/_pileup/$sample.snps.vcf -l ERROR --filterExpression "(MQ0/DP) > 0.20 || BaseQRankSum < -10.0 || QD < 1.0 || MQRankSum < -4.0 || ReadPosRankSum < -4.0" --filterName LowConfidence --filterExpression "DP < 20 || (DP*(1-ABHet)) < 5" --filterName LowCoverage --filterExpression "ABHet > 0.95" --filterName LowFrequency --filterExpression "SB < -1.0e+09 || SB > 1.0e+09 || FS > 5" --filterName StrandBias
done
done

# GET SUMMARY STATS FROM INDIVIDUAL SAMPLES
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.st "python /idi/sabeti-scratch/kandersen/bin/scripts/getcov.py $directory/_pileup/$sample.coverage.txt > $directory/_pileup/$sample.cov.txt"
done
done

#-------- SHORT META ANALYSIS @ BROAD --------#
# SUBSAMPLE READS FOR LATER RETRIEVAL
for sample in
do
for directory in
do
bsub -R "rusage[mem=6]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_meta -J $sample.ss "python /idi/sabeti-scratch/kandersen/bin/scripts/subsampler.py -n 1000000 -mode p -in $directory/_reads/$sample.reads1.fastq $directory/_reads/$sample.reads2.fastq -out $directory/_temp/$sample.reads1.sub.fastq $directory/_temp/$sample.reads2.sub.fastq"
done
done

# COUNT SPIKE-INS
for sample in
do
for directory in
do
bsub -W 4:00 -q hour -R "rusage[mem=2]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.si "/idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign -f $directory/_temp/$sample.reads1.sub.fastq $directory/_temp/$sample.reads2.sub.fastq -r Random -F STDFQ -d /idi/sabeti-scratch/kandersen/references/other/ercc_spike-ins.nix -o SAM 2> $directory/_temp/$sample.log.spike_novo.txt | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.spikes.bam CREATE_INDEX=true && samtools view -b -q 1 -u $directory/_temp/$sample.spikes.bam | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=$directory/_temp/$sample.spikes_mapped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT && /idi/sabeti-scratch/kandersen/bin/scripts/CountAlignmentsByDescriptionLine -bam $directory/_temp/$sample.spikes_mapped.bam > $directory/_logs/$sample.log.spike_count.txt"
done
done

# ALIGN TO HUMAN OR OTHER REFERENCE USING BWA
for sample in
do
for directory in
do
bsub -R "rusage[mem=8]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.hg "perl /idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -lc_method dust -lc_threshold 7 -line_width 0 -fastq $directory/_temp/$sample.reads1.sub.fastq -out_good stdout -out_bad null | /idi/sabeti-scratch/kandersen/bin/bwa/bwa aln -q 5 -l 32 /idi/sabeti-scratch/kandersen/references/human/hg19 - | /idi/sabeti-scratch/kandersen/bin/bwa/bwa samse /idi/sabeti-scratch/kandersen/references/human/hg19 - $directory/_temp/$sample.reads1.sub.fastq | java -Xmx2g -jar /seq/software/picard/current/bin/SortSam.jar SO=coordinate I=/dev/stdin O=/dev/stdout VALIDATION_STRINGENCY=SILENT | samtools flagstat - > $directory/_logs/$sample.log.hg19.txt"
done
done

# PERFORM BLASTN ANALYSIS
for sample in
do
for directory in
do
bsub -W 4:00 -R "rusage[mem=12]" -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.fa "/idi/sabeti-scratch/kandersen/bin/prinseq/prinseq-lite.pl -min_len 40 -seq_num 10000 -out_format 1 -line_width 0 -fastq $directory/_temp/$sample.reads1.sub.fastq -out_good $directory/_temp/$sample -out_bad null && blastn -db /idi/sabeti-scratch/kandersen/references/blast/nt -outfmt 6 -evalue 1e-10 -num_descriptions 30 -query $directory/_temp/$sample.fasta -out $directory/_temp/$sample.blastn.txt"
done
done

# CREATE MEGAN FILES
# You will have to be running an xhost in order to do this step: open /Applications/Utilities/X11.app/ xhost + and export DISPLAY=:0
for sample in
do
for directory in
do
bsub -R "rusage[mem=4]" -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.mg /idi/sabeti-scratch/kandersen/bin/megan/MEGAN -E +g -p "/idi/sabeti-scratch/kandersen/bin/megan/Megan.def" -x "load gi2taxfile="/idi/sabeti-scratch/kandersen/bin/megan/class/resources/files/gi_taxid_nucl.bin"; import blastfile="$directory/_temp/$sample.blastn.txt" fastafile="$directory/_temp/$sample.fasta" meganfile="$directory/_meta/$sample.rma" useseed=false usekegg=false; select nodes=all; uncollapse subtrees; update; exportimage file="$directory/_meta/$sample.pdf" format=pdf replace=true; quit;"
done
done
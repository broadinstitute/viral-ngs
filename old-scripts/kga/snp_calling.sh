#-------- SNP CALLING IN VIRAL SAMPLES @ BROAD --------#
# NO TRAILING / IN directory
use Perl-5.10
use Samtools
use Python-2.7
use BamTools
use Java-1.6

# MAKE REQUIRED SUB-DIRECTORIES - DON'T DO THIS IF YOU ALREADY CREATED THESE WITH ANOTHER PIPELINE
for directory in
do
bsub -o ~/log.txt -P sabeti_align -J Directories "mkdir $directory/_logs $directory/_temp $directory/_reads $directory/_bams $directory/_reports $directory/_pileup $directory/_meta"
done

#-------- VARSCAN --------#
# CALCULATE SNPS
for sample in
do
for directory in
do
for species in lassa
do
for reference in $sample
do
for suffix in realigned.bam
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.st "samtools mpileup -Q 0 -B -q 1 -d 10000 -f /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta $directory/_bams/$sample.$suffix | java -jar /idi/sabeti-scratch/kandersen/bin/varscan/varscan.jar pileup2snp --min-reads2 5 --min-var-freq 0.01 --p-value 0.1 --min-coverage 5 --min-avg-qual 5 > $directory/_temp/$sample.snps.txt"
done
done
done
done
done

# FILTER SNPS
for sample in
do
for directory in
do
for min_reads in 5
do
for min_frequency in 0.05
do
for strand_bias in 10
do
for min_quality in 25
do
for min_quality_difference in 5
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.st "perl /idi/sabeti-data/jesse/src/parseVarscan_bothStrands.pl $directory/_temp/$sample.snps.txt $min_reads $min_frequency $strand_bias $min_quality $min_quality_difference $sample"
done
done
done
done
done
done
done

# MOVE SNPS
for sample in
do
for directory in
do
mv $directory/*$sample.*.txt $directory/_pileup/$sample.varscan.snps.txt
done
done

#-------- GATK --------#
# CALCULATE SNPS WITH GATK
for sample in
do
for directory in
do
for species in lassa
do
for reference in $sample
do
for in_directory in $directory/_bams # /idi/sabeti-data/kandersen/sequencing_storage/lasv/positive/bam_mapped
do
for suffix in realigned.bam # merged.bam
do
for output_mode in EMIT_VARIANTS_ONLY # EMIT_ALL_SITES
do
for output_name in . # .merged
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.gk1 "java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta -I $in_directory/$sample.$suffix -o $directory/_temp/$sample$output_name.snps.vcf --baq OFF --useOriginalQualities -out_mode $output_mode -dt NONE --num_threads 1 --min_base_quality_score 20 -ploidy 10 -stand_call_conf 80.0 -stand_emit_conf 30.0 -A AlleleBalance"
done
done
done
done
done
done
done
done

# FILTER SNPS
for sample in
do
for directory in
do
for species in lassa
do
for reference in $sample
do
for output_name in . # .merged
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.gk2 java -Xmx2g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar -T VariantFiltration -R /idi/sabeti-scratch/kandersen/references/$species/$reference.fasta -V $directory/_temp/$sample$output_name.snps.vcf -o $directory/_temp/$sample$output_name.filtered_snps.vcf -l ERROR --filterExpression "(MQ0/DP) > 0.20 || BaseQRankSum < -10.0 || QD < 1.0 || MQRankSum < -4.0 || ReadPosRankSum < -4.0" --filterName LowConfidence --filterExpression "DP < 20 || (DP*(1-ABHet)) < 5" --filterName LowCoverage --filterExpression "ABHet > 0.95" --filterName LowFrequency --filterExpression "SB < -1.0e+09 || SB > 1.0e+09 || FS > 5" --filterName StrandBias
done
done
done
done
done

# MOVE SNPS
for sample in
do
for directory in
do
for output_name in . # .merged
do
for snp_name in .snps.vcf # .merged.snps.vcf
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.gk3 "mv $directory/_temp/$sample$output_name.filtered_snps.vcf $directory/_pileup/$sample$snp_name"
done
done
done
done

#-------- VPHASER --------#
# PREREQUISITES
bash
export LD_LIBRARY_PATH=/seq/viral/analysis/xyang/programs/Library/pezmaster31-bamtools-e235c55/lib:$LD_LIBRARY_PATH

# CALCULATE SNPS WITH VPHASER V2
for sample in
do
for directory in
do
for in_directory in $directory/_bams # /idi/sabeti-data/kandersen/sequencing_storage/lasv/positive/bam_mapped
do
for suffix in realigned.bam
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.vp1 "mkdir $directory/_temp/$sample && /seq/viral/analysis/xyang/programs/VariantCaller/bin/variant_caller -i $in_directory/$sample.$suffix -o $directory/_temp/$sample"
done
done
done
done

# MOVE RELEVANT FILES
for sample in
do
for directory in
do
bsub -W 4:00 -o $directory/_logs/$sample.log.bsub.txt -P sabeti_align -J $sample.vp2 "cat $directory/_temp/$sample/*.fdr.var.txt > $directory/_pileup/$sample.vp_snps.txt"
done
done

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


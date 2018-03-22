#task ref_guided_consensus {
#  String sample
#
#  File referenceGenome #fasta
#  File inBam
#  # TODO: input settings must compute chr_names
#  Array[String] chrNames # sampleName-1, sampleName-2, ...
#
#  Int? threads
#  Int? minCoverage
#
#  String? novoalignOptions
#
#  command {
#    assembly.py refine_assembly \
#      "${referenceGenome}" \
#      "${inBam}" \
#      "${sample}.fasta" \
#      --outBam "${sample}.realigned.bam" \
#      --outVcf "${sample}.vcf.gz" \
#      "${'--min_coverage' + minCoverage || '3'}" \
#      "${'--novo_params' + novoalignOptions || '-r Random -l 40 -g 40 -x 20 -t 100'}" \
#      --keep_all_reads \
#      --chr_names "${sep=' ' chrNames+}" \
#      "${'--threads' + threads}"
#  }
#
#  output {
#    File assembly = "${sample}.fasta"
#    File realignedBam = "${sample}.realigned.bam"
#    File variantCalls = "${sample}.vcf.gz"
#  }
#  runtime {
#    memory: "4 GB"
#    docker: "quay.io/broadinstitute/viral-ngs"
#  }
#}

task ref_guided_consensus_aligned_with_dups {
  String sample

  File realignedBam

  command {
    samtools view -b -F 4 -o "${sample}.realigned.only_aligned.bam" "${realignedBam}"
  }

  output {
    File onlyAlignedReadsBam = "${sample}.realigned.only_aligned.bam"
  }
  runtime {
    memory: "8 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}

#task ref_guided_diversity {
#  String sample
#
#  File assembly # .fasta
#  File variantCalls # .vcf.gz
#  File referenceGenome # .fasta, but also .fasta.fai and .dict
#
#  command {
#    GenomeAnalysisTK.jar
#      -T CombineVariants -R "${referenceGenome}" {inFilesString} -o {outFile}
#      --genotypemergeoption UNIQUIFY
#  }
#
#  output {
#
#  }
#  runtime {
#    docker: "quay.io/broadinstitute/viral-ngs"
#  }
#}

task multi_align_mafft_ref {
  File           reference_fasta
  Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
  String         out_prefix = basename(reference_fasta, '.fasta')
  Int?           mafft_maxIters
  Float?         mafft_ep
  Float?         mafft_gapOpeningPenalty

  command {
    interhost.py multichr_mafft \
      ${reference_fasta} ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep' + mafft_ep} \
      ${'--gapOpeningPenalty' + mafft_gapOpeningPenalty} \
      ${'--maxiters' + mafft_maxIters} \
      --outFilePrefix ${out_prefix} \
      --preservecase \
      --localpair \
      --sampleNameListFile ${out_prefix}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    File         sampleNamesFile   = "${out_prefix}-sample_names.txt"
    Array[File]+ alignments_by_chr = glob("${out_prefix}*.fasta")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "7 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x8"
  }
}

task multi_align_mafft {
  Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
  String         out_prefix
  Int?           mafft_maxIters
  Float?         mafft_ep
  Float?         mafft_gapOpeningPenalty

  command {
    interhost.py multichr_mafft \
      ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep' + mafft_ep} \
      ${'--gapOpeningPenalty' + mafft_gapOpeningPenalty} \
      ${'--maxiters' + mafft_maxIters} \
      --outFilePrefix ${out_prefix} \
      --preservecase \
      --localpair \
      --sampleNameListFile ${out_prefix}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    File        sampleNamesFile   = "${out_prefix}-sample_names.txt"
    Array[File] alignments_by_chr = glob("${out_prefix}*.fasta")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "7 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x8"
  }
}

task index_ref {
  File referenceGenome

  command {
    read_utils.py novoindex "${referenceGenome}"
    read_utils.py index_fasta_samtools "${referenceGenome}"
    read_utils.py index_fasta_picard "${referenceGenome}"
  }

  output {
    File referenceNix = "*.nix"
    File referenceFai = "*.fasta.fai"
    File referenceDict = "*.dict"
  }
  runtime {
    memory: "4 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}
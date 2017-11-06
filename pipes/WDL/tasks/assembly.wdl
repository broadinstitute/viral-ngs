# ======================================================================
# assemble_trinity: 
#   Run the Trinity assembler.
#   First trim reads with trimmomatic, rmdup with prinseq, and
#   random subsample to no more than 100k reads
# ======================================================================

task assemble_trinity {
  String sample

  File inputBam
  File trim_clip_db #fasta

  Int? trinity_n_reads
  Int? threads

  command {
    assembly.py assemble_trinity \
      "${inputBam}" \
      "${trim_clip_db}" \
      "${sample}.assembly1-trinity.fasta" \
      "${'--n_reads' + trinity_n_reads}" \
      --outReads="${sample}.subsamp.bam" \
      "${'--threads' + threads}"
  }

  output {
    File trinityAssembly = "${sample}.assembly1-trinity.fasta"
    File subsampBam = "${sample}.subsamp.bam"
  }

  runtime {
    memory: "4GB"
    docker: "broadinstitute/viral-ngs"
  }

}

task orient_and_impute {
  String sample

  File inputAssembly # fasta
  File referenceGenome # fasta

  Float minUnambig
  Float minLengthFraction # min fraction of reference length to be considered acceptible 
  Int? replaceLength
  Int? assembly_nucmer_max_gap
  Int? assembly_nucmer_min_match
  Int? assembly_nucmer_min_cluster
  Int? scaffold_min_pct_contig_aligned

  command {
    assembly.py order_and_orient \
      "${inputAssembly}" \
      "${referenceGenome}" \
      "${sample}.assembly2-scaffolded.fasta" \
      "${'--maxgap=' + assembly_nucmer_max_gap}"  \
      "${'--minmatch=' + assembly_nucmer_min_match}" \
      "${'--mincluster=' + assembly_nucmer_min_cluster}" \
      "${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned}" \
      --alternate_fasta "${sample}.assembly2-alternate_sequences.fasta"

    assembly.py impute_from_reference \
      "${sample}.assembly2-scaffolded.fasta" \
      "${referenceGenome}" \
      "${sample}.assembly3-modify.fasta" \
      --newName "${sample}" \
      "${'--replaceLength' + replaceLength}" \
      "${'--minLengthFraction' + minLengthFraction}" \
      "${'--minUnambig' + minUnambig}" \
      --index
  }

  output {
    File refinedAssembly = "${sample}.assembly3-modify.fasta"
    File scaffoldedFasta = "${sample}.assembly2-scaffolded.fasta"
    File alternate_fasta = "${sample}.assembly2-alternate_sequences.fasta"
  }

  runtime {
    memory: "12GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task refine_assembly_1 {
  String sample

  File inputAssembly # fasta
  File cleanedBam

  Int? threads
  Int? minCoverage

  String? novoalignOptions

  command {
    assembly.py refine_assembly \
      "${inputAssembly}" \
      "${sample}.assembly4-refined.fasta" \
      --outVcf "${sample}.assembly3.vcf.gz" \
      "${'--min_coverage' + minCoverage || '2'}" \
      "${'--novo_params' + novoalignOptions || '-r Random -l 30 -g 40 -x 20 -t 502'}"  \ # default notation may be wrong...
      "${'--threads' + threads}"
  }

  output {
    File assembly = "${sample}.assembly4-refined.fasta"
    File variantCalls = "${sample}.assembly3.vcf.gz"
  }

  runtime {
    memory: "7GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task refine_assembly_2 {
  String sample

  File inputAssembly #fasta
  File cleanedBam

  Int? threads
  Int? minCoverage

  String? novoalignOptions

  command {
    assembly.py refine_assembly \
      "${inputAssembly}" \
      "${sample}.fasta" \
      --outVcf "${sample}.assembly4.vcf.gz" \
      "${'--min_coverage' + minCoverage || '3'}" \
      "${'--novo_params' + novoalignOptions || '-r Random -l 40 -g 40 -x 20 -t 100'}" \
      "${'--threads' + threads}"
  }

  output {
    File assembly = "${sample}.fasta"
    File variantCalls = "${sample}.assembly4.vcf.gz"
  }

  runtime {
    memory: "7GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task map_reads_to_self {
  String sample

  File inputAssembly
  File cleanedBam

  Int? threads
  String? aligner
  String? alignerOptions

  command {
    read_utils.py align_and_fix \
      "${inputAssembly}" \
      --outBamAll "${sample}.bam" \
      --outBamFiltered "${sample}.mapped.bam" \
      "${'--aligner' + aligner || 'novoalign'}" \
      # set aligner options. The default, if novoalign, should be '-r Random -l 40 -g 40 -x 20 -t 100 -k'
      "${'--aligner_options' + alignerOptions}" \
      "${'--threads' + threads}" \
  }

  output {
    File outBam = "${sample}.bam"
    File outBamFiltered = "${sample}.mapped.bam"
  }

  runtime {
    memory: "4GB"
    docker: "broadinstitute/viral-ngs"
  }
}




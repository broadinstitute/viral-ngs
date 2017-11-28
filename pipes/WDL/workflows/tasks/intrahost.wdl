task isnvs_per_sample {
  String sample

  File mappedBam
  File assembly # fasta

  Int? threads
  Int? minReadsPerStrand
  Int? maxBias

  command {
    intrahost.py vphaser_one_sample \
        "${mappedBam}" \
        "${assembly}" \
         "vphaser2.${sample}.txt.gz" \
         "${'--vphaserNumThreads' + threads}" \
         --removeDoublyMappedReads \
         "${'--minReadsEach' + minReadsPerStrand}" \
         "${'--maxBias' + maxBias}"
  }

  output {
    File isnvsFile = "vphaser2.${sample}.txt.gz"
  }
  runtime {
    memory: "7 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}

task isnvs_vcf {
  Array[File] vphaser2Calls # vphaser output; ex. vphaser2.${sample}.txt.gz
  Array[File] perSegmentMultiAlignments # aligned_##.fasta, where ## is segment number
  File referenceGenome #fasta
  File sampleNameList

  Array[String] snpEffRef # list of accessions to build/find snpEff database
  Array[String] sampleNames # list of sample names
  String emailAddress # email address passed to NCBI if we need to download reference sequences

  command {
    intrahost.py merge_to_vcf  \
        "${referenceGenome}" \
        "isnvs.vcf.gz" \
        --samples "${sampleNameList}" \
        --isnvs "${sep=' ' vphaser2Calls}" \
        --alignments "${sep=' ' perSegmentMultiAlignments}" \
        --strip_chr_version \
        --parse_accession

    interhost.py snpEff \
        "isnvs.vcf.gz" \
        "${sep=' ' snpEffRef}" \
        "isnvs.annot.vcf.gz" \
        "${emailAddress}"

    intrahost.py iSNV_table \
        "isnvs.annot.vcf.gz" \
        "isnvs.annot.txt.gz"        
  }

  output {
    Array[File] isnvFiles = ["isnvs.vcf.gz", "isnvs.vcf.gz.tbi", "isnvs.annot.vcf.gz", "isnvs.annot.txt.gz", "isnvs.annot.vcf.gz.tbi"]
  }
  runtime {
    memory: "4 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}

task isnvs_vcf_filtered {
  Array[File] vphaser2Calls # vphaser output; ex. vphaser2.${sample}.txt.gz
  Array[File] perSegmentMultiAlignments # aligned_##.fasta, where ## is segment number
  File referenceGenome #fasta
  File sampleNameList

  Array[String] snpEffRef # list of accessions to build/find snpEff database
  Array[String] sampleNames # list of sample names
  String emailAddress # email address passed to NCBI if we need to download reference sequences
  Boolean naiveFilter

  command {
    intrahost.py merge_to_vcf \
        "${referenceGenome}" \
        "isnvs.vcf.gz" \
        --samples "${sampleNameList}" \
        --isnvs "${sep=' ' vphaser2Calls}" \
        --alignments "${sep=' ' perSegmentMultiAlignments}" \
        --strip_chr_version \
        "${'--naive_filter' + naiveFilter}" \
        --parse_accession

    interhost.py snpEff \
        "isnvs.vcf.gz" \
        "${sep=' ' snpEffRef}" \
        "isnvs.annot.vcf.gz" \
        "${emailAddress}"

    intrahost.py iSNV_table \
        "isnvs.annot.vcf.gz" \
        "isnvs.annot.txt.gz" \
  }

  output {
    Array[File] isnvFiles = ["isnvs.vcf.gz", "isnvs.vcf.gz.tbi", "isnvs.annot.vcf.gz", "isnvs.annot.txt.gz", "isnvs.annot.vcf.gz.tbi"]
  }
  runtime {
    memory: "4 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}
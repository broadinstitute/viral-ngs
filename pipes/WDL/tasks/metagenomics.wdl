task isnvs_vcf {
  Array[File] vphaser2Calls # vphaser output; ex. vphaser2.${sample}.txt.gz
  Array[File] perSegmentMultiAlignments # aligned_##.fasta, where ## is segment number
  File referenceGenome #fasta
  File sampleNameList

  Array[String] snpEffRef # list of accessions to build/find snpEff database
  Array[String] sampleNames # list of sample names
  String emailAddress # email address passed to NCBI if we need to download reference sequences

  command {
    intrahost.py merge_to_vcf \
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
    Array[File] = ["isnvs.vcf.gz", "isnvs.vcf.gz.tbi", "isnvs.annot.vcf.gz", "isnvs.annot.txt.gz", "isnvs.annot.vcf.gz.tbi"]
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# old ===============================

task diamond {
  String sample

  File inputBam

  command {
    metagenomics.py diamond
  }

  output {
    File report = "${sample}.diamond.report"
    File lca    = "${sample}.diamond.lca.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}

task kraken {
  String sample

  File inputBam

  command {
    metagenomics.py kraken
  }

  output {
    File report = "${sample}.kraken.report"
    File reads  = "${sample}.kraken.reads.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}

task align_rna {
  String sample

  File inputBam

  command {
    metagenomics.py align_rna
  }

  output {
    File report      = "{sample}.{adjective,raw|cleaned}.rna_bwa.report"
    File dupe_report = "{sample}.{adjective,raw|cleaned}.rna_bwa.dupes.report"
    File bam         = "{sample}.{adjective,raw|cleaned}.rna_bwa.bam"
    File lca         = "{sample}.{adjective,raw|cleaned}.rna_bwa.lca.gz"
    File dupes_lca   = "{sample}.{adjective,raw|cleaned}.rna_bwa.lca_dupes.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}

task krona {
  String sample

  # input...

  command {
    metagenomics.py krona 
      
  }

  output {

  }

  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}
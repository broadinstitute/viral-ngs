
task downsample_bams {
  Array[File]  reads_bam
  Int?         readCount
  Boolean?     deduplicateBefore=false
  Boolean?     deduplicateAfter=false

  command {
    if [[ "${deduplicateBefore}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateBefore"
    elif [[ "${deduplicateAfter}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateAfter"
    fi

    if [[ "${deduplicateBefore}" == "true" && "${deduplicateAfter}" == "true" ]]; then
      echo "deduplicateBefore and deduplicateAfter are mutually exclusive. Only one can be used."
      exit 1
    fi
    
    read_utils.py downsample_bams \
        ${sep=' ' reads_bam} \
        --outPath ./output \
        ${'--readCount=' + readCount} \
        $DEDUP_OPTION \
        --JVMmemory "1g"
  }

  output {
    Array[File] downsampled_bam  = glob("output/*.downsampled-*.bam")
    String      viralngs_version = "viral-ngs_version_unknown"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}

task fastq_to_bam {
  File  in_fastq1
  File  in_fastq2
  File? header

  # also used for outfile name
  String? sample_name = basename(basename(basename(basename(basename(basename(in_fastq1, ".1.fastq"),".1.fq"),".fastq"),".fq"),".fq1"),".fastq1")

  command {
    read_utils.py fastq_to_bam \
        ${in_fastq1} \
        ${in_fastq2} \
        ${sample_name}.bam \
        ${'--sampleName=' + sample_name} \
        ${'--header ' + header} \
        --JVMmemory "1g"
  }

  output {
    File out_bam            = "${sample_name}.bam"
    String viralngs_version = "viral-ngs_version_unknown"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x4"
    preemptible: 0
  }
}

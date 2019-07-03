
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
    dx_instance_type: "mem1_ssd1_x4"
  }
}

task dedup_bam {
  File  in_bam
  Int?  max_mismatches=3

  String sample_name = basename(in_bam, ".bam")

  command {
    read_utils.py rmdup_clumpify_bam \
        in_bam \
        ${sample_name}.dedup.bam \
        ${'--maxMismatches=' + max_mismatches} \
        --JVMmemory "8g"

    reports.py fastqc ${sample_name}.dedup.bam ${sample_name}.deduped_fastqc.html
  }

  output {
    File dedup_bam               = "${sample_name}.dedup.bam"
    File dedup_only_reads_fastqc = "${sample_name}.deduped_fastqc.html"
    String      viralngs_version        = "viral-ngs_version_unknown"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "52 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x32"
    preemptible: 0
  }
}
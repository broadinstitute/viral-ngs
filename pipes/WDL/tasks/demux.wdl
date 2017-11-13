
task illumina_demux {

  File flowcell_tgz # pipeable
  Int lane
  File? samplesheet
  String? sequencingCenter

  String? flowcell
  Int? minimumBaseQuality
  Int? maxMismatches = 1
  Int? minMismatchDelta
  Int? maxNoCalls
  String? readStructure
  Int? minimumQuality = 10
  String? runStartDate

  command {
    set -ex -o pipefail

    cat ${flowcell_tgz} |
      read_utils.py extract_tarball \
        - /mnt/tmp/flowcell \
        --pipe_hint=${flowcell_tgz} \
        --loglevel=DEBUG

    illumina.py illumina_demux \
      /mnt/tmp/flowcell \
      ${lane} \
      . \
      ${'--sampleSheet=' + samplesheet} \
      ${'--sequencing_center=' + sequencingCenter} \
      --outMetrics="metrics.txt" \
      --commonBarcodes="barcodes.txt" \
      ${'--flowcell=' + flowcell} \
      ${'--minimum_base_quality=' + minimumBaseQuality} \
      ${'--max_mismatches=' + maxMismatches} \
      ${'--min_mismatch_delta=' + minMismatchDelta} \
      ${'--max_no_calls=' + maxNoCalls} \
      ${'--read_structure=' + readStructure} \
      ${'--minimum_quality=' + minimumQuality} \
      ${'--run_start_date=' + runStartDate} \
      --JVMmemory=7g \
      --threads=32 \ # yes, we are overloading this on purpose, seems to speed things up
      --compression_level=5 \
      --loglevel=DEBUG \
      --tmp_dir=/mnt/tmp

    rm -f Unmatched.bam
  }

  output {
    File metrics = "metrics.txt"
    File commonBarcodes = "barcodes.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "8GB"
    cpu: 16
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
  }
}
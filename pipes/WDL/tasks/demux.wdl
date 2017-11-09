
task illumina_demux {

  File flowcell_tgz
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
    illumina.py illumina_demux \
      ${flowcell_tgz} \
      ${lane} \
      /mnt/output \
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
      --JVMmemory=15g \
      --threads=32 \
      --compression_level=5 \
      --max_records_in_ram=1000000 \
      --loglevel=DEBUG \
      --tmp_dir=/mnt/tmp
  }

  output {
    File metrics = "metrics.txt"
    File commonBarcodes = "barcodes.txt"
    Array[File] raw_reads_unaligned_bams = glob("/mnt/output/*.bam")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    #memory: "52GB"
    #cpu: 8
    memory: "14GB"
    cpu: 16
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL, /mnt/output 375 LOCAL"
  }
}
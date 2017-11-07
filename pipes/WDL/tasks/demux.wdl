
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
    set -e -o pipefail
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
    --tmp_dir=/mnt/tmp
  }

  output {
    File metrics = "metrics.txt"
    File commonBarcodes = "barcodes.txt"
    Array[File] raw_reads_unaligned_bams = glob("/mnt/output/*.bam")
  }

  runtime {
    docker: "broadinstitute/viral-ngs:latest"
    memory: "52GB"
    cpu: "8"
    preemptible: 0
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL, /mnt/output 375 LOCAL"
  }
}
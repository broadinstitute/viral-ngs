task illumina_demux {

  File inDir
  Int lane
  File samplesheet
  String sequencingCenter

  String flowcell

  Int? minimumBaseQuality
  Int? maxMismatches
  Int? minMismatchDelta
  Int? maxNoCalls
  String? readStructure
  Int? minimumQuality
  String? runStartDate

  command {
    illumina.py illumina_demux \
    "${inDir}" \
    "${lane}" \
    "${flowcell}.${lane}" \

    --sampleSheet="${samplesheet}" \
    --sequencing_center="${sequencingCenter}" \
    --outMetrics="metrics.txt" \
    --commonBarcodes="barcodes.txt" \
    "${'--flowcell' + flowcell}" \

    "${'--minimum_base_quality' + minimumBaseQuality}" \
    "${'--max_mismatches' + maxMismatches}" \
    "${'--min_mismatch_delta' + minMismatchDelta}" \
    "${'--max_no_calls' + maxNoCalls}" \
    "${'--read_structure' + readStructure}" \
    "${'--minimum_quality' + minimumQuality}" \
    "${'--run_start_date' + runStartDate}"
  }

  output {
    File metrics = "metrics.txt"
    File commonBarcodes = "barcodes.txt"
    Array[File] outputBams = glob("${flowcell}.${lane}/*.bam")
  }
  runtime {
    memory: "60GB"
    docker: "broadinstitute/viral-ngs"
  }
}
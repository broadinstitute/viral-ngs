
task illumina_demux {

  File    flowcell_tgz
  Int?    lane=1
  File?   samplesheet
  String? sequencingCenter

  String? flowcell
  Int?    minimumBaseQuality
  Int?    maxMismatches = 1
  Int?    minMismatchDelta
  Int?    maxNoCalls
  String? readStructure
  Int?    minimumQuality = 10
  String? runStartDate

  parameter_meta {
    flowcell_tgz : "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v illumina.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi
    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    FLOWCELL_DIR=$(mktemp -d)

    cat ${flowcell_tgz} |
      read_utils.py extract_tarball \
        - $FLOWCELL_DIR \
        --pipe_hint=${flowcell_tgz} \
        --loglevel=DEBUG

    # note that we are intentionally setting --threads to about 2x the core
    # count. seems to still provide speed benefit (over 1x) when doing so.
    illumina.py illumina_demux \
      $FLOWCELL_DIR \
      ${lane} \
      . \
      ${'--sampleSheet=' + samplesheet} \
      ${'--sequencing_center=' + sequencingCenter} \
      --outMetrics=metrics.txt \
      --commonBarcodes=barcodes.txt \
      ${'--flowcell=' + flowcell} \
      ${'--minimum_base_quality=' + minimumBaseQuality} \
      ${'--max_mismatches=' + maxMismatches} \
      ${'--min_mismatch_delta=' + minMismatchDelta} \
      ${'--max_no_calls=' + maxNoCalls} \
      ${'--read_structure=' + readStructure} \
      ${'--minimum_quality=' + minimumQuality} \
      ${'--run_start_date=' + runStartDate} \
      --JVMmemory=14g \
      --threads=64 \
      --compression_level=5 \
      --loglevel=DEBUG

    rm -f Unmatched.bam
  }

  output {
    File        metrics                  = "metrics.txt"
    File        commonBarcodes           = "barcodes.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "16 GB"
    cpu: 32
    dx_instance_type: "mem1_ssd2_x36"
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
  }
}

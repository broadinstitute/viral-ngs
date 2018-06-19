
#task merge_tar_chunks {
#  Array[File]+  tar_chunks
#  String        out_filename
#  command {
#    set -ex -o pipefail
#
#    if [ -d /mnt/tmp ]; then
#      TMPDIR=/mnt/tmp
#    fi
#    FLOWCELL_DIR=$(mktemp -d)
#
#    read_utils.py extract_tarball \
#      ${flowcell_tgz} $FLOWCELL_DIR \
#      --loglevel=DEBUG
#
#  }
#  output {
#    File tar_lz4=
#  }
#  runtime {
#    docker: "quay.io/broadinstitute/viral-ngs"
#    memory: "7 GB"
#    cpu: 4
#    dx_instance_type: "mem1_hdd2_x8"
#    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
#  }
#}

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

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    FLOWCELL_DIR=$(mktemp -d)

    read_utils.py extract_tarball \
      ${flowcell_tgz} $FLOWCELL_DIR \
      --loglevel=DEBUG

    # find 95% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_95.sh`

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
      --JVMmemory="$mem_in_mb"m \
      --threads=64 \
      --compression_level=5 \
      --loglevel=DEBUG

    rm -f Unmatched.bam
    for bam in *.bam; do
      fastqc_out=$(basename $bam .bam)_fastqc.html
      reports.py fastqc $bam $fastqc_out
    done
  }

  output {
    File        metrics                  = "metrics.txt"
    File        commonBarcodes           = "barcodes.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
    Array[File] raw_reads_fastqc         = glob("*_fastqc.html")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "16 GB"
    cpu: 32
    dx_instance_type: "mem1_ssd2_x36"
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
  }
}

task merge_and_reheader_bams {
  Array[File]+  in_bams
  File?         reheader_table # tsv with 3 cols: field, old value, new value
  String        out_basename

  command {
    set -ex -o pipefail

    if [ ${length(in_bams)} -gt 1 ]; then
      read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --loglevel DEBUG
    else
      echo "Skipping merge, only one input file"
      ln -s ${select_first(in_bams)} merged.bam
    fi    

    if [[ -f "${reheader_table}" ]]; then
      read_utils.py reheader_bam merged.bam ${reheader_table} ${out_basename}.bam --loglevel DEBUG
    else
      echo "Skipping reheader, no mapping table specified"
      ln -s merged.bam ${out_basename}.bam
    fi
  }

  output {
    File  out_bam = "${out_basename}.bam"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2000 MB"
    cpu: 2
    dx_instance_type: "mem1_ssd2_x4"
  }
}

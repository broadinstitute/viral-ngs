# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete_taxa {
  File         raw_reads_unmapped_bam
  Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Int?         query_chunk_size

  String      bam_basename = basename(raw_reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi

    # find 50% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_50.sh`

    # bmtagger and blast db args
    DBS_BMTAGGER="${sep=' ' bmtaggerDbs}"
    DBS_BLAST="${sep=' ' blastDbs}"
    if [ -n "$DBS_BMTAGGER" ]; then DBS_BMTAGGER="--bmtaggerDbs $DBS_BMTAGGER"; fi
    if [ -n "$DBS_BLAST" ]; then DBS_BLAST="--blastDbs $DBS_BLAST"; fi

    # run depletion
    taxon_filter.py deplete_human \
      ${raw_reads_unmapped_bam} \
      tmpfile.raw.bam \
      tmpfile.bmtagger_depleted.bam \
      tmpfile.rmdup.bam \
      ${bam_basename}.cleaned.bam \
      $DBS_BMTAGGER $DBS_BLAST \
      ${'--chunkSize=' + query_chunk_size} \
      --JVMmemory="$mem_in_mb"m \
      --srprismMemory=$mem_in_mb \
      --loglevel=DEBUG

    samtools view -c ${raw_reads_unmapped_bam} | tee depletion_read_count_pre
    samtools view -c ${bam_basename}.cleaned.bam | tee depletion_read_count_post
    reports.py fastqc ${bam_basename}.cleaned.bam ${bam_basename}.cleaned_fastqc.html
  }

  output {
    File cleaned_bam               = "${bam_basename}.cleaned.bam"
    File cleaned_fastqc            = "${bam_basename}.cleaned_fastqc.html"
    Int  depletion_read_count_pre  = read_int("depletion_read_count_pre")
    Int  depletion_read_count_post = read_int("depletion_read_count_post")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "14 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x8"
    preemptible: 0
  }
}


# ======================================================================
# filter_to_taxon: 
#   This step reduces the read set to a specific taxon (usually the genus
#   level or greater for the virus of interest)
# ======================================================================
task filter_to_taxon {
  File reads_unmapped_bam
  File lastal_db_fasta

  # do this in two steps in case the input doesn't actually have "cleaned" in the name
  String bam_basename = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_90.sh`

    taxon_filter.py filter_lastal_bam \
      ${reads_unmapped_bam} \
      ${lastal_db_fasta} \
      ${bam_basename}.taxfilt.bam \
      --JVMmemory="$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c ${bam_basename}.taxfilt.bam | tee filter_read_count_post
    reports.py fastqc ${bam_basename}.taxfilt.bam ${bam_basename}.taxfilt_fastqc.html
  }

  output {
    File taxfilt_bam            = "${bam_basename}.taxfilt.bam"
    File taxfilt_fastqc         = "${bam_basename}.taxfilt_fastqc.html"
    Int  filter_read_count_post = read_int("filter_read_count_post")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "14 GB"
    cpu: 16
    dx_instance_type: "mem1_ssd1_x8"
  }
}


task merge_one_per_sample {
  String       out_bam_basename
  Array[File]+ inputBams
  Boolean?     rmdup=false

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_90.sh`

    read_utils.py merge_bams \
      "${sep=' ' inputBams}" \
      "${out_bam_basename}.bam" \
      --picardOptions SORT_ORDER=queryname \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    if [[ "${rmdup}" == "true" ]]; then
      mv "${out_bam_basename}.bam" tmp.bam
      read_utils.py rmdup_mvicuna_bam \
        tmp.bam \
        ${out_bam_basename}.bam \
        --JVMmemory "$mem_in_mb"m \
        --loglevel=DEBUG
    fi
  }

  output {
    File mergedBam = "${out_bam_basename}.bam"
  }

  runtime{
    memory: "7 GB"
    cpu: 4
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}

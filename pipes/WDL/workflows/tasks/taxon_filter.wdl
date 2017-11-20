# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete_taxa {
  File        raw_reads_unmapped_bam
  Array[File] bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Array[File] blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Int?        query_chunk_size=0

  String      bam_basename = basename(raw_reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail

    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v taxon_filter.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi
    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi

    taxon_filter.py deplete_human \
      ${raw_reads_unmapped_bam} \
      tmpfile.raw.bam \
      tmpfile.bmtagger_depleted.bam \
      tmpfile.rmdup.bam \
      ${bam_basename}.cleaned.bam \
      --bmtaggerDbs ${sep=' ' bmtaggerDbs} \
      --blastDbs ${sep=' ' blastDbs} \
      --chunkSize ${query_chunk_size} \
      --JVMmemory=14g \
      --loglevel=DEBUG

    samtools view -c ${raw_reads_unmapped_bam} | tee depletion_read_count_pre
    samtools view -c ${bam_basename}.cleaned.bam | tee depletion_read_count_post
  }

  output {
    File cleaned_bam               = "${bam_basename}.cleaned.bam"
    Int  depletion_read_count_pre  = read_int("depletion_read_count_pre")
    Int  depletion_read_count_post = read_int("depletion_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14 GB"
    cpu: 8
    preemptible: 1
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

    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v taxon_filter.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi

    taxon_filter.py filter_lastal_bam \
      ${reads_unmapped_bam} \
      ${lastal_db_fasta} \
      ${bam_basename}.taxfilt.bam \
      --JVMmemory=14g \
      --loglevel=DEBUG

    samtools view -c ${bam_basename}.taxfilt.bam | tee filter_read_count_post
  }

  output {
    File taxfilt_bam            = "${bam_basename}.taxfilt.bam"
    Int  filter_read_count_post = read_int("filter_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14 GB"
    cpu: 16
  }
}


task merge_one_per_sample {
  String       out_bam_basename
  Array[File]+ inputBams
  Boolean?     rmdup=false

  command {
    set -ex -o pipefail

    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v read_utils.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi

    read_utils.py merge_bams \
      "${sep=' ' inputBams}" \
      "${out_bam_basename}.bam" \
      --picardOptions SORT_ORDER=queryname \
      --JVMmemory 7g \
      --loglevel=DEBUG

    if [[ "${rmdup}" == "true" ]]; then
      mv "${out_bam_basename}.bam" tmp.bam
      read_utils.py rmdup_mvicuna_bam \
        tmp.bam \
        ${out_bam_basename}.bam \
        --JVMmemory 7g \
        --loglevel=DEBUG
    fi
  }

  output {
    File mergedBam = "${out_bam_basename}.bam"
  }

  runtime{
    memory: "7 GB"
    cpu: 4
    docker: "broadinstitute/viral-ngs"
  }
}

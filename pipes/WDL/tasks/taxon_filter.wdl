# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete_taxa {
  File raw_reads_unmapped_bam
  Array[File] bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Array[File] blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz

  command {
    set -ex -o pipefail

    # set inputs/outputs
    echo "$(basename ${raw_reads_unmapped_bam} .bam).cleaned.bam" > fname-out-cleaned.txt

    taxon_filter.py deplete_human \
      ${raw_reads_unmapped_bam} \
      /mnt/tmp/tmpfile.raw.bam \
      /mnt/tmp/tmpfile.bmtagger_depleted.bam \
      /mnt/tmp/tmpfile.rmdup.bam \
      `cat fname-out-cleaned.txt` \
      --bmtaggerDbs `cat ${write_lines(bmtaggerDbs)}` \
      --blastDbs `cat ${write_lines(blastDbs)}` \
      --chunkSize=0 \
      --JVMmemory=14g \
      --tmp_dir=/mnt/tmp \
      --loglevel=DEBUG

    samtools view -c "${raw_reads_unmapped_bam}" | tee depletion_read_count_pre
    samtools view -c `cat fname-out-cleaned.txt` | tee depletion_read_count_post
  }

  output {
    File cleaned_bam               = read_string("fname-out-cleaned.txt")
    Int  depletion_read_count_pre  = read_int("depletion_read_count_pre")
    Int  depletion_read_count_post = read_int("depletion_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14GB"
    cpu: 8
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
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

  command {
    set -ex -o pipefail

    # set inputs/outputs
    # do this in two steps in case the input doesn't actually have "cleaned" in the name
    BASE_NAME="$(basename ${reads_unmapped_bam} .bam)"
    echo "$(basename $BASE_NAME .cleaned).taxfilt.bam" > fname-out-taxfilt.txt

    taxon_filter.py filter_lastal_bam \
      ${reads_unmapped_bam} \
      ${lastal_db_fasta} \
      `cat fname-out-taxfilt.txt` \
      --JVMmemory=14g \
      --tmp_dir=/mnt/tmp \
      --loglevel=DEBUG

    samtools view -c `cat fname-out-taxfilt.txt` | tee filter_read_count_post
  }

  output {
    File taxfilt_bam            = read_string("fname-out-taxfilt.txt")
    Int  filter_read_count_post = read_int("filter_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14GB"
    cpu: 16
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
  }
}


# ======================================================================
# merge_one_per_sample:
#   All of the above depletion steps are run once per flowcell per lane per
#   multiplexed sample.  This reduces recomputation on the same data when
#   additional sequencing runs are performed on the same samples.  This
#   step merges reads down to one-per-sample, which all downstream
#   analysis steps require.  For cleaned and taxfilt outputs, we re-run
#   the rmdup step on a per-library basis after merging.
# ======================================================================

task merge_one_per_sample {
  String out_bam_name
  Array[File] inputBams

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      ${out_bam_name}.bam \
      --picardOptions SORT_ORDER=queryname \
      --JVMmemory 7g \
      --tmp_dir=/mnt/tmp \
      --loglevel=DEBUG
  }

  output {
    File mergedBam = "${out_bam_name}.bam"
  }

  runtime{
    memory: "7GB"
    cpu: 4
    docker: "broadinstitute/viral-ngs"
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
  }
}

task merge_one_per_sample_rmdup {
  String out_bam_name
  Array[File] inputBams

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      /mnt/tmp/temp_merged-${out_bam_name}.bam \
      --picardOptions SORT_ORDER=queryname \
      --JVMmemory 7g \
      --tmp_dir=/mnt/tmp \
      --loglevel=DEBUG

    read_utils.py rmdup_mvicuna_bam \
      /mnt/tmp/temp_merged-${out_bam_name}.bam \
      ${out_bam_name}.bam \
      --JVMmemory 7g \
      --tmp_dir=/mnt/tmp \
      --loglevel=DEBUG
  }

  output {
    File mergedBam = "${out_bam_name}.bam"
  }

  runtime{
    memory: "7GB"
    cpu: 4
    docker: "broadinstitute/viral-ngs"
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
  }
}
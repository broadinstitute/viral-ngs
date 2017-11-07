# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete {
  String sample

  File raw_reads_unmapped_bam
  Array[File] bmtaggerDbs
  Array[File] blastDbs

  command {
    set -ex -o pipefail

    # stage the databases for BMTagger and BLAST
    # assumptions: each database is stored in a tarball. If the database name
    # is X then the tarball is named X.bmtagger_db.tar.gz or X.blastndb.tar.gz.
    # The tarball contains the database files in the root (NOT in subdirectory
    # X/). The individual database files have X as their base name, e.g.
    # X.srprism.amp, X.nin
    stage_db() {
        local dbname=$(basename "$2" .$1.tar.gz)
        mkdir -p "$TMP_DIR/$1/$dbname"
        cat "$2" | gzip -dc | tar xv -C "$TMP_DIR/$1/$dbname"
        rm "$2"
    }
    export -f stage_db

    taxon_filter.py deplete_human \
      ${inputBam} \
      /mnt/tmp/tmpfile-${sample}.raw.bam \
      /mnt/tmp/tmpfile-${sample}.bmtagger_depleted.bam \
      /mnt/tmp/tmpfile-${sample}.rmdup.bam \
      /mnt/output/${sample}.cleaned.bam \
      --bmtaggerDbs $(ls -1 $TMP_DIR/bmtagger_db | xargs -i -n 1 printf " $TMP_DIR/bmtagger_db/%s/%s " {} {}) \
      --blastDbs $(ls -1 $TMP_DIR/blastndb | xargs -i -n 1 printf " $TMP_DIR/blastndb/%s/%s " {} {})
      --chunkSize=0 \
      --JVMmemory=14g \
      --tmp_dir=/mnt/tmp

    samtools view -c "${reads_unmapped_bam}" | tee depletion_read_count_pre
    samtools view -c "${sample_name}.cleaned.bam" | tee depletion_read_count_post
  }

  output {
    File cleaned_bam          = "/mnt/output/${sample}.cleaned.bam"
    Int depletion_read_count_pre = read_int("depletion_read_count_pre")
    Int depletion_read_count_post = read_int("depletion_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14GB"
    cpu: "8"
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL, /mnt/output 375 LOCAL"
  }
}


# ======================================================================
# filter_to_taxon: 
#   This step reduces the read set to a specific taxon (usually the genus
#   level or greater for the virus of interest)
# ======================================================================
task filter_to_taxon {
  String sample_name

  File reads # unmapped bam
  File lastal_db # fasta

  command {
    set -ex -o pipefail

    taxon_filter.py filter_lastal_bam \
      "${inputBam}" \
      "${lastalDbPath}" \
      "${sample_name}.taxfilt.bam" \
      --JVMmemory=7g \
      tmp_dir=/mnt/tmp

    samtools view -c "${sample_name}.taxfilt.bam" | tee filter_read_count_post
  }

  output {
    File taxfilt_bam = "${sample_name}.taxfilt.bam"
    Int filter_read_count_post = read_int("filter_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "8"
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
      --tmp_dir=/mnt/tmp
  }

  output {
    File mergedBam = ${out_bam_name}.bam
  }

  runtime{
    memory: "7GB"
    cpu: "4"
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
      --tmp_dir=/mnt/tmp

    read_utils.py rmdup_mvicuna_bam \
      /mnt/tmp/temp_merged-${out_bam_name}.bam \
      ${sample}.${adjective}.bam \
      --JVMmemory 7g \
      --tmp_dir=/mnt/tmp
  }

  output {
    File mergedBam = ${out_bam_name}.bam
  }

  runtime{
    memory: "7GB"
    cpu: "4"
    docker: "broadinstitute/viral-ngs"
    disks: "local-disk 375 LOCAL, /mnt/tmp 375 LOCAL"
  }
}
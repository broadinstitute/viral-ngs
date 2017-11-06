# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete {
  String sample

  File inputBam
  Array[File] bmtaggerDbs
  Array[File] blastDbs

  Int? JVMmemory = 14
  
  command <<<
    set -ex -o pipefail

    TMP_DIR=/tmp

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
      "${inputBam}" \
      "tmpfile-${sample}.raw.bam" \
      "tmpfile-${sample}.bmtagger_depleted.bam" \
      "tmpfile-${sample}.rmdup.bam" \
      "${sample}.cleaned.bam" \
      --bmtaggerDbs $(ls -1 $TMP_DIR/bmtagger_db | xargs -i -n 1 printf " $TMP_DIR/bmtagger_db/%s/%s " {} {}) \
      --blastDbs $(ls -1 $TMP_DIR/blastndb | xargs -i -n 1 printf " $TMP_DIR/blastndb/%s/%s " {} {})
      --chunkSize=0 \
      --threads $(nproc) \
      "${'--JVMmemory' + JVMmemory + 'g'}"

    samtools view -c "${reads_unmapped_bam}" | tee depletion_read_count_pre
    samtools view -c "${sample_name}.cleaned.bam" | tee depletion_read_count_post
  >>>

  output {
    File cleaned_bam          = "${sample}.cleaned.bam"
    Int depletion_read_count_pre = read_int("depletion_read_count_pre")
    Int depletion_read_count_post = read_int("depletion_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "14GB"
    cpu: "8"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}


# ======================================================================
# filterToTaxon: 
#   This step reduces the read set to a specific taxon (usually the genus
#   level or greater for the virus of interest)
# ======================================================================
task filterToTaxon {
  String sample_name

  File reads
  File lastal_db_tgz

  command <<<
    set -ex -o pipefail

    taxon_filter.py filter_lastal_bam \
      "${inputBam}" \
      "${lastalDbPath}" \
      "${sample_name}.taxfilt.bam"

    samtools view -c "${sample_name}.taxfilt.bam" | tee filter_read_count_post
  >>>

  output {
    File taxfilt_bam = "${sample_name}.taxfilt.bam"
    Int filter_read_count_post = read_int("filter_read_count_post")
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "8"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
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
  String sample
  String adjective

  Array[File] inputBams

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      "${sample}.${adjective}.bam" \
      --picardOptions SORT_ORDER=queryname
  }

  output {
    File mergedBam = "${sample}.${adjective}.bam"
  }

  runtime{
    memory: "10GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task merge_one_per_sample_rmdup {
  String sample
  String adjective

  Array[File] inputBams

  Int? JVMmemory

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      "temp_merged_${sample}.bam" \
      --picardOptions SORT_ORDER=queryname;
    read_utils.py rmdup_mvicuna_bam \
      "temp_merged_${sample}.bam" \
      "${sample}.${adjective}.bam" \
      "${'--JVMmemory' + JVMmemory || '8g'}"
  }

  output {
    File mergedBam = "${sample}.${adjective}.bam"
  }

  runtime{
    memory: "10GB"
    docker: "broadinstitute/viral-ngs"
  }
}
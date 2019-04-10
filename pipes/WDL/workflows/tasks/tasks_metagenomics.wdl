task krakenuniq {
  Array[File] reads_unmapped_bam
  File        krakenuniq_db_tar_lz4  # krakenuniq/{database.kdb,taxonomy}
  File        krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

  parameter_meta {
    krakenuniq_db_tar_lz4:  "stream" # for DNAnexus, until WDL implements the File| type
    krona_taxonomy_db_tgz : "stream" # for DNAnexus, until WDL implements the File| type
    reads_unmapped_bam: "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${krakenuniq_db_tar_lz4} $DB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG &  # we don't need this until later

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    OUT_BASENAME=basenames_reports.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "$(basename $bam .bam).krakenuniq-reads.txt.gz" >> $OUT_READS
      echo "$(basename $bam .bam).krakenuniq" >> $OUT_BASENAME
      echo "$(basename $bam .bam).krakenuniq-summary_report.txt" >> $OUT_REPORTS
    done

    # execute on all inputs and outputs serially, but with a single
    # database load into ram
    metagenomics.py krakenuniq \
      $DB_DIR/krakenuniq \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract

    # run single-threaded krona on up to nproc samples at once
    parallel -I ,, \
      "metagenomics.py krona \
        ,,-summary_report.txt \
        taxonomy \
        ,,.krona.html \
        --noRank --noHits --inputType krakenuniq \
        --loglevel=DEBUG" \
      ::: `cat $OUT_BASENAME`
  }

  output {
    Array[File] krakenuniq_classified_reads = glob("*.krakenuniq-reads.txt.gz")
    Array[File] krakenuniq_summary_report   = glob("*.krakenuniq-summary_report.txt")
    Array[File] krona_report_html       = glob("*.krakenuniq.krona.html")
    String      viralngs_version        = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "200 GB"
    cpu: 32
    dx_instance_type: "mem3_ssd1_x32"
    preemptible: 0
  }
}

task krona {
  File  classified_reads_txt_gz
  File  krona_taxonomy_db_tgz

  String input_basename = basename(classified_reads_txt_gz, ".txt.gz")

  parameter_meta {
    krona_taxonomy_db_tgz : "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      taxonomy \
      ${input_basename}.html \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File krona_report_html = "${input_basename}.html"
    String viralngs_version   = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "4 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd2_x2"
  }
}

task filter_bam_to_taxa {
  File classified_bam
  File classified_reads_txt_gz
  File ncbi_taxonomy_db_tgz # nodes.dmp names.dmp
  Array[String]? taxonomic_names
  Array[Int]? taxonomic_ids
  Boolean? withoutChildren=false

  String input_basename = basename(classified_bam, ".bam")

  parameter_meta {
    ncbi_taxonomy_db_tgz              : "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    TAX_NAMES="${sep=' ' taxonomic_names}"
    if [ -n "$TAX_NAMES" ]; then TAX_NAMES="--taxNames $TAX_NAMES"; fi

    TAX_IDs="${sep=' ' taxonomic_ids}"
    if [ -n "$TAX_IDs" ]; then TAX_IDs="--taxIDs $TAX_IDs"; fi

    metagenomics.py filter_bam_to_taxa \
      ${classified_bam} \
      ${classified_reads_txt_gz} \
      "${input_basename}_filtered.bam" \
      taxonomy/nodes.dmp \
      taxonomy/names.dmp \
      $TAX_NAMES \
      $TAX_IDs \
      ${true='--without-children' false='' withoutChildren} \
      --loglevel=DEBUG

      samtools view -c ${classified_bam} | tee classified_taxonomic_filter_read_count_pre
      samtools view -c "${input_basename}_filtered.bam" | tee classified_taxonomic_filter_read_count_post
  }

  output {
    File   bam_filtered_to_taxa                        = "${input_basename}_filtered.bam"
    Int    classified_taxonomic_filter_read_count_pre  = read_int("classified_taxonomic_filter_read_count_pre")
    Int    classified_taxonomic_filter_read_count_post = read_int("classified_taxonomic_filter_read_count_post")
    String viralngs_version                            = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "4 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd2_x2"
  }

}

task kaiju {
  File  reads_unmapped_bam
  File  kaiju_db_lz4  # <something>.fmi
  File  ncbi_taxonomy_db_tgz # taxonomy/{nodes.dmp, names.dmp}
  File  krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

  String input_basename = basename(reads_unmapped_bam, ".bam")

  parameter_meta {
    kaiju_db_lz4            : "stream" # for DNAnexus, until WDL implements the File| type
    ncbi_taxonomy_db_tgz    : "stream"
    krona_taxonomy_db_tgz   : "stream"
  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    lz4 -d ${kaiju_db_lz4} $DB_DIR/kaiju.fmi

    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} $DB_DIR \
      --loglevel=DEBUG

    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG

    # classify contigs
    metagenomics.py kaiju \
      ${reads_unmapped_bam} \
      $DB_DIR/kaiju.fmi \
      $DB_DIR/taxonomy \
      ${input_basename}.kaiju.report.txt \
      --outReads ${input_basename}.kaiju.reads.gz \
      --loglevel=DEBUG

    # run krona
    metagenomics.py krona \
      ${input_basename}.kaiju.report.txt \
      taxonomy \
      ${input_basename}.kaiju.html \
      --inputType kaiju \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File kaiju_report = "${input_basename}.kaiju.report.txt"
    File kaiju_reads = "${input_basename}.kaiju.reads.gz"
    File krona_report_html = "${input_basename}.kaiju.html"
    String viralngs_version                 = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "100 GB"
    cpu: 16
    dx_instance_type: "mem3_ssd1_x16"
  }
}

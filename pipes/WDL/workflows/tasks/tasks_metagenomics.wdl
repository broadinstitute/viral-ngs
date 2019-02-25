
# TO DO:
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken {
  Array[File] reads_unmapped_bam
  File        kraken_db_tar_lz4
  File        krona_taxonomy_db_tgz

  parameter_meta {
    kraken_db_tar_lz4:  "stream" # for DNAnexus, until WDL implements the File| type
    krona_taxonomy_db_tgz : "stream" # for DNAnexus, until WDL implements the File| type
    #reads_unmapped_bam: "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${kraken_db_tar_lz4} $DB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} . \
      --loglevel=DEBUG &  # we don't need this until later

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    OUT_BASENAME=basenames_reads.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "$(basename $bam .bam).kraken-reads" >> $OUT_BASENAME
      echo "$(basename $bam .bam).kraken-reads.txt.gz" >> $OUT_READS
      echo "$(basename $bam .bam).kraken-summary_report.txt" >> $OUT_REPORTS
    done

    # execute on all inputs and outputs serially, but with a single
    # database load into ram
    metagenomics.py kraken \
      $DB_DIR \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract

    # run single-threaded krona on up to nproc samples at once
    parallel -I ,, \
      "metagenomics.py krona \
        ,,.txt.gz \
        taxonomy \
        ,,.html \
        --noRank --noHits \
        --loglevel=DEBUG" \
      ::: `cat $OUT_BASENAME`
    # run single-threaded gzip on up to nproc samples at once
    parallel -I ,, "tar czf ,,.krona.tar.gz ,,.html*" ::: `cat $OUT_BASENAME`
  }

  output {
    Array[File] kraken_classified_reads = glob("*.kraken-reads.txt.gz")
    Array[File] kraken_summary_report   = glob("*.kraken-summary_report.txt")
    Array[File] krona_report_html       = glob("*.kraken-reads.html")
    Array[File] krona_report_tgz        = glob("*.kraken-reads.krona.tar.gz")
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

    tar czf ${input_basename}.krona.tar.gz ${input_basename}.html*
  }

  output {
    File   krona_report_html  = "${input_basename}.html"
    File   krona_report_tgz   = "${input_basename}.krona.tar.gz"
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

task diamond_contigs {
  File  contigs_fasta
  File  reads_unmapped_bam
  File  diamond_db_lz4
  File  diamond_taxonomy_db_tar_lz4
  File  krona_taxonomy_db_tar_lz4

  String contigs_basename = basename(contigs_fasta, ".fasta")

  parameter_meta {
    diamond_db_lz4              : "stream" # for DNAnexus, until WDL implements the File| type
    diamond_taxonomy_db_tar_lz4 : "stream" # for DNAnexus, until WDL implements the File| type
    krona_taxonomy_db_tar_lz4   : "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DIAMOND_TAXDB_DIR=$(mktemp -d)

    # find 90% memory
    mem_in_gb=`/opt/viral-ngs/source/docker/calc_mem.py gb 90`

    # decompress DBs to /mnt/db
    cat ${diamond_db_lz4} | lz4 -d > $TMPDIR/diamond_db.dmnd &
    read_utils.py extract_tarball \
      ${diamond_taxonomy_db_tar_lz4} $DIAMOND_TAXDB_DIR \
      --loglevel=DEBUG &
    wait
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tar_lz4} . \
      --loglevel=DEBUG &  # we don't need this until later

    # classify contigs
    metagenomics.py diamond_fasta \
      ${contigs_fasta} \
      $TMPDIR/diamond_db.dmnd \
      $DIAMOND_TAXDB_DIR/taxonomy/ \
      ${contigs_basename}.diamond.fasta \
      --memLimitGb $mem_in_gb \
      --loglevel=DEBUG

    # map reads to contigs & create kraken-like read report
    bwa index ${contigs_basename}.diamond.fasta
    metagenomics.py align_rna \
      ${reads_unmapped_bam} \
      ${contigs_basename}.diamond.fasta \
      $DIAMOND_TAXDB_DIR/taxonomy/ \
      ${contigs_basename}.diamond.summary_report.txt \
      --outReads ${contigs_basename}.diamond.reads.txt.gz \
      --dupeReads ${contigs_basename}.diamond.reads_w_dupes.txt.gz \
      --outBam ${contigs_basename}.diamond.bam \
      --loglevel=DEBUG

    # run krona
    wait # for krona_taxonomy_db_tgz to download and extract
    metagenomics.py krona \
      ${contigs_basename}.diamond.reads.txt.gz \
      taxonomy \
      ${contigs_basename}.diamond.html \
      --noRank --noHits \
      --loglevel=DEBUG
    tar czf ${contigs_basename}.diamond.krona.tar.gz ${contigs_basename}.diamond.html*
  }

  output {
    File   diamond_contigs                  = "${contigs_basename}.diamond.fasta"
    File   reads_mapped_to_contigs          = "${contigs_basename}.diamond.bam"
    File   diamond_summary_report           = "${contigs_basename}.diamond.summary_report.txt"
    File   diamond_classified_reads         = "${contigs_basename}.diamond.reads.txt.gz"
    File   diamond_classified_reads_w_dupes = "${contigs_basename}.diamond.reads_w_dupes.txt.gz"
    File   krona_report_html                = "${contigs_basename}.diamond.html"
    File   krona_report_tgz                 = "${contigs_basename}.diamond.krona.tar.gz"
    String viralngs_version                 = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "100 GB"
    cpu: 16
    dx_instance_type: "mem3_ssd1_x16"
  }
}



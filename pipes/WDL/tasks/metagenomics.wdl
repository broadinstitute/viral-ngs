
# TO DO:
# kraken_multi (Array[File] bam inputs and reads/reports outputs)
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken_single {
  File reads_unmapped_bam
  File kraken_db_tar_lz4 # pipeable

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    cat ${kraken_db_tar_lz4} |
      read_utils.py extract_tarball \
        - /mnt/db \
        --pipe_hint=${kraken_db_tar_lz4} \
        --loglevel=DEBUG

    time metagenomics.py kraken \
      /mnt/db \
      ${reads_unmapped_bam} \
      --outReads=kraken-reads.txt.gz \
      --outReport=kraken-report.txt \
      --loglevel=DEBUG
  }

  output {
    File kraken_classified_reads = "kraken-reads.txt.gz"
    File kraken_summary_report = "kraken-report.txt"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "200GB"
    cpu: 32
    disks: "local-disk 375 LOCAL, /mnt/db 375 LOCAL"
  }
}

task krona {
  File classified_reads_txt_gz
  File? krona_taxonomy_db_tgz = "gs://sabeti-public-dbs/krona/krona_taxonomy_20160502.tar.lz4" # pipeable

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    cat ${krona_taxonomy_db_tgz} |
      read_utils.py extract_tarball \
        - . \
        --pipe_hint=${krona_taxonomy_db_tgz} \
        --loglevel=DEBUG

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      taxonomy \
      krona-report.html \
      --noRank --noHits \
      --loglevel=DEBUG

    tar czf krona-report.tar.gz krona-report.html*
  }

  output {
    File krona_report_html = "krona-report.html"
    File krona_report_tgz = "krona-report.tar.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "2GB"
    cpu: 1
  }
}


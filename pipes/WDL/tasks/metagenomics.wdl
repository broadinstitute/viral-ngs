
# TO DO:
# kraken_multi (Array[File] bam inputs and reads/reports outputs)
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken_single {
  File reads_unmapped_bam
  File kraken_db_tar_lz4

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball ${kraken_db_tar_lz4} /mnt/db
    date >&2

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
  File? krona_taxonomy_db_tgz = "gs://sabeti-public-dbs/krona/krona_taxonomy_20160502.tar.lz4"

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    read_utils.py extract_tarball ${krona_taxonomy_db_tgz} .

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


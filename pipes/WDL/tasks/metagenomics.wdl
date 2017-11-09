
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
    date >&2
    echo "decompressing ${kraken_db_tar_lz4}"
    decompressor="pigz -dc"
    if [[ "${kraken_db_tar_lz4}" == *.lz4 ]]; then
      decompressor="lz4 -d"
    elif [[ "${kraken_db_tar_lz4}" == *.bz2 ]]; then
      decompressor="bzcat -d"
    fi
    cat ${kraken_db_tar_lz4} | $decompressor | tar -C /mnt/db -xv >&2
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
    docker: "broadinstitute/viral-ngs-dev:dp_wdl"
    memory: "200GB"
    cpu: 32
    disks: "local-disk 375 LOCAL, /mnt/db 375 LOCAL"
  }
}

task krona {
  File classified_reads_txt_gz
  File krona_taxonomy_db_tgz

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    date >&2
    echo "decompressing ${krona_taxonomy_db_tgz}" >&2
    decompressor="pigz -dc"
    if [[ "${krona_taxonomy_db_tgz}" == *.lz4 ]]; then
      decompressor="lz4 -d"
    elif [[ "${krona_taxonomy_db_tgz}" == *.bz2 ]]; then
      decompressor="bzcat -d"
    fi
    cat ${krona_taxonomy_db_tgz} | $decompressor | tar -xv >&2

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
    memory: "3GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
  }
}


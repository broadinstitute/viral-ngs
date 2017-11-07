

task kraken_single {
  File reads_unmapped_bam
  File kraken_db_tar_lz4

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    date
    echo "decompressing $kraken_db_tar_lz4"
    decompressor="pigz -dc"
    if [[ "$kraken_db_tar_lz4" == *.lz4 ]]; then
      decompressor="lz4 -d"
    elif [[ "$kraken_db_tar_lz4" == *.bz2 ]]; then
      decompressor="bzcat -d"
    fi
    time cat $kraken_db_tar_lz4 | $decompressor | tar -C /mnt/db -xvf -

    time metagenomics.py kraken \
      ${reads_unmapped_bam} \
      /mnt/db \
      --outReads=/mnt/output/kraken-reads.txt.gz \
      --outReport=/mnt/output/kraken-report.txt
  }

  output {
    File kraken_classified_reads = /mnt/output/kraken-reads.txt.gz
    File kraken_summary_report = /mnt/output/kraken-report.txt
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "200GB"
    cpu: "32"
    disks: "local-disk 375 LOCAL, /mnt/output 375 LOCAL, /mnt/db 375 LOCAL"
  }
}

task krona {
  File classified_reads_txt_gz
  File krona_taxonomy_db_tgz

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    date
    echo "decompressing $krona_taxonomy_db_tgz"
    decompressor="pigz -dc"
    if [[ "$krona_taxonomy_db_tgz" == *.lz4 ]]; then
      decompressor="lz4 -d"
    elif [[ "$krona_taxonomy_db_tgz" == *.bz2 ]]; then
      decompressor="bzcat -d"
    fi
    mkdir -p krona_db krona_out
    cat $krona_taxonomy_db_tgz | $decompressor | tar -C krona_db -xvf -

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      krona_db \
      krona_out/krona-report.html \
      --noRank --noHits

    tar czf krona-report.tar.gz -C krona_out .
  }

  output {
    File krona_report_html = krona_out/krona-report.html
    File krona_report_tgz = krona-report.tar.gz
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "3GB"
    cpu: "2"
    disks: "local-disk 375 LOCAL"
  }
}


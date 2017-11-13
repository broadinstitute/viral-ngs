
# TO DO:
# kraken_multi (Array[File] bam inputs and reads/reports outputs)
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken {
  Array[File]+ reads_unmapped_bam
  File kraken_db_tar_lz4 # pipeable

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    cat ${kraken_db_tar_lz4} |
      read_utils.py extract_tarball \
        - /mnt/db \
        --pipe_hint=${kraken_db_tar_lz4} \
        --loglevel=DEBUG

    # prep input and output file names
    IN_BAMS=${write_lines(reads_unmapped_bam)}
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    for bam in `cat $IN_BAMS`; do
      echo "kraken-reads-$(basename $bam .bam).txt.gz" >> $OUT_READS
      echo "kraken-report-$(basename $bam .bam).txt" >> $OUT_REPORTS
    done

    time metagenomics.py kraken \
      /mnt/db \
      `cat $IN_BAMS` \
      --outReads=`cat $OUT_READS` \
      --outReport=`cat $OUT_REPORTS` \
      --loglevel=DEBUG
  }

  output {
    Array[File] kraken_classified_reads = read_lines("fnames_outreads.txt")
    Array[File] kraken_summary_report   = read_lines("fnames_outreports.txt")
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

    # prep output file names
    OUT_HTML="fname-out_html.txt"
    OUT_TGZ="fname-out_tgz.txt"
    echo "$(basename ${classified_reads_txt_gz} .txt.gz).html" > $OUT_HTML
    echo "$(basename ${classified_reads_txt_gz} .txt.gz).krona.tar.gz" > $OUT_TGZ

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      taxonomy \
      `cat $OUT_HTML` \
      --noRank --noHits \
      --loglevel=DEBUG

    tar czf `cat $OUT_TGZ` `cat $OUT_HTML`*
  }

  output {
    File krona_report_html = read_string("fname-out_html.txt")
    File krona_report_tgz  = read_string("fname-out_tgz.txt")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "2GB"
    cpu: 1
  }
}


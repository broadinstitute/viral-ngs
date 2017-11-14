
# TO DO:
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken {
  Array[File]+ reads_unmapped_bam
  File kraken_db_tar_lz4

  parameter_meta {
    kraken_db_tar_lz4:  "stream" # for DNAnexus, until WDL implements the File| type
    reads_unmapped_bam: "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    cat ${kraken_db_tar_lz4} |
      read_utils.py extract_tarball \
        - /mnt/db \
        --pipe_hint=${kraken_db_tar_lz4} \
        --loglevel=DEBUG

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "kraken-reads-$(basename $bam .bam).txt.gz" >> $OUT_READS
      echo "kraken-report-$(basename $bam .bam).txt" >> $OUT_REPORTS
    done

    metagenomics.py kraken \
      /mnt/db \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    ls -alF `cat $OUT_READS $OUT_REPORTS`
  }

  output {
    Array[File] kraken_classified_reads = glob("kraken-reads-*.txt.gz")
    Array[File] kraken_summary_report   = glob("kraken-report-*.txt")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "200 GB"
    cpu: 32
    disks: "local-disk 375 LOCAL, /mnt/db 375 LOCAL"
    preemptible: 0
  }
}

task krona {
  File classified_reads_txt_gz
  File? krona_taxonomy_db_tgz = "gs://sabeti-public-dbs/krona/krona_taxonomy_20160502.tar.lz4" # pipeable

  parameter_meta {
    krona_taxonomy_db_tgz : "stream" # for DNAnexus, until WDL implements the File| type
  }

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
#    File krona_report_html = "${base_outfname}.html"
#    File krona_report_tgz  = "${base_outfname}.krona.tar.gz"
    File krona_report_html = select_first(glob("*.html"))
    File krona_report_tgz  = select_first(glob("*.krona.tar.gz"))
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
  }
}


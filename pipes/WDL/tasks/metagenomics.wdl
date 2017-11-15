
# TO DO:
# kraken_build (input & output tarballs)
# diamond, bwa, etc

task kraken {
  Array[File]+ reads_unmapped_bam
  File kraken_db_tar_lz4

  parameter_meta {
    kraken_db_tar_lz4:  "stream" # for DNAnexus, until WDL implements the File| type
    #reads_unmapped_bam: "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    # decompress DB to /mnt/db
    cat ${kraken_db_tar_lz4} |
      read_utils.py extract_tarball \
        - /mnt/db \
        --pipe_hint=${kraken_db_tar_lz4} \
        --loglevel=DEBUG

    ## prep input and output file names
    #OUT_READS=fnames_outreads.txt
    #OUT_REPORTS=fnames_outreports.txt
    #for bam in ${sep=' ' reads_unmapped_bam}; do
    #  echo "kraken-reads-$(basename $bam .bam).txt.gz" >> $OUT_READS
    #  echo "kraken-report-$(basename $bam .bam).txt" >> $OUT_REPORTS
    #done
    #
    #metagenomics.py kraken \
    #  /mnt/db \
    #  ${sep=' ' reads_unmapped_bam} \
    #  --outReads `cat $OUT_READS` \
    #  --outReport `cat $OUT_REPORTS` \
    #  --loglevel=DEBUG

    # execute on each bam sequentially
    for bam in ${sep=' ' reads_unmapped_bam}; do
      metagenomics.py kraken \
        /mnt/db \
        $bam
        --outReads kraken-reads-$(basename $bam .bam).txt.gz \
        --outReport kraken-report-$(basename $bam .bam).txt \
        --loglevel=DEBUG
    done
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
  File   classified_reads_txt_gz
  File?  krona_taxonomy_db_tgz = "gs://sabeti-public-dbs/krona/krona_taxonomy_20160502.tar.lz4" # pipeable

  String input_basename = basename(classified_reads_txt_gz, ".txt.gz")

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

    metagenomics.py krona \
      ${classified_reads_txt_gz} \
      taxonomy \
      ${input_basename}.html
      --noRank --noHits \
      --loglevel=DEBUG

    tar czf ${input_basename}.krona.tar.gz ${input_basename}.html*
  }

  output {
    File krona_report_html = "${input_basename}.html"
    File krona_report_tgz  = "${input_basename}.krona.tar.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
  }
}



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
      --loglevel=DEBUG

    # prep input and output file names
    OUT_READS=fnames_outreads.txt
    OUT_REPORTS=fnames_outreports.txt
    for bam in ${sep=' ' reads_unmapped_bam}; do
      echo "kraken-reads-$(basename $bam .bam).txt.gz" >> $OUT_READS
      echo "kraken-report-$(basename $bam .bam).txt" >> $OUT_REPORTS
    done

    # execute on all inputs and outputs at once
    metagenomics.py kraken \
      $DB_DIR \
      ${sep=' ' reads_unmapped_bam} \
      --outReads `cat $OUT_READS` \
      --outReport `cat $OUT_REPORTS` \
      --loglevel=DEBUG

    for bam in ${sep=' ' reads_unmapped_bam}; do
      report_basename="kraken-reads-$(basename $bam .bam)"
      metagenomics.py krona \
        $report_basename.txt.gz \
        taxonomy \
        $report_basename.html \
        --noRank --noHits \
        --loglevel=DEBUG
      tar czf $report_basename.krona.tar.gz $report_basename.html*
    done
  }

  output {
    Array[File] kraken_classified_reads = glob("kraken-reads-*.txt.gz")
    Array[File] kraken_summary_report   = glob("kraken-report-*.txt")
    Array[File] krona_report_html       = glob("kraken-reads-*.html")
    Array[File] krona_report_tgz        = glob("kraken-reads-*.krona.tar.gz")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
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
    File krona_report_html = "${input_basename}.html"
    File krona_report_tgz  = "${input_basename}.krona.tar.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
  }
}


task diamond_contigs {
  File  contigs_fasta
  File  reads_unmapped_bam
  File  diamond_db_tar_lz4
  File  diamond_taxonomy_db_tar_lz4
  File  bwa_taxonomy_db_tar_lz4
  File  krona_taxonomy_db_tar_lz4

  String contigs_basename = basename(contigs_fasta, ".fasta")

  parameter_meta {
    diamond_db_tar_lz4          : "stream" # for DNAnexus, until WDL implements the File| type
    diamond_taxonomy_db_tar_lz4 : "stream" # for DNAnexus, until WDL implements the File| type
    bwa_taxonomy_db_tar_lz4     : "stream" # for DNAnexus, until WDL implements the File| type
    krona_taxonomy_db_tar_lz4   : "stream" # for DNAnexus, until WDL implements the File| type
  }

  command {
    set -ex -o pipefail

    echo "TO DO: this is not yet implemented"
    exit 1

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DIAMOND_DB_DIR=$(mktemp -d)
    DIAMOND_TAXDB_DIR=$(mktemp -d)
    BWA_TAXDB_DIR=$(mktemp -d)

    # decompress DBs to /mnt/db
    read_utils.py extract_tarball \
      ${diamond_db_tar_lz4} $DIAMOND_DB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${diamond_taxonomy_db_tar_lz4} $DIAMOND_TAXDB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${bwa_taxonomy_db_tar_lz4} $BWA_TAXDB_DIR \
      --loglevel=DEBUG
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tar_lz4} . \
      --loglevel=DEBUG

    # classify contigs
    metagenomics.py diamond_fasta \
      ${contigs_fasta} \
      $DIAMOND_DB_DIR \
      $DIAMOND_TAXDB_DIR \
      ${contigs_basename}.diamond.fasta \
      --loglevel=DEBUG

    # map reads to contigs & create kraken-like read report
    metagenomics.py align_rna \
      ${reads_unmapped_bam} \
      ${contigs_basename}.diamond.fasta \
      $BWA_TAXDB_DIR \
      ${contigs_basename}.diamond.summary_report.txt \
      --outReads ${contigs_basename}.diamond.reads.txt.gz \
      --dupeReads ${contigs_basename}.diamond.reads_w_dupes.txt.gz \
      --outBam ${contigs_basename}.diamond.bam \
      --loglevel=DEBUG

    # run krona
    metagenomics.py krona \
      ${contigs_basename}.diamond.reads.txt.gz \
      taxonomy \
      ${contigs_basename}.diamond.html \
      --noRank --noHits \
      --loglevel=DEBUG
    tar czf ${contigs_basename}.diamond.krona.tar.gz ${contigs_basename}.diamond.html*
  }

  output {
    File diamond_contigs = "${contigs_basename}.diamond.fasta"
    File reads_mapped_to_contigs = "${contigs_basename}.diamond.bam"
    File diamond_summary_report = "${contigs_basename}.diamond.summary_report.txt"
    File diamond_classified_reads = "${contigs_basename}.diamond.reads.txt.gz"
    File diamond_classified_reads_w_dupes = "${contigs_basename}.diamond.reads_w_dupes.txt.gz"
    File krona_report_html = "${contigs_basename}.diamond.html"
    File krona_report_tgz  = "${contigs_basename}.diamond.krona.tar.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "100 GB"
    cpu: 16
    dx_instance_type: "mem3_ssd1_x16"
  }
}

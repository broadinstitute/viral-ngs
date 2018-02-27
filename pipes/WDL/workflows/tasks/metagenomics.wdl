
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
    File krona_report_html = "${input_basename}.html"
    File krona_report_tgz  = "${input_basename}.krona.tar.gz"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "4 GB"
    cpu: 1
    dx_instance_type: "mem2_hdd2_x2"
  }
}


task download_blastn_nt_db {
  File? nt_tgz
  File? tax_db_tgz
  File? blast_tax_db_tgz

  command {
    set -ex -o pipefail
    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)
    TAX_DB_DIR=$(mktemp -d)

    DB=${nt_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      lftp -c "mirror --exclude-glob * --include-glob nt.*.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/ $DB_DIR"
      pigz -d $DB_DIR/*.gz
    else
      tar -zx -C "$DB_DIR" -f "$DB"
    fi

    DB=${blast_tax_db_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar -zx -C $DB_DIR
    else
      tar -zx -C "$DB_DIR" -f "$DB"
    fi

    DB=${tax_db_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -zx -C $TAX_DB_DIR
    else
      tar -zx -C "$TAX_DB_DIR" -f "$DB"
    fi

    find $DB_DIR -name "*.???" > blast_db_files.txt
    find $TAX_DB_DIR -type f > tax_db_files.txt
  }

  output {
    Array[File]+ db = read_lines("blast_db_files.txt")
    Array[File]+ tax_db = read_lines("tax_db_files.txt")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_hdd2_x1"
  }
}

task blastn_contigs {
  File contigs_fasta
  Array[File]+ blast_db_files
  Array[File]+ tax_db_files

  String contigs_basename = basename(contigs_fasta, ".fasta")
  String dirname = sub(blast_db_files[0], "/[^/]*$", "")
  String tax_dirname = sub(tax_db_files[0], "/[^/]*$", "")
  String db_basename = "nt"

  command {

    metagenomics.py blast_taxonomy ${contigs_fasta} \
      --taxDb ${tax_dirname} \
      --ntDb ${dirname + "/" + db_basename} \
      --outBlastn ${contigs_basename}.blastn.m8.gz \
      --outBlastnLca ${contigs_basename}.blastn.lca.tsv.gz \
      --outBlastnReport ${contigs_basename}.blastn.report
  }

  output {
    File blastn = "${contigs_basename}.blastn.m8.gz"
    File blastn_lca = "${contigs_basename}.blastn.lca.tsv.gz"
    File blastn_report = "${contigs_basename}.blastn.report"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "32 GB"
    cpu: 16
    dx_instance_type: "mem1_hdd2_x16"
  }
}

task download_blastx_nr_db {
  File? nr_tgz
  File? tax_db_tgz
  File? blast_tax_db_tgz

  command {
    set -ex -o pipefail
    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)
    TAX_DB_DIR=$(mktemp -d)

    DB=${nr_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      lftp -c "mirror --exclude-glob * --include-glob nr.*.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/db/ $DB_DIR"
      pigz -d $DB_DIR/*.gz
    else
      tar -zx -C "$DB_DIR" -f "$DB"
    fi

    DB=${blast_tax_db_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar -zx -C $DB_DIR
    else
      tar -zx -C "$DB_DIR" -f "$DB"
    fi

    DB=${tax_db_tgz}
    if [ -z "$DB" ]; then
      # Download latest databases from NCBI FTP
      curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -zx -C $TAX_DB_DIR
    else
      tar -zx -C "$TAX_DB_DIR" -f "$DB"
    fi

    find $DB_DIR -name "*.???" > blast_db_files.txt
    find $TAX_DB_DIR -type f > tax_db_files.txt
  }

  output {
    Array[File]+ db = read_lines("blast_db_files.txt")
    Array[File]+ tax_db = read_lines("tax_db_files.txt")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_hdd2_x1"
  }
}


task blastx_contigs {
  File contigs_fasta
  Array[File]+ blast_db_files
  Array[File]+ tax_db_files

  String contigs_basename = basename(contigs_fasta, ".fasta")
  String dirname = sub(blast_db_files[0], "/[^/]*$", "")
  String tax_dirname = sub(tax_db_files[0], "/[^/]*$", "")
  String db_basename = "nr"

  command {
    metagenomics.py blast_taxonomy ${contigs_fasta} \
      --taxDb ${tax_dirname} \
      --nrDb ${dirname + "/" + db_basename} \
      --outBlastx ${contigs_basename}.blastx.m8.gz \
      --outBlastxLca ${contigs_basename}.blastx.lca.tsv.gz \
      --outBlastxReport ${contigs_basename}.blastx.report
  }

  output {
    File blastx = "${contigs_basename}.blastx.m8.gz"
    File blastx_lca = "${contigs_basename}.blastx.lca.tsv.gz"
    File blastx_report = "${contigs_basename}.blastx.report"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "32 GB"
    cpu: 16
    dx_instance_type: "mem1_hdd2_x16"
  }
}

task download_cdd_db {
  File? cdd_db_tgz
  command {
    set -ex -o pipefail
    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    DB_DIR=$(mktemp -d)

    DB=${cdd_db_tgz}
    if [ -z "$DB" ]; then
      curl -s ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz | tar -zx -C "$DB_DIR"
    else
      tar -zx -C "$DB_DIR" -f "$DB"
    fi

    echo "$DB_DIR"/Cdd
  }

  output {
    String db = read_string(stdout())
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_hdd2_x1"
  }
}


task download_rfam_db {
  File? rfam_db_tgz
  command {
    DB=${rfam_db_tgz}
    if [ -z "$DB" ]; then
      curl -s ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz | unpigz > Rfam.cm
      cmpress Rfam.cm 1>&2
    else
      tar -zxf "$DB"
    fi
  }

  output {
    Array[File] db_files = glob("Rfam.cm*")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_hdd2_x1"
  }
}

task rpsblast_models {
  File contigs_fasta
  String cdd_db

  String contigs_basename = basename(contigs_fasta, ".fasta")
  command {
    metagenomics.py rpsblast_models \
      ${cdd_db} \
      ${contigs_fasta} \
      ${contigs_basename}.cdd.report
  }

  output {
    File cdd_report = "${contigs_basename}.cdd.report"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "32 GB"
    cpu: 16
    dx_instance_type: "mem1_hdd2_x16"
  }
}


task infernal_rfam {
  File contigs_fasta
  Array[File] rfam_db_files

  String rfam_db = rfam_db_files[0]

  String contigs_basename = basename(contigs_fasta, ".fasta")
  command {
    metagenomics.py infernal_contigs \
      ${rfam_db} \
      ${contigs_fasta} \
      ${contigs_basename}.infernal.tbl
  }

  output {
    File infernal_tbl = "${contigs_basename}.infernal.tbl"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "32 GB"
    cpu: 16
    dx_instance_type: "mem1_hdd2_x16"
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
    mem_in_gb=`/opt/viral-ngs/source/docker/mem_in_gb_90.sh`

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
    File diamond_contigs = "${contigs_basename}.diamond.fasta"
    File reads_mapped_to_contigs = "${contigs_basename}.diamond.bam"
    File diamond_summary_report = "${contigs_basename}.diamond.summary_report.txt"
    File diamond_classified_reads = "${contigs_basename}.diamond.reads.txt.gz"
    File diamond_classified_reads_w_dupes = "${contigs_basename}.diamond.reads_w_dupes.txt.gz"
    File krona_report_html = "${contigs_basename}.diamond.html"
    File krona_report_tgz  = "${contigs_basename}.diamond.krona.tar.gz"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "100 GB"
    cpu: 16
    dx_instance_type: "mem3_ssd1_x16"
  }
}

import "tasks/metagenomics.wdl" as metagenomics


workflow classify_kraken_single {
  File reads_unmapped_bam
  File kraken_db_tarball
  File krona_db_tarball

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = [reads_unmapped_bam],
      kraken_db_tar_lz4 = kraken_db_tarball
  }

  call metagenomics.krona as krona {
    input:
      classified_reads_txt_gz = select_first(kraken.kraken_classified_reads),
      krona_taxonomy_db_tgz = krona_db_tarball
  }
}

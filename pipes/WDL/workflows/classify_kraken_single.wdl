import "tasks/metagenomics.wdl" as metagenomics

# DX_SKIP_WORKFLOW
workflow classify_kraken_single {
  File reads_unmapped_bam

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = [reads_unmapped_bam]
  }

  call metagenomics.krona as krona {
    input:
      classified_reads_txt_gz = select_first(kraken.kraken_classified_reads)
  }
}

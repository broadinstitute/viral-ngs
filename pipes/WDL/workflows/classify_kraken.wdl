import "metagenomics.wdl" as metagenomics

workflow classify_kraken {
    Array[File] reads_unmapped_bam
    File        kraken_db_tar_lz4
    File        krona_taxonomy_db_tgz
    call metagenomics.kraken {
        input:
            reads_unmapped_bam = reads_unmapped_bam,
            kraken_db_tar_lz4 = kraken_db_tar_lz4,
            krona_taxonomy_db_tgz = krona_taxonomy_db_tgz
    }
}

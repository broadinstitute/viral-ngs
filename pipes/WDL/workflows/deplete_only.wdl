import "taxon_filter.wdl" as taxon_filter

workflow deplete_only {
    File raw_reads_unmapped_bam
    Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? bwaDbs
    call taxon_filter.deplete_taxa {
        input:
            raw_reads_unmapped_bam = raw_reads_unmapped_bam,
            bmtaggerDbs = bmtaggerDbs,
            blastDbs = blastDbs,
            bwaDbs = bwaDbs
    }
}

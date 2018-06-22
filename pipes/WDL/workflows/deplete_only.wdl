import "taxon_filter.wdl" as taxon_filter

workflow deplete_only {
    File raw_reads_unmapped_bam
    call taxon_filter.deplete_taxa {
        input:
            raw_reads_unmapped_bam = raw_reads_unmapped_bam
    }
}

import "tasks_metagenomics.wdl" as metagenomics

workflow filter_classified_bam_to_taxa {
    call metagenomics.filter_bam_to_taxa
}

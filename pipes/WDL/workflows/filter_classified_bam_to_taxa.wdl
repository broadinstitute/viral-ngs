import "tasks_metagenomics.wdl" as metagenomics
import "tasks_reports.wdl" as reports

workflow filter_classified_bam_to_taxa {
    call metagenomics.filter_bam_to_taxa
    call reports.software_version
}

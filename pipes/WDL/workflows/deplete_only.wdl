import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_reports.wdl" as reports

workflow deplete_only {
    call taxon_filter.deplete_taxa
    call reports.software_version
}

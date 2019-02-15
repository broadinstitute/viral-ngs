import "tasks_metagenomics.wdl" as metagenomics
import "tasks_reports.wdl" as reports

workflow classify_kraken {
    call metagenomics.kraken
    call reports.software_version
}

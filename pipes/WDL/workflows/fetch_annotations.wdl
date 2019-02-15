import "tasks_ncbi.wdl" as ncbi
import "tasks_reports.wdl" as reports

workflow fetch_annotations {
    call ncbi.download_annotations
    call reports.software_version
}

import "tasks_ncbi.wdl" as ncbi

workflow fetch_annotations {
    call ncbi.download_annotations
}

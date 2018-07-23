import "ncbi.wdl" as ncbi

workflow fetch_annotations {
    call ncbi.download_annotations
}

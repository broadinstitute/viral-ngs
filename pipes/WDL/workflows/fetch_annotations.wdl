import "ncbi.wdl" as ncbi

workflow fetch_annotations {
    Array[String]+ accessions
    String         emailAddress
    String         combined_out_prefix
    call ncbi.download_annotations as download {
        input:
            accessions = accessions,
            emailAddress = emailAddress,
            combined_out_prefix = combined_out_prefix
    }
}

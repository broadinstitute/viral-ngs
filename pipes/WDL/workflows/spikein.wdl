import "reports.wdl" as reports

workflow spikein {
    File reads_bam
    File  spikein_db
    call reports.spikein_report as spikein_report {
        input:
            reads_bam = reads_bam,
            spikein_db = spikein_db
    }
}

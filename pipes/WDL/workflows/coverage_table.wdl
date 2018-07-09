import "reports.wdl" as reports

workflow coverage_table {
    Array[File]+ mapped_bams
    Array[File] mapped_bam_idx # optional.. speeds it up if you provide it, otherwise we auto-index
    call reports.coverage_report {
        input:
            mapped_bams = mapped_bams,
            mapped_bam_idx = mapped_bam_idx
    }
}

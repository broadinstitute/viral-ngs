import "tasks_demux.wdl" as demux
import "tasks_reports.wdl" as reports

workflow merge_bams {
    call demux.merge_and_reheader_bams
    call reports.software_version
}

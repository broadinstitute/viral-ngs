import "tasks_read_utils.wdl" as reads
import "tasks_reports.wdl" as reports

workflow downsample {
    call reads.downsample_bams
    call reports.software_version
}

import "tasks_demux.wdl" as tasks_demux
import "tasks_reports.wdl" as reports

workflow demux_only {
    call tasks_demux.illumina_demux
    call reports.software_version
}

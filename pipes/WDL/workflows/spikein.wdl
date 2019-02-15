import "tasks_reports.wdl" as reports
import "tasks_reports.wdl" as reports

workflow spikein {
    call reports.spikein_report

    call reports.software_version
}

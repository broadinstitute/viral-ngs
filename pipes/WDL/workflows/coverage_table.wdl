import "tasks_reports.wdl" as reports

workflow coverage_table {
    call reports.coverage_report
    call reports.software_version
}

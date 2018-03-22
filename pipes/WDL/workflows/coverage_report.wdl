import "reports.wdl" as reports

workflow coverage_report {
  call reports.coverage_report as coverage
}

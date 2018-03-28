import "reports.wdl" as reports

workflow spikein {

  call reports.spikein_report as spikein_report

}

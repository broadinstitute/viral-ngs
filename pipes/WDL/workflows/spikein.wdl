import "tasks/reports.wdl" as reports

workflow spikein {

  call reports.spikein_report as spikein

}

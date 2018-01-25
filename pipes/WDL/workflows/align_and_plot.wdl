import "tasks/reports.wdl" as reports

workflow align_and_plot {
  call reports.plot_coverage
}

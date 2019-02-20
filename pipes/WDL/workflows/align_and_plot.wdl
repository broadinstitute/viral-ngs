import "tasks_reports.wdl" as reports

workflow align_and_plot {
    call reports.plot_coverage {
        input:
            aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502"
    }
}

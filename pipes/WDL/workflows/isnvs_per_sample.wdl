import "tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_per_sample {
    call tasks_intrahost.isnvs_per_sample
}

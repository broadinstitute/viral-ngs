import "tasks_intrahost.wdl" as intrahost

workflow isnvs_one_sample {
    call intrahost.isnvs_per_sample
}

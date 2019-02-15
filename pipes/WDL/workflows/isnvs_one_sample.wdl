import "tasks_intrahost.wdl" as intrahost
import "tasks_reports.wdl" as reports

workflow isnvs_one_sample {
    call intrahost.isnvs_per_sample
    call reports.software_version
}

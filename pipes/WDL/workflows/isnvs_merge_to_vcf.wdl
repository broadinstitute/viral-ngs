import "tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_merge_to_vcf {
    call tasks_intrahost.isnvs_vcf
}

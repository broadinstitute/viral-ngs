import "tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_merge_to_vcf_filtered {
    call tasks_intrahost.isnvs_vcf_filtered
}

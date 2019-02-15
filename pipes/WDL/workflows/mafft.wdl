import "tasks_interhost.wdl" as interhost
import "tasks_reports.wdl" as reports

workflow mafft {
    call interhost.multi_align_mafft
    call reports.software_version
}

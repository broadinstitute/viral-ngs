import "tasks_interhost.wdl" as interhost

workflow mafft {
    call interhost.multi_align_mafft
}

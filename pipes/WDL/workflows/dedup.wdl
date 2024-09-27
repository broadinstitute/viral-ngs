import "tasks_read_utils.wdl" as reads

workflow dedup {
    call reads.dedup_bam
}

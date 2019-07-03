import "tasks_read_utils.wdl" as reads

workflow dedup_bam {
    call reads.dedup_bam
}

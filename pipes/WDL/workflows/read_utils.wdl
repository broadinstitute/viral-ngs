import "tasks_read_utils.wdl" as reads

workflow downsample {
    call reads.downsample_bams
}

workflow dedup_bam {
    call reads.dedup_bam
}

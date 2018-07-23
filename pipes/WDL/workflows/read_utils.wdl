import "tasks_read_utils.wdl" as reads

workflow downsample {
    call reads.downsample_bams
}

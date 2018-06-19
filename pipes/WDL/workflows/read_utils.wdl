import "tasks_read_utils.wdl" as reads

workflow downsample {
    Array[File] reads_bam

    call reads.downsample_bams {
        input:
            reads_bam = reads_bam
    }
}

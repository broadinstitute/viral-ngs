import "read_utils.wdl" as read_utils

workflow downsample {
    Array[File] reads_bam

    call read_utils.downsample_bams {
        input:
            reads_bam = reads_bam
    }
}

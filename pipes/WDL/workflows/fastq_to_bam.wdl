import "tasks_read_utils.wdl" as read_utils

workflow paired_fastqs_to_bam {
    call read_utils.fastq_to_bam
}

import "tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_one_sample {
    File mapped_bam
    File assembly_fasta
    call tasks_intrahost.isnvs_per_sample {
        input:
            mapped_bam = mapped_bam,
            assembly_fasta = assembly_fasta
    }
}

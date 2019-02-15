import "tasks_assembly.wdl" as assembly
import "tasks_reports.wdl" as reports

workflow scaffold_and_refine {
    File reads_unmapped_bam

    call assembly.scaffold {
        input:
            reads_bam = reads_unmapped_bam
    }

    call assembly.refine_2x_and_plot {
        input:
            assembly_fasta = scaffold.scaffold_fasta,
            reads_unmapped_bam = reads_unmapped_bam
    }

    call reports.software_version
}

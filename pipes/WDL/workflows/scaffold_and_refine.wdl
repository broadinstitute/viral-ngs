import "assembly.wdl" as assembly

workflow scaffold_and_refine {
    File reads_unmapped_bam
    File contigs_fasta
    Array[File]+ reference_genome_fasta
    call assembly.scaffold {
        input:
            reads_bam = reads_unmapped_bam,
            contigs_fasta = contigs_fasta,
            reference_genome_fasta = reference_genome_fasta
    }

    File gatk_jar
    call assembly.refine_2x_and_plot {
        input:
            assembly_fasta = scaffold.scaffold_fasta,
            reads_unmapped_bam = reads_unmapped_bam,
            gatk_jar = gatk_jar
    }
}

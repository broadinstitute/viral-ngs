import "reports.wdl" as reports

workflow align_and_plot {
    String sample_name

    File assembly_fasta
    File reads_unmapped_bam

    File gatk_jar

    call reports.plot_coverage {
        input:
            aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502",
            sample_name = sample_name,
            assembly_fasta = assembly_fasta,
            reads_unmapped_bam = reads_unmapped_bam,
            gatk_jar = gatk_jar
    }
}

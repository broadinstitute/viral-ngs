import "assembly.wdl" as assembly

workflow assemble_refbased {

  File reads_unmapped_bam
  File ref_fasta
  File gatk_jar
  call assembly.refine_2x_and_plot{
    input:
        assembly_fasta = ref_fasta,
        reads_unmapped_bam = reads_unmapped_bam,
        gatk_jar = gatk_jar
  }
}

import "tasks_assembly.wdl" as assembly

workflow assemble_denovo_from_taxfilt {

  meta { 
     description: "Assemble denovo where cleaned reads have already been taxon-filtered"
  }
  
  File reads_unmapped_bam
  File taxfilt_bam

  call assembly.assemble {
    input:
      reads_unmapped_bam = taxfilt_bam
  }

  call assembly.scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = taxfilt_bam
  }

  call assembly.refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = reads_unmapped_bam
  }
}
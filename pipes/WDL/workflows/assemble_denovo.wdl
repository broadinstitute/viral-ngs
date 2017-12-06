import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/assembly.wdl" as assembly

workflow assemble_denovo {
  File reads_unmapped_bam

  call taxon_filter.filter_to_taxon {
    input:
      reads_unmapped_bam = reads_unmapped_bam
  }

  call assembly.assemble {
    input:
      reads_unmapped_bam = filter_to_taxon.taxfilt_bam
  }

  call assembly.scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = filter_to_taxon.taxfilt_bam
  }

  call assembly.refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = reads_unmapped_bam
  }
}
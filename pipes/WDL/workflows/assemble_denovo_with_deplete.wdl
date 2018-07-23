import "tasks_taxon_filter.wdl" as filter
import "tasks_assembly.wdl" as assembly

workflow assemble_denovo_with_deplete {
  
  call filter.deplete_taxa

  call filter.filter_to_taxon {
    input:
      reads_unmapped_bam = deplete_taxa.cleaned_bam
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
      reads_unmapped_bam = deplete_taxa.cleaned_bam
  }
}

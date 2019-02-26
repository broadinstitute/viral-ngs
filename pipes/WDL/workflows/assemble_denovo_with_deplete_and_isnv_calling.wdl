import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_assembly.wdl" as assembly
import "tasks_intrahost.wdl" as intrahost

workflow assemble_denovo_with_deplete_and_isnv_calling {
  
  call taxon_filter.deplete_taxa

  call taxon_filter.filter_to_taxon {
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

  call intrahost.isnvs_per_sample {
    input:
      assembly_fasta = refine_2x_and_plot.final_assembly_fasta,
      mapped_bam = refine_2x_and_plot.aligned_bam
  }

}

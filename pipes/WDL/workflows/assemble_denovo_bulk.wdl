import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_assembly.wdl" as assembly

workflow assemble_denovo_bulk {
  
  File reads_unmapped_bam
  Array[File]+ testing_scatter_files
  
  scatter(testing_scatter in testing_scatter_files)
  {
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
}
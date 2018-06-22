import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly

workflow assemble_denovo {
  File reads_unmapped_bam
  File lastal_db_fasta

  call taxon_filter.filter_to_taxon {
    input:
      reads_unmapped_bam = reads_unmapped_bam,
      lastal_db_fasta = lastal_db_fasta
  }

  File trim_clip_db
  call assembly.assemble {
    input:
      reads_unmapped_bam = filter_to_taxon.taxfilt_bam,
      trim_clip_db = trim_clip_db
  }

  call assembly.scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = filter_to_taxon.taxfilt_bam
  }

  File gatk_jar
  call assembly.refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = reads_unmapped_bam,
      gatk_jar = gatk_jar
  }
}
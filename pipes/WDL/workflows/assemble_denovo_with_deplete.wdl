import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly

workflow assemble_denovo_with_deplete {
  
  File raw_reads_unmapped_bam
  Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
  Array[File]? bwaDbs
  call taxon_filter.deplete_taxa {
    input:
      raw_reads_unmapped_bam = raw_reads_unmapped_bam,
      bmtaggerDbs = bmtaggerDbs,
      blastDbs = blastDbs,
      bwaDbs = bwaDbs
  }

  File lastal_db_fasta
  call taxon_filter.filter_to_taxon {
    input:
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
      lastal_db_fasta = lastal_db_fasta
  }

  File trim_clip_db
  call assembly.assemble {
    input:
      reads_unmapped_bam = filter_to_taxon.taxfilt_bam,
      trim_clip_db = trim_clip_db
  }

  Array[File]+ reference_genome_fasta
  call assembly.scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = filter_to_taxon.taxfilt_bam,
      reference_genome_fasta = reference_genome_fasta
  }

  File gatk_jar
  call assembly.refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
      gatk_jar = gatk_jar
  }
}
import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/assembly.wdl" as assembly
import "tasks/reports.wdl" as reports

workflow assemble_denovo_with_deplete {
  File gatk_tar_bz2
  File? novocraft_license

  call taxon_filter.deplete_taxa

  call taxon_filter.filter_to_taxon {
    input:
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
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

  call assembly.refine as refine1 {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 30 -g 40 -x 20 -t 502",
      novocraft_license = novocraft_license,
      min_coverage = 2
  }

  call assembly.refine as refine2 {
    input:
      assembly_fasta = refine1.refined_assembly_fasta,
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100",
      novocraft_license = novocraft_license,
      min_coverage = 3
  }

  call reports.plot_coverage {
    input:
      assembly_fasta = refine2.refined_assembly_fasta,
      reads_unmapped_bam = deplete_taxa.cleaned_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novocraft_license = novocraft_license
  }
}
import "../../tasks/depletion.wdl" as depletion
import "../../tasks/assembly.wdl" as assembly


workflow assemble_denovo {
  String sample_name
  File reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  File trim_clip_db # fasta
  File lastal_db  # fasta

  call depletion.filter_to_taxon {
    input:
      sample_name = sample_name,
      reads = reads_unmapped_bam,
      lastal_db = lastal_db
  }

  call assembly.assemble_denovo {
    input:
      sample_name = sample_name,
      reads_unmapped_bam = filter_to_taxon.taxfilt_bam,
      trim_clip_db = trim_clip_db
  }

  call assembly.scaffold {
    input:
      sample_name = sample_name,
      trinity_contigs_fasta = assemble_denovo.contigs_fasta,
      trinity_reads_unmapped_bam = filter_to_taxon.taxfilt_bam
  }

  call assembly.refine as refine1 {
    input:
      sample_name = sample_name,
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 30 -g 40 -x 20 -t 502",
      novocraft_license = novocrat_license,
      min_coverage = 2
  }

  call assembly.refine as refine2 {
    input:
      sample_name = sample_name,
      assembly_fasta = refine1.refined_assembly_fasta,
      reads_unmapped_bam = reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100",
      novocraft_license = novocrat_license,
      min_coverage = 3
  }

  call assembly.analysis {
    input:
      sample_name = sample_name,
      assembly_fasta = refine2.refined_assembly_fasta,
      reads_unmapped_bam = reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100 -k",
      novocraft_license = novocrat_license,
  }
}
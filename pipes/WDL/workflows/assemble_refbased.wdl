import "tasks/assembly.wdl" as assembly
import "tasks/reports.wdl" as reports

workflow assemble_refbased {
  String sample_name
  File reads_unmapped_bam
  File scaffold_fasta

  File gatk_tar_bz2
  File? novocraft_license

  call assembly.refine as refine1 {
    input:
      sample_name = sample_name,
      assembly_fasta = scaffold_fasta,
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

  call reports.plot_coverage {
    input:
      sample_name = sample_name,
      assembly_fasta = refine2.refined_assembly_fasta,
      reads_unmapped_bam = reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novocraft_license = novocrat_license
  }
}
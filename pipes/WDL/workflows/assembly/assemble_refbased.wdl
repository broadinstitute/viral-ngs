import "../../tasks/depletion.wdl" as depletion
import "../../tasks/assembly.wdl" as assembly


workflow assemble_denovo_with_deplete {
  String sample_name
  File raw_reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  call assembly.refine as refine1 {
    input:
      sample_name = sample_name,
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = raw_reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 30 -g 40 -x 20 -t 502",
      novocraft_license = novocrat_license,
      min_coverage = 2
  }

  call assembly.refine as refine2 {
    input:
      sample_name = sample_name,
      assembly_fasta = refine1.refined_assembly_fasta,
      reads_unmapped_bam = raw_reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100",
      novocraft_license = novocrat_license,
      min_coverage = 3
  }

  call assembly.analysis {
    input:
      sample_name = sample_name,
      assembly_fasta = refine2.refined_assembly_fasta,
      reads_unmapped_bam = raw_reads_unmapped_bam,
      gatk_tar_bz2 = gatk_tar_bz2,
      novoalign_options = "-r Random -l 40 -g 40 -x 20 -t 100 -k",
      novocraft_license = novocrat_license,
  }
}
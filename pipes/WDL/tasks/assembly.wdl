# ======================================================================
# assemble_trinity: 
#   Run the Trinity assembler.
#   First trim reads with trimmomatic, rmdup with prinseq, and
#   random subsample to no more than 100k reads
# ======================================================================

task assemble_denovo {
  # TO DO: update this to be able to perform trinity, spades or trinity+spades

  String sample_name

  File reads_unmapped_bam
  File trim_clip_db #fasta

  Int? trinity_n_reads=200000

  command {
    set -ex -o pipefail

    assembly.py assemble_trinity \
      ${reads_unmapped_bam} \
      ${trim_clip_db} \
      ${sample_name}.assembly1-trinity.fasta \
      ${'--n_reads=' + trinity_n_reads} \
      --JVMmemory 7g \
      --outReads=${sample_name}.subsamp.bam \
      --loglevel=DEBUG

    samtools view -c ${sample_name}.subsamp.bam | tee subsample_read_count
  }

  output {
    File contigs_fasta = "${sample_name}.assembly1-trinity.fasta"
    File subsampBam = "${sample_name}.subsamp.bam"
    Int subsample_read_count = read_int("subsample_read_count")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: 4
    disks: "local-disk 375 LOCAL"
  }

}

task scaffold {
  String sample_name
  File contigs_fasta
  File reads_bam
  File reference_genome_fasta

  String? aligner="muscle"
  Float? min_length_fraction=0.5
  Float? min_unambig=0.5
  Int? replace_length=55

  Int? assembly_nucmer_max_gap
  Int? assembly_nucmer_min_match
  Int? assembly_nucmer_min_cluster
  Int? scaffold_min_pct_contig_aligned

  command {
    set -ex -o pipefail

    assembly.py order_and_orient \
      ${contigs_fasta} \
      ${reference_genome_fasta} \
      ${sample_name}.intermediate_scaffold.fasta \
      ${'--maxgap=' + assembly_nucmer_max_gap} \
      ${'--minmatch=' + assembly_nucmer_min_match} \
      ${'--mincluster=' + assembly_nucmer_min_cluster} \
      ${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
      --loglevel=DEBUG

    assembly.py gapfill_gap2seq \
      ${sample_name}.intermediate_scaffold.fasta \
      ${reads_bam} \
      ${sample_name}.intermediate_gapfill.fasta \
      --memLimitGb 12 \
      --maskErrors \
      --loglevel=DEBUG

    assembly.py impute_from_reference \
      ${sample_name}.intermediate_gapfill.fasta \
      ${reference_genome_fasta} \
      ${sample_name}.scaffold.fasta \
      --newName ${sample_name} \
      --replaceLength ${replace_length} \
      --minLengthFraction ${min_length_fraction} \
      --minUnambig ${min_unambig} \
      --aligner ${aligner} \
      --loglevel=DEBUG
  }

  output {
    File scaffold_fasta = "${sample_name}.scaffold.fasta"
    File intermediate_scaffold_fasta = "${sample_name}.intermediate_scaffold.fasta"
    File intermediate_gapfill_fasta = "${sample_name}.intermediate_gapfill.fasta"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "12GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
  }
}

task refine {
  String sample_name
  File assembly_fasta
  File reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  String? novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
  Float? major_cutoff=0.5
  Int? min_coverage=1

  command {
    set -ex -o pipefail

    mkdir gatk/
    tar jxf ${gatk_tar_bz2} -C gatk/
    if [ -n "${novocraft_license}" ]; then
      cp ${novocraft_license} /tmp/novocraft.lic
      export NOVOALIGN_LICENSE_PATH=/tmp/novocraft.lic
    fi
    mv ${assembly_fasta} assembly.fasta
    read_utils.py novoindex assembly.fasta

    assembly.py refine_assembly \
      assembly.fasta \
      ${reads_unmapped_bam} \
      ${sample_name}.refined_assembly.fasta \
      --outVcf ${sample_name}.sites.vcf.gz \
      --min_coverage ${min_coverage} \
      --major_cutoff ${major_cutoff} \
      --GATK_PATH gatk/ \
      ${'--novo_params=' + novoalign_options} \
      --JVMmemory 7g \
      --loglevel=DEBUG
  }

  output {
    File refined_assembly_fasta = "${sample_name}.refined_assembly.fasta"
    File sites_vcf_gz = "${sample_name}.sites.vcf.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: 8
    disks: "local-disk 375 LOCAL"
  }
}


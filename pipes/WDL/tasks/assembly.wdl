# ======================================================================
# assemble_trinity: 
#   Run the Trinity assembler.
#   First trim reads with trimmomatic, rmdup with prinseq, and
#   random subsample to no more than 100k reads
# ======================================================================

task assemble_denovo {
  String sample

  File inputBam
  File trim_clip_db #fasta

  Int? trinity_n_reads=200000
  Int? threads

  command {
    assembly.py assemble_trinity \
      "${inputBam}" \
      "${trim_clip_db}" \
      "${sample}.assembly1-trinity.fasta" \
      "${'--n_reads' + trinity_n_reads}" \
      --threads $(nproc) \
      --JVMmemory 7g \
      --outReads="${sample}.subsamp.bam"
  }

  output {
    File trinityAssembly = "${sample}.assembly1-trinity.fasta"
    File subsampBam = "${sample}.subsamp.bam"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "4"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }

}

task scaffold {
  String sample_name
  File contigs_fasta
  File reads_bam
  File reference_genome_fasta

  String? aligner=muscle
  Float? min_length_fraction=0.5
  Float? min_unambig=0.5
  Int? replace_length=55

  Int? assembly_nucmer_max_gap
  Int? assembly_nucmer_min_match
  Int? assembly_nucmer_min_cluster
  Int? scaffold_min_pct_contig_aligned

  command <<<
    set -ex -o pipefail

    assembly.py order_and_orient \
      "${contigs_fasta}" \
      "${reference_genome_fasta}" \ "${sample_name}.intermediate_scaffold.fasta" \
      "${'--maxgap=' + assembly_nucmer_max_gap}"  \
      "${'--minmatch=' + assembly_nucmer_min_match}" \
      "${'--mincluster=' + assembly_nucmer_min_cluster}" \
      "${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned}"

    assembly.py gapfill_gap2seq \
      "${sample_name}.intermediate_scaffold.fasta" \
      "${reads_bam}" \
      "${sample_name}.intermediate_gapfill.fasta" \
      --memLimitGb 12 \
      --threads $(nproc) \
      --maskErrors

    assembly.py impute_from_reference \
      "${sample_name}.intermediate_gapfill.fasta" \
      "${reference_genome_fasta}" \
      "${sample_name}.scaffold.fasta" \
      --newName "${sample_name}" \
      --replaceLength ${replace_length} \
      --minLengthFraction ${min_length_fraction} \
      --minUnambig ${min_unambig} \
      --aligner ${aligner}
  >>>

  output {
    File scaffold_fasta = "${sample_name}.scaffold.fasta"
    File intermediate_scaffold_fasta = "${sample_name}.intermediate_scaffold.fasta"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "12GB"
    cpu: "4"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}

task refine {
  String sample_name
  File assembly_fasta
  File reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  String? novoalign_options
  Float? major_cutoff=0.5
  Int? min_coverage=1

  command <<<
    set -ex -o pipefail
    if [ -n "${novocraft_license}" ]; then
     cp ${novocraft_license} /tmp/novocraft.lic
     NOVO_LICENSE=" --NOVOALIGN_LICENSE_PATH /tmp/novocraft.lic"
    else
     NOVO_LICENSE=" "
    fi
    mkdir gatk/
    tar jxf "${gatk_tar_bz2}" -C gatk/
    mv "${assembly_fasta}" assembly.fasta
    novoindex assembly.nix assembly.fasta
    assembly.py refine_assembly assembly.fasta "${reads_unmapped_bam}" "${sample_name}.refined_assembly.fasta" \
      --outVcf "${sample_name}.sites.vcf.gz" --min_coverage ${min_coverage} --major_cutoff ${major_cutoff} \
      --threads $(nproc) --GATK_PATH gatk/ \
      --novo_params '${default="-r Random -l 40 -g 40 -x 20 -t 100" novoalign_options}' $NOVO_LICENSE
  >>>

  output {
    File refined_assembly_fasta = "${sample_name}.refined_assembly.fasta"
    File sites_vcf_gz = "${sample_name}.sites.vcf.gz"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "15GB"
    cpu: "8"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}

task map_reads_to_self {
  String sample

  File inputAssembly
  File cleanedBam

  Int? threads
  String? aligner
  String? alignerOptions

  command {
    read_utils.py align_and_fix \
      "${inputAssembly}" \
      --outBamAll "${sample}.bam" \
      --outBamFiltered "${sample}.mapped.bam" \
      "${'--aligner' + aligner || 'novoalign'}" \
      # set aligner options. The default, if novoalign, should be '-r Random -l 40 -g 40 -x 20 -t 100 -k'
      "${'--aligner_options' + alignerOptions}" \
      "${'--threads' + threads}" \
  }

  output {
    File outBam = "${sample}.bam"
    File outBamFiltered = "${sample}.mapped.bam"
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "8"
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}


task analysis {
  String sample_name

  File assembly_fasta
  File reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  String? novoalign_options

  command <<<
    set -ex -o pipefail
    if [ -n "${novocraft_license}" ]; then
     cp ${novocraft_license} /tmp/novocraft.lic
     NOVO_LICENSE=" --NOVOALIGN_LICENSE_PATH /tmp/novocraft.lic"
    else
     NOVO_LICENSE=" "
    fi
    mkdir gatk/
    tar jxf "${gatk_tar_bz2}" -C gatk/
    mv "${assembly_fasta}" assembly.fasta

    novoindex assembly.nix assembly.fasta
    read_utils.py index_fasta_picard assembly.fasta
    read_utils.py index_fasta_samtools assembly.fasta

    read_utils.py align_and_fix "${reads_unmapped_bam}" assembly.fasta \
      --outBamAll "${sample_name}.bam" --outBamFiltered "${sample_name}.mapped.bam" \
      --GATK_PATH gatk/ \
      --aligner_options "$novoalign_options" $NOVO_LICENSE
    samtools index "${sample_name}.mapped.bam"

    reports.py plot_coverage "${sample_name}.mapped.bam" "${sample_name}.coverage_plot.pdf" --plotFormat pdf --plotWidth 1100 --plotHeight 850 --plotDPI 100

    # collect figures of merit
    tail -n +1 assembly.fasta | tr -d '\n' | wc -c | tee assembly_length
    samtools view -c "${sample_name}.mapped.bam" | tee reads_aligned
    samtools flagstat "${sample_name}.bam" | tee "${sample_name}.bam.flagstat.txt"
    grep properly "${sample_name}.bam.flagstat.txt" | awk '{print $1}' | tee read_pairs_aligned
    samtools view "${sample_name}.mapped.bam" | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    expr $(cat bases_aligned) / $(cat assembly_length) | tee mean_coverage
  >>>

  output {
    File reads_bam = "${sample_name}.bam"
    File reads_bam_flagstat = "${sample_name}.bam.flagstat.txt"
    File coverage_plot = "${sample_name}.coverage_plot.pdf"
    Int assembly_length = read_int("assembly_length")
    Int reads_aligned = read_int("reads_aligned")
    Int read_pairs_aligned = read_int("read_pairs_aligned")
    Int bases_aligned = read_int("bases_aligned")
    Int mean_coverage = read_int("mean_coverage")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "4"
    preemptible: 1
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}



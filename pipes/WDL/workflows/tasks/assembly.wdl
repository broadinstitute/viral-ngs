
task assemble {

  File    reads_unmapped_bam
  File    trim_clip_db

  Int?    trinity_n_reads=250000
  Int?    spades_n_reads=10000000

  String? assembler="trinity"  # trinity, spades, or trinity-spades
  String  cleaned_assembler = select_first([assembler, ""]) # workaround for https://gatkforums.broadinstitute.org/wdl/discussion/10462/string-type-in-output-section

  String  sample_name = basename(reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_90.sh`
    mem_in_gb=`/opt/viral-ngs/source/docker/mem_in_gb_90.sh`

    if [[ "${assembler}" == "trinity" ]]; then
      assembly.py assemble_trinity \
        ${reads_unmapped_bam} \
        ${trim_clip_db} \
        ${sample_name}.assembly1-trinity.fasta \
        ${'--n_reads=' + trinity_n_reads} \
        --JVMmemory "$mem_in_mb"m \
        --outReads=${sample_name}.subsamp.bam \
        --loglevel=DEBUG

    elif [[ "${assembler}" == "spades" ]]; then
      assembly.py assemble_spades \
        ${reads_unmapped_bam} \
        ${trim_clip_db} \
        ${sample_name}.assembly1-spades.fasta \
        ${'--nReads=' + spades_n_reads} \
        --memLimitGb $mem_in_gb \
        --outReads=${sample_name}.subsamp.bam \
        --loglevel=DEBUG

    elif [[ "${assembler}" == "trinity-spades" ]]; then
      assembly.py assemble_trinity \
        ${reads_unmapped_bam} \
        ${trim_clip_db} \
        ${sample_name}.assembly1-trinity.fasta \
        ${'--n_reads=' + trinity_n_reads} \
        --JVMmemory "$mem_in_mb"m \
        --outReads=${sample_name}.subsamp.bam \
        --loglevel=DEBUG
      assembly.py assemble_spades \
        ${reads_unmapped_bam} \
        ${trim_clip_db} \
        ${sample_name}.assembly1-spades.fasta \
        --contigsUntrusted=${sample_name}.assembly1-trinity.fasta \
        ${'--nReads=' + spades_n_reads} \
        --memLimitGb $mem_in_gb \
        --loglevel=DEBUG

    else
      echo "unrecognized assembler ${assembler}" >&2
      exit 1
    fi

    samtools view -c ${sample_name}.subsamp.bam | tee subsample_read_count >&2
  }

  output {
    File contigs_fasta        = "${sample_name}.assembly1-${cleaned_assembler}.fasta"
    File subsampBam           = "${sample_name}.subsamp.bam"
    Int  subsample_read_count = read_int("subsample_read_count")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "15 GB"
    cpu: 4
    dx_instance_type: "mem1_ssd1_x8"
  }

}

task scaffold {
  File         contigs_fasta
  File         reads_bam
  Array[File]+ reference_genome_fasta

  String? aligner
  Float?  min_length_fraction
  Float?  min_unambig
  Int?    replace_length=55

  Int?    nucmer_max_gap
  Int?    nucmer_min_match
  Int?    nucmer_min_cluster
  Int?    scaffold_min_pct_contig_aligned

  String  sample_name = basename(contigs_fasta, ".fasta")

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_gb=`/opt/viral-ngs/source/docker/mem_in_gb_90.sh`

    assembly.py order_and_orient \
      ${contigs_fasta} \
      ${sep=' ' reference_genome_fasta} \
      ${sample_name}.intermediate_scaffold.fasta \
      ${'--maxgap=' + nucmer_max_gap} \
      ${'--minmatch=' + nucmer_min_match} \
      ${'--mincluster=' + nucmer_min_cluster} \
      ${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
      --outReference ${sample_name}.scaffolding_chosen_ref.fasta \
      --outStats ${sample_name}.scaffolding_stats.txt \
      --outAlternateContigs ${sample_name}.scaffolding_alt_contigs.fasta \
      --loglevel=DEBUG

    assembly.py gapfill_gap2seq \
      ${sample_name}.intermediate_scaffold.fasta \
      ${reads_bam} \
      ${sample_name}.intermediate_gapfill.fasta \
      --memLimitGb $mem_in_gb \
      --maskErrors \
      --loglevel=DEBUG

    grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\n' | wc -c | tee assembly_preimpute_length
    grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\nNn' | wc -c | tee assembly_preimpute_length_unambiguous

    assembly.py impute_from_reference \
      ${sample_name}.intermediate_gapfill.fasta \
      ${sample_name}.scaffolding_chosen_ref.fasta \
      ${sample_name}.scaffold.fasta \
      --newName ${sample_name} \
      ${'--replaceLength=' + replace_length} \
      ${'--minLengthFraction=' + min_length_fraction} \
      ${'--minUnambig=' + min_unambig} \
      ${'--aligner=' + aligner} \
      --loglevel=DEBUG
  }

  output {
    File scaffold_fasta              = "${sample_name}.scaffold.fasta"
    File intermediate_scaffold_fasta = "${sample_name}.intermediate_scaffold.fasta"
    File intermediate_gapfill_fasta  = "${sample_name}.intermediate_gapfill.fasta"
    Int  assembly_preimpute_length             = read_int("assembly_preimpute_length")
    Int  assembly_preimpute_length_unambiguous = read_int("assembly_preimpute_length_unambiguous")
    File scaffolding_chosen_ref      = "${sample_name}.scaffolding_chosen_ref.fasta"
    File scaffolding_stats           = "${sample_name}.scaffolding_stats.txt"
    File scaffolding_alt_contigs     = "${sample_name}.scaffolding_alt_contigs.fasta"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "15 GB"
    cpu: 4
    dx_instance_type: "mem1_ssd1_x8"
  }
}

task refine {
  File    assembly_fasta
  File    reads_unmapped_bam

  File    gatk_jar
  File?   novocraft_license

  String? novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
  Float?  major_cutoff=0.5
  Int?    min_coverage=1

  String  assembly_basename=basename(assembly_fasta, ".fasta")

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_90.sh`

    # prep GATK
    mkdir gatk
    if [[ ${gatk_jar} == *.tar.bz2 ]]; then
      tar -xjvf ${gatk_jar} -C gatk
    else
      ln -s ${gatk_jar} gatk/GenomeAnalysisTK.jar
    fi

    ln -s ${assembly_fasta} assembly.fasta
    read_utils.py novoindex assembly.fasta --loglevel=DEBUG

    assembly.py refine_assembly \
      assembly.fasta \
      ${reads_unmapped_bam} \
      ${assembly_basename}.refined.fasta \
      --outVcf ${assembly_basename}.sites.vcf.gz \
      --min_coverage ${min_coverage} \
      --major_cutoff ${major_cutoff} \
      --GATK_PATH gatk/ \
      --novo_params="${novoalign_options}" \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG
  }

  output {
    File refined_assembly_fasta = "${assembly_basename}.refined.fasta"
    File sites_vcf_gz           = "${assembly_basename}.sites.vcf.gz"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "7 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x8"
  }
}


task refine_2x_and_plot {
  # This combined task exists just to streamline the two calls to
  # assembly.refine and one call to reports.plot_coverage that almost
  # every assembly workflow uses. It saves on instance spin up and
  # docker pull times, file staging time, and all steps contained
  # here have similar hardware requirements. It is also extremely
  # rare for analyses to branch off of intermediate products between
  # these three steps.
  # The more atomic WDL tasks are still available for custom workflows.

  File    assembly_fasta
  File    reads_unmapped_bam

  File    gatk_jar  # can alternatively be the .tar.bz2
  File?   novocraft_license

  String? refine1_novoalign_options="-r Random -l 30 -g 40 -x 20 -t 502"
  Float?  refine1_major_cutoff=0.5
  Int?    refine1_min_coverage=2

  String? refine2_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
  Float?  refine2_major_cutoff=0.5
  Int?    refine2_min_coverage=3

  String? plot_coverage_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100 -k"

  String  assembly_basename = basename(assembly_fasta, ".fasta")
  String  sample_name       = basename(reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/mem_in_mb_90.sh`

    # prep GATK
    mkdir gatk
    if [[ ${gatk_jar} == *.tar.bz2 ]]; then
      tar -xjvf ${gatk_jar} -C gatk
    else
      ln -s ${gatk_jar} gatk/GenomeAnalysisTK.jar
    fi

    ln -s ${assembly_fasta} assembly.fasta
    read_utils.py novoindex assembly.fasta --loglevel=DEBUG

    # refine 1
    assembly.py refine_assembly \
      assembly.fasta \
      ${reads_unmapped_bam} \
      ${assembly_basename}.refine1.fasta \
      --outVcf ${assembly_basename}.refine1.pre_fasta.vcf.gz \
      --min_coverage ${refine1_min_coverage} \
      --major_cutoff ${refine1_major_cutoff} \
      --GATK_PATH gatk/ \
      --novo_params="${refine1_novoalign_options}" \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    # refine 2
    assembly.py refine_assembly \
      ${assembly_basename}.refine1.fasta \
      ${reads_unmapped_bam} \
      ${assembly_basename}.refine2.fasta \
      --outVcf ${assembly_basename}.refine2.pre_fasta.vcf.gz \
      --min_coverage ${refine2_min_coverage} \
      --major_cutoff ${refine2_major_cutoff} \
      --GATK_PATH gatk/ \
      --novo_params="${refine2_novoalign_options}" \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    # final alignment
    read_utils.py align_and_fix \
      ${reads_unmapped_bam} \
      ${assembly_basename}.refine2.fasta \
      --outBamAll ${sample_name}.bam \
      --outBamFiltered ${sample_name}.mapped.bam \
      --GATK_PATH gatk/ \
      --aligner_options "${plot_coverage_novoalign_options}" \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    # plot_coverage
    reports.py plot_coverage \
      ${sample_name}.mapped.bam \
      ${sample_name}.coverage_plot.pdf \
      --plotFormat pdf \
      --plotWidth 1100 \
      --plotHeight 850 \
      --plotDPI 100 \
      --loglevel=DEBUG

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
    samtools flagstat ${sample_name}.bam | tee ${sample_name}.bam.flagstat.txt
    grep properly ${sample_name}.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    expr $(cat bases_aligned) / $(cat assembly_length) | tee mean_coverage

    # fastqc mapped bam
    reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html
  }

  output {
    File refine1_sites_vcf_gz        = "${assembly_basename}.refine1.pre_fasta.vcf.gz"
    File refine1_assembly_fasta      = "${assembly_basename}.refine1.fasta"
    File refine2_sites_vcf_gz        = "${assembly_basename}.refine2.pre_fasta.vcf.gz"
    File refine2_assembly_fasta      = "${assembly_basename}.refine2.fasta"
    File aligned_bam                 = "${sample_name}.bam"
    File aligned_bam_flagstat        = "${sample_name}.bam.flagstat.txt"
    File aligned_only_reads_bam      = "${sample_name}.mapped.bam"
    File aligned_only_reads_fastqc   = "${sample_name}.mapped_fastqc.html"
    File coverage_plot               = "${sample_name}.coverage_plot.pdf"
    Int  assembly_length             = read_int("assembly_length")
    Int  assembly_length_unambiguous = read_int("assembly_length_unambiguous")
    Int  reads_aligned               = read_int("reads_aligned")
    Int  read_pairs_aligned          = read_int("read_pairs_aligned")
    Int  bases_aligned               = read_int("bases_aligned")
    Int  mean_coverage               = read_int("mean_coverage")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "7 GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_x8"
  }
}

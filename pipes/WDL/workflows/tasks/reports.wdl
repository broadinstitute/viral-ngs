


task plot_coverage {
  # TO DO: add a BWA option
  # TO DO: make GATK indel-realigner optional
  String sample_name

  File assembly_fasta
  File reads_unmapped_bam

  File gatk_tar_bz2
  File? novocraft_license

  String? aligner="novoalign" # novoalign or bwa
  String? aligner_options="-r Random -l 40 -g 40 -x 20 -t 100 -k"

  command {
    set -ex -o pipefail

    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v read_utils.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi

    read_utils.py extract_tarball ${gatk_tar_bz2} gatk
    cp ${assembly_fasta} assembly.fasta

    read_utils.py novoindex assembly.fasta --loglevel=DEBUG
    read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
    read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

    read_utils.py align_and_fix \
      ${reads_unmapped_bam} \
      assembly.fasta \
      --outBamAll ${sample_name}.bam \
      --outBamFiltered ${sample_name}.mapped.bam \
      --GATK_PATH gatk/ \
      --aligner ${aligner} \
      --aligner_options "${aligner_options}" \
      --JVMmemory=3g \
      --loglevel=DEBUG

    samtools index ${sample_name}.mapped.bam

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
  }

  output {
    File reads_bam = "${sample_name}.bam"
    File reads_bam_flagstat = "${sample_name}.bam.flagstat.txt"
    File mapped_reads_bam = "${sample_name}.mapped.bam"
    File coverage_plot = "${sample_name}.coverage_plot.pdf"
    Int assembly_length = read_int("assembly_length")
    Int assembly_length_unambiguous = read_int("assembly_length_unambiguous")
    Int reads_aligned = read_int("reads_aligned")
    Int read_pairs_aligned = read_int("read_pairs_aligned")
    Int bases_aligned = read_int("bases_aligned")
    Int mean_coverage = read_int("mean_coverage")
  }

  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "3500 MB"
    cpu: 4
  }
}


task fastqc {
  File reads_bam

  String reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v reports.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi
    ln -s ${reads_bam} ${reads_basename}.bam
    fastqc -t `nproc` ${reads_basename}.bam
  }

  output {
    File fastqc_html = "${reads_basename}_fastqc.html"
    File fastqc_zip  = "${reads_basename}_fastqc.zip"
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "broadinstitute/viral-ngs"
  }
}


task spikein_report {
  File  reads_bam
  File? spikein_db = "gs://sabeti-public-dbs-gz/spikeins/ercc_spike-ins.fasta"
  Int?  minScoreToFilter = 60

  String reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    # for those backends that prefer to override our Docker ENTRYPOINT
    if [ -z "$(command -v read_utils.py)" ]; then
      source /opt/viral-ngs/source/docker/container_environment.sh
    fi
    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py bwamem_idxstats \
      ${reads_basename}.bam \
      ${spikein_db} \
      --outStats ${reads_basename}.spike_count.txt \
      --minScoreToFilter=${minScoreToFilter} \
      --loglevel=DEBUG
  }

  output {
    File report = "${reads_basename}.spike_count.txt"
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "broadinstitute/viral-ngs"
  }
}


task plot_coverage {
  # TO DO: add a BWA option
  # TO DO: make GATK indel-realigner optional
  String sample_name

  File assembly_fasta
  File reads_unmapped_bam

  File gatk_jar
  File? novocraft_license

  String? aligner="novoalign" # novoalign or bwa
  String? aligner_options="-r Random -l 40 -g 40 -x 20 -t 100 -k"

  Boolean? skip_mark_dupes=false
  Boolean? plot_only_non_duplicates=false

  command {
    set -ex -o pipefail

    # prep GATK
    mkdir gatk
    if [[ ${gatk_jar} == *.tar.bz2 ]]; then
      tar -xjvof ${gatk_jar} -C gatk
    else
      ln -s ${gatk_jar} gatk/GenomeAnalysisTK.jar
    fi

    cp ${assembly_fasta} assembly.fasta
    read_utils.py novoindex assembly.fasta --loglevel=DEBUG
    read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
    read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

    read_utils.py align_and_fix \
      ${reads_unmapped_bam} \
      assembly.fasta \
      --outBamAll ${sample_name}.all.bam \
      --outBamFiltered ${sample_name}.mapped.bam \
      --GATK_PATH gatk/ \
      --aligner ${aligner} \
      --aligner_options "${aligner_options}" \
      ${true='--skipMarkDupes' false="" skip_mark_dupes} \
      --JVMmemory=3g \
      --loglevel=DEBUG

    samtools index ${sample_name}.mapped.bam

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
    grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    #echo $(( $(cat bases_aligned) / $(cat assembly_length) )) | tee mean_coverage
    python -c "print (float("`cat bases_aligned`")/"`cat assembly_length`") if "`cat assembly_length`">0 else 0" > mean_coverage

    # fastqc mapped bam
    reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html

    PLOT_DUPE_OPTION=""
    if [[ "${skip_mark_dupes}" != "true" ]]; then
      PLOT_DUPE_OPTION="${true='--plotOnlyNonDuplicates' false="" plot_only_non_duplicates}"
    fi

    # plot coverage
    if [ $(cat reads_aligned) != 0 ]; then
      reports.py plot_coverage \
        ${sample_name}.mapped.bam \
        ${sample_name}.coverage_plot.pdf \
        --plotFormat pdf \
        --plotWidth 1100 \
        --plotHeight 850 \
        --plotDPI 100 \
        $PLOT_DUPE_OPTION \
        --plotTitle "${sample_name} coverage plot" \
        --loglevel=DEBUG
    else
      touch ${sample_name}.coverage_plot.pdf
    fi
  }

  output {
    File   aligned_bam                 = "${sample_name}.all.bam"
    File   aligned_bam_idx             = "${sample_name}.all.bai"
    File   aligned_bam_flagstat        = "${sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam      = "${sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx  = "${sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc   = "${sample_name}.mapped_fastqc.html"
    File   coverage_plot               = "${sample_name}.coverage_plot.pdf"
    Int    assembly_length             = read_int("assembly_length")
    Int    assembly_length_unambiguous = read_int("assembly_length_unambiguous")
    Int    reads_aligned               = read_int("reads_aligned")
    Int    read_pairs_aligned          = read_int("read_pairs_aligned")
    Int    bases_aligned               = read_int("bases_aligned")
    Float  mean_coverage               = read_float("mean_coverage")
    String viralngs_version            = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3500 MB"
    cpu: 4
    dx_instance_type: "mem1_ssd1_x8"
  }
}


task coverage_report {
  Array[File]+ mapped_bams
  Array[File]  mapped_bam_idx # optional.. speeds it up if you provide it, otherwise we auto-index
  String       out_report_name="coverage_report.txt"

  command {
    reports.py coverage_only \
      ${sep=' ' mapped_bams} \
      ${out_report_name} \
      --loglevel DEBUG
  }

  output {
    File  coverage_report = "${out_report_name}"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2000 MB"
    cpu: 2
    dx_instance_type: "mem1_ssd2_x4"
  }
}


task fastqc {
  File reads_bam

  String reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    reports.py fastqc ${reads_bam} ${reads_basename}_fastqc.html
  }

  output {
    File   fastqc_html      = "${reads_basename}_fastqc.html"
    String viralngs_version = "viral-ngs_version_unknown"
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "quay.io/broadinstitute/viral-ngs"
    dx_instance_type: "mem1_ssd1_x4"
  }
}


task spikein_report {
  File  reads_bam
  File  spikein_db
  Int?  minScoreToFilter = 60
  Int?  topNHits = 3

  String reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail

    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py bwamem_idxstats \
      ${reads_basename}.bam \
      ${spikein_db} \
      --outStats ${reads_basename}.spike_count.txt.unsorted \
      --minScoreToFilter=${minScoreToFilter} \
      --loglevel=DEBUG

      sort -b -r -n -k3 ${reads_basename}.spike_count.txt.unsorted > ${reads_basename}.spike_count.txt
      head -n ${topNHits} ${reads_basename}.spike_count.txt > ${reads_basename}.spike_count.top_${topNHits}_hits.txt
  }

  output {
    File   report           = "${reads_basename}.spike_count.txt"
    File   report_top_hits  = "${reads_basename}.spike_count.top_${topNHits}_hits.txt"
    String viralngs_version = "viral-ngs_version_unknown"
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "quay.io/broadinstitute/viral-ngs"
    dx_instance_type: "mem1_ssd1_x4"
  }
}

task spikein_summary {
  Array[File]+  spikein_count_txt

  command {
    set -ex -o pipefail

    mkdir spike_summaries
    cp ${sep=' ' spikein_count_txt} spike_summaries/

    reports.py aggregate_spike_count spike_summaries/ spikein_summary.tsv \
      --loglevel=DEBUG
  }

  output {
    File   spikein_summary  = "spikein_summary.tsv"
    String viralngs_version = "viral-ngs_version_unknown"
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "quay.io/broadinstitute/viral-ngs"
    dx_instance_type: "mem1_ssd1_x4"
  }
}


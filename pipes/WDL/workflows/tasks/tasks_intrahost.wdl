
task isnvs_per_sample {
  File mapped_bam
  File assembly_fasta

  Int? threads
  Int? minReadsPerStrand
  Int? maxBias

  String sample_name = basename(basename(basename(mapped_bam, ".bam"), ".all"), ".mapped")

  command {
    intrahost.py vphaser_one_sample \
        ${mapped_bam} \
        ${assembly_fasta} \
        vphaser2.${sample_name}.txt.gz \
        ${'--vphaserNumThreads' + threads} \
        --removeDoublyMappedReads \
        ${'--minReadsEach' + minReadsPerStrand} \
        ${'--maxBias' + maxBias}
  }

  output {
    File   isnvsFile        = "vphaser2.${sample_name}.txt.gz"
    String viralngs_version = "viral-ngs_version_unknown"
  }
  runtime {
    memory: "7 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}


task isnvs_vcf {
  Array[File] vphaser2Calls # vphaser output; ex. vphaser2.${sample}.txt.gz
  Array[File] perSegmentMultiAlignments # aligned_##.fasta, where ## is segment number
  File reference_fasta

  Array[String]? snpEffRef # list of accessions to build/find snpEff database
  Array[String]? sampleNames # list of sample names
  String? emailAddress # email address passed to NCBI if we need to download reference sequences
  Boolean naiveFilter=false

  command {
    set -ex -o pipefail

    SAMPLES="${sep=' ' sampleNames}"
    if [ -n "$SAMPLES" ]; then SAMPLES="--samples $SAMPLES"; fi

    providedSnpRefAccessions="${sep=' ' snpEffRef}"
    if [ -n "$providedSnpRefAccessions" ]; then 
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('${reference_fasta}', 'fasta'))))")"
    fi

    echo "snpRefAccessions: $snpRefAccessions"

    intrahost.py merge_to_vcf \
        ${reference_fasta} \
        isnvs.vcf.gz \
        $SAMPLES \
        --isnvs ${sep=' ' vphaser2Calls} \
        --alignments ${sep=' ' perSegmentMultiAlignments} \
        --strip_chr_version \
        ${true="--naive_filter" false="" naiveFilter} \
        --parse_accession
        
    interhost.py snpEff \
        isnvs.vcf.gz \
        $snpRefAccessions \
        isnvs.annot.vcf.gz \
        ${'--emailAddress=' + emailAddress}

    intrahost.py iSNV_table \
        isnvs.annot.vcf.gz \
        isnvs.annot.txt.gz
  }

  output {
    Array[File] isnvFiles        = ["isnvs.vcf.gz", "isnvs.vcf.gz.tbi", "isnvs.annot.vcf.gz", "isnvs.annot.txt.gz", "isnvs.annot.vcf.gz.tbi"]
    String      viralngs_version = "viral-ngs_version_unknown"
  }
  runtime {
    memory: "4 GB"
    docker: "quay.io/broadinstitute/viral-ngs"
  }
}



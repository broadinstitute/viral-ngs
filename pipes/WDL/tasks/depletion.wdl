# ======================================================================
# deplete: 
#   Runs a full human read depletion pipeline and removes PCR duplicates
# ======================================================================
task deplete {
  String sample

  File inputBam
  Array[File] bmtaggerDbPaths
  Array[File] blastDbPaths

  Int? threads = 8
  Int? JVMmemory = 15
  
  command {
    taxon_filter.py deplete_human \
      "${inputBam}" \
      "tmpfile-${sample}.raw.bam" \
      "tmpfile-${sample}.bmtagger_depleted.bam" \
      "tmpfile-${sample}.rmdup.bam" \
      "${sample}.cleaned.bam" \
      --bmtaggerDbs "${sep=' ' bmtaggerDbPaths+}" \
      --blastDbs "${sep=' ' blastDbPaths+}" \
      --chunkSize=0 \
      "${'--threads' + threads}" \
      "${'--JVMmemory' + JVMmemory + 'g'}"
  }
  output {
    File cleaned          = "${sample}.cleaned.bam"
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "15GB"
    cpu: "8"
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}


# ======================================================================
# filterToTaxon: 
#   This step reduces the read set to a specific taxon (usually the genus
#   level or greater for the virus of interest)
# ======================================================================
task filterToTaxon {
  String sample

  File inputBam
  File lastalDbPath

  command {
    taxon_filter.py filter_lastal_bam \
      "${inputBam}" \
      "${lastalDbPath}" \
      "${sample}.taxfilt.bam"
  }

  output {
    File taxfiltBam = "${sample}.taxfilt.bam"
  }
  runtime {
    docker: "broadinstitute/viral-ngs"
    memory: "7GB"
    cpu: "8"
    zones: "us-east1-b us-east1-c us-east1-d"
    disks: "local-disk 375 LOCAL"
  }
}


# ======================================================================
# merge_one_per_sample:
#   All of the above depletion steps are run once per flowcell per lane per
#   multiplexed sample.  This reduces recomputation on the same data when
#   additional sequencing runs are performed on the same samples.  This
#   step merges reads down to one-per-sample, which all downstream
#   analysis steps require.  For cleaned and taxfilt outputs, we re-run
#   the rmdup step on a per-library basis after merging.
# ======================================================================

task merge_one_per_sample {
  String sample
  String adjective

  Array[File] inputBams

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      "${sample}.${adjective}.bam" \
      --picardOptions SORT_ORDER=queryname
  }

  output {
    File mergedBam = "${sample}.${adjective}.bam"
  }

  runtime{
    memory: "10GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task merge_one_per_sample_rmdup {
  String sample
  String adjective

  Array[File] inputBams

  Int? JVMmemory

  command {
    read_utils.py merge_bams \
      "${sep=' ' inputBams+}" \
      "temp_merged_${sample}.bam" \
      --picardOptions SORT_ORDER=queryname;
    read_utils.py rmdup_mvicuna_bam \
      "temp_merged_${sample}.bam" \
      "${sample}.${adjective}.bam" \
      "${'--JVMmemory' + JVMmemory || '8g'}"
  }

  output {
    File mergedBam = "${sample}.${adjective}.bam"
  }

  runtime{
    memory: "10GB"
    docker: "broadinstitute/viral-ngs"
  }
}
task download_reference_genome {
  String referenceName
  Array[String] accessions # NCBI accessions to include in the reference
  String emailAddress

  command {
    ncbi.py fetch_fastas \
        "${emailAddress}" \
        "./" \
        "${sep=' ' accessions}" \
        --combinedFilePrefix "${referenceName}" \
        --removeSeparateFiles \
        --forceOverwrite
    ncbi.py fetch_feature_tables \
        "${emailAddress}" \
        "./" \
        "${sep=' ' accessions}" \
        --forceOverwrite
  }

  output {
    File referenceFasta = "${referenceName}.fasta"
    Array[File] featureTables = glob("*.tbl")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task download_lastal_sources {
  String referenceName
  Array[String] accessions # NCBI accessions to include in the lastal db
  String emailAddress

  command {
    ncbi.py fetch_fastas \
        "${emailAddress}" \
        "./" \
        "${sep=' ' accessions}" \
        --combinedFilePrefix lastal \
        --removeSeparateFiles \
        --forceOverwrite \
        --chunkSize 300
  }

  output {
    File referenceFasta = "lastal.fasta"
    Array[File] featureTables = glob("*.tbl")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task build_lastal_db {
  File    sequences_fasta

  String  db_name = basename(sequences_fasta, ".fasta")

  command {
    set -ex -o pipefail
    taxon_filter.py lastal_build_db ${sequences_fasta} ./ --loglevel=DEBUG
    tar -c ${db_name}* | lz4 -9 > ${db_name}.tar.lz4
  }

  output {
    File lastal_db = "${db_name}.tar.lz4"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x4"
  }
}

task download_annotation {
  String        referenceName
  Array[String] accessions
  String        emailAddress

  command {
    set -ex -o pipefail
    ncbi.py fetch_feature_tables \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --combinedFilePrefix ${referenceName} \
        --loglevel DEBUG
  }

  output {
    File        referenceFasta  = "${referenceName}.fasta"
    File        featureTable    = "${referenceName}.tbl"
    Array[File] featureTables   = glob("*.tbl")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task annot_transfer {
  # TODO: Iterate over chr-specifc MSAs in workflow rather than in task

  File chrMultipleAlignment # fasta; multiple alignments of sample sequences
  File referenceFeatureTable # feature table corresponding to the chr in the alignment
  File referenceGenome # fasta

  command {
    ncbi.py tbl_transfer_prealigned \
        "${chrMultipleAlignment}" \
        "${referenceGenome}" \
        "${referenceFeatureTable}" \
        "./" \
        --oob_clip
  }

  output {
    Array[File] featureTables = glob(".tbl")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task prepare_genbank {
  Array[File] fastaFiles
  File assemblySummary # summary.assembly.txt
  File featureTableDir

  String genbankTemplate
  String genbankSourceTable
  String biosampleMap
  String sequencingTech
  String comment

  command {
    ncbi.py prep_genbank_files \
        "${genbankTemplate}" \
        "${sep=' ' fastaFiles}" \
        "${featureTableDir}" \
        --master_source_table "${genbankSourceTable}" \
        --sequencing_tech "${sequencingTech}" \
        --biosample_map "${biosampleMap}" \
        --coverage_table "${assemblySummary}" \
        --comment "${comment}"
  }

  output {
    File errorSummary = "${featureTableDir}/errorsummary.val"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}
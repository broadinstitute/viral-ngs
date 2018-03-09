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
    ncbi.py fetch_feature_tables \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --combinedFilePrefix ${referenceName} \
        --loglevel DEBUG
  }

  output {
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
  File chr_mutli_aln_fasta # fasta; multiple alignments of sample sequences
  File reference_fasta # fasta
  File reference_feature_table # feature table corresponding to the chr in the alignment

  command {
    ncbi.py tbl_transfer_prealigned \
        ${chr_mutli_aln_fasta} \
        ${reference_fasta} \
        ${reference_feature_table} \
        . \
        --oob_clip \
        --loglevel DEBUG
  }

  output {
    Array[File] featureTables = glob("*.tbl")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task prepare_genbank {
  Array[File]+ assemblies_fasta
  Array[File]+ annotations_tbl
  File         authors_sbt
  File         assemblySummary # summary.assembly.txt
  File         genbankSourceTable
  File         biosampleMap
  String       sequencingTech
  String       comment

  command {
    set -ex -o pipefail
    cp ${sep=' ' annotations_tbl} .
    ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        --master_source_table ${genbankSourceTable} \
        --sequencing_tech ${sequencingTech} \
        --biosample_map ${biosampleMap} \
        --coverage_table ${assemblySummary} \
        --comment ${comment} \
        --loglevel DEBUG
    tar -czpvf ncbi_package.tar.gz *.val *.cmt *.fsa *.gbf *.sqn *.src *.tbl
  }

  output {
    Array[File] sequin_files = glob("*.sqn")
    File        ncbi_package = "ncbi_package.tar.gz"
    File        errorSummary = "errorsummary.val"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}
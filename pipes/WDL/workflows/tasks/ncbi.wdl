
task download_fasta {
  String         out_prefix
  Array[String]+ accessions
  String         emailAddress

  command {
    ncbi.py fetch_fastas \
        ${emailAddress} \
        . \
        ${sep=' ' accessions} \
        --combinedFilePrefix ${out_prefix} \
  }

  output {
    File sequences_fasta = "${out_prefix}.fasta"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task download_annotations {
  Array[String]+ accessions
  String         emailAddress
  String         combined_out_prefix

  command {
    set -ex -o pipefail
    ncbi.py fetch_feature_tables \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --loglevel DEBUG
    ncbi.py fetch_fastas \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --combinedFilePrefix "${combined_out_prefix}" \
        --loglevel DEBUG
  }

  output {
    File        combined_fasta = "${combined_out_prefix}.fasta"
    Array[File] genomes_fasta  = glob("*.fasta")
    Array[File] features_tbl   = glob("*.tbl")
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

task annot_transfer {
  File chr_mutli_aln_fasta # fasta; multiple alignments of sample sequences for a single chr
  File reference_fasta # fasta (may contain multiple chrs, only one with the same name as reference_feature_table will be used)
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
    Array[File] transferred_feature_tables = glob("*.tbl")
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
  File?        coverage_table # summary.assembly.txt
  File?        genbankSourceTable
  File?        biosampleMap
  String?      sequencingTech
  String?      comment
  String       out_prefix = "ncbi_package"

  command {
    set -ex -o pipefail
    cp ${sep=' ' annotations_tbl} .
    ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        ${'--master_source_table=' + genbankSourceTable} \
        ${'--sequencing_tech=' + sequencingTech} \
        ${'--biosample_map=' + biosampleMap} \
        ${'--coverage_table=' + coverage_table} \
        ${'--comment=' + comment} \
        --loglevel DEBUG
    tar -czpvf ${out_prefix}.tar.gz *.val *.cmt *.fsa *.gbf *.sqn *.src *.tbl
  }

  output {
    Array[File] sequin_files = glob("*.sqn")
    File        ncbi_package = "${out_prefix}.tar.gz"
    File        errorSummary = "errorsummary.val"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_x2"
  }
}

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
    File   sequences_fasta  = "${out_prefix}.fasta"
    String viralngs_version = "viral-ngs_version_unknown"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
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
    File        combined_fasta   = "${combined_out_prefix}.fasta"
    Array[File] genomes_fasta    = glob("*.fasta")
    Array[File] features_tbl     = glob("*.tbl")
    String      viralngs_version = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task annot_transfer {
  Array[File]+ multi_aln_fasta         # fasta; multiple alignments of sample sequences for each chromosome
  File         reference_fasta         # fasta; all chromosomes in one file
  Array[File]+ reference_feature_table # tbl; feature table corresponding to each chromosome in the alignment

  Array[Int]   chr_nums=range(length(multi_aln_fasta))

  command {
    set -ex -o pipefail
    echo ${sep=' ' multi_aln_fasta} > alignments.txt
    echo ${sep=' ' reference_feature_table} > tbls.txt
    for i in ${sep=' ' chr_nums}; do
      _alignment_fasta=`cat alignments.txt | cut -f $(($i+1)) -d ' '`
      _feature_tbl=`cat tbls.txt | cut -f $(($i+1)) -d ' '`
      ncbi.py tbl_transfer_prealigned \
          $_alignment_fasta \
          ${reference_fasta} \
          $_feature_tbl \
          . \
          --oob_clip \
          --loglevel DEBUG
    done
  }

  output {
    Array[File]+ transferred_feature_tables = glob("*.tbl")
    String       viralngs_version           = "viral-ngs_version_unknown"
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task prepare_genbank {
  Array[File]+ assemblies_fasta
  Array[File]+ annotations_tbl
  File         authors_sbt
  File         biosampleMap
  File         genbankSourceTable
  File?        coverage_table # summary.assembly.txt (from Snakemake) -- change this to accept a list of mapped bam files and we can create this table ourselves
  String       sequencingTech
  String       comment # TO DO: make this optional
  String       organism
  String       molType = "cRNA"

  command {
    set -ex -o pipefail
    cp ${sep=' ' annotations_tbl} .
    ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        --mol_type ${molType} \
        --organism "${organism}" \
        --biosample_map ${biosampleMap} \
        --master_source_table ${genbankSourceTable} \
        ${'--coverage_table ' + coverage_table} \
        --comment "${comment}" \
        --sequencing_tech "${sequencingTech}" \
        --loglevel DEBUG
    mv errorsummary.val errorsummary.val.txt # to keep it separate from the glob
  }

  output {
    Array[File] sequin_files             = glob("*.sqn")
    Array[File] structured_comment_files = glob("*.cmt")
    Array[File] genbank_preview_files    = glob("*.gbf")
    Array[File] source_table_files       = glob("*.src")
    Array[File] fasta_per_chr_files      = glob("*.fsa")
    Array[File] validation_files         = glob("*.val")
    File        errorSummary             = "errorsummary.val.txt"
    String      viralngs_version         = "viral-ngs_version_unknown"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}



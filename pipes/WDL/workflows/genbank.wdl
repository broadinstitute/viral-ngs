import "interhost.wdl" as interhost
import "ncbi.wdl" as ncbi

workflow genbank {

  File          reference_fasta
  Array[File]+  assemblies_fasta     # one per genome
  Array[File]+  ref_annotations_tbl  # one per chromosome

  call interhost.multi_align_mafft_ref as mafft {
    input:
      reference_fasta = reference_fasta,
      assemblies_fasta = assemblies_fasta
  }

  call ncbi.annot_transfer as annot {
    input:
      multi_aln_fasta = mafft.alignments_by_chr,
      reference_fasta = reference_fasta,
      reference_feature_table = ref_annotations_tbl
  }
 
  call ncbi.prepare_genbank as prep_genbank {
    input:
      assemblies_fasta = assemblies_fasta,
      annotations_tbl = annot.transferred_feature_tables
  }

}

import "tasks/interhost.wdl" as interhost
import "tasks/ncbi.wdl" as ncbi

workflow align_and_annot {

  File          reference_fasta
  Array[File]+  assemblies_fasta
  Array[File]+  annotations_tbl

  call interhost.multi_align_mafft_ref as mafft {
    input:
      reference_fasta = reference_fasta,
      assemblies_fasta = assemblies_fasta
  }

  scatter(aln_by_chr, ref_annot_by_chr in zip(mafft.alignments_by_chr, annotations_tbl)) {
    call ncbi.annot_transfer as annot {
      input:
        chr_mutli_aln_fasta = aln_by_chr,
        reference_fasta = reference_fasta,
        reference_feature_table = ref_annot_by_chr
    }
    call ncbi.prepare_genbank as genbank {
      input:
        assemblies_fasta = assemblies_fasta,
        annotations_tbl = annot.transferred_feature_tables,
        out_prefix = basename(annotations_tbl, '.tbl')
    }
  }
}

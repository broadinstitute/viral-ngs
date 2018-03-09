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

  scatter(chr_num in range(len(mafft.alignments_by_chr))) {
    call ncbi.annot_transfer as annot {
      input:
        chr_mutli_aln_fasta = mafft.alignments_by_chr[chr_num],
        reference_fasta = reference_fasta,
        reference_feature_table = annotations_tbl[chr_num]
    }
    call ncbi.prepare_genbank as genbank {
      input:
        assemblies_fasta = assemblies_fasta,
        annotations_tbl = annot.featureTables # I'm worried that the order got messed up and we'll have to remap?
    }
  }
}

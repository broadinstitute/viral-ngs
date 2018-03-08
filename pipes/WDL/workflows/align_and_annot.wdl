import "tasks/interhost.wdl" as interhost
import "tasks/ncbi.wdl" as ncbi
import "tasks/reports.wdl" as reports


workflow align_and_annot {

  File        reference_fasta
  Array[File] assemblies_fasta

  call interhost.multi_align_mafft_ref as mafft {
    input:
      reference_fasta = reference_fasta,
      assemblies_fasta = assemblies_fasta
  }

  scatter(mutli_aln_fasta in mafft.alignments_by_chr) {
    call ncbi.annot_transfer as annot {
      input:
        chr_mutli_aln_fasta = mutli_aln_fasta,
        reference_fasta = reference_fasta
    }
    call ncbi.prepare_genbank as genbank {
      input:
        assemblies_fasta = assemblies_fasta,
        annotations_tbl = annot.featureTables # I'm worried that the order got messed up
    }
  }
}

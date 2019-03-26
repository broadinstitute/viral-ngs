import "tasks_interhost.wdl" as interhost
import "tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_merge_to_vcf {
    File          reference_fasta
    Array[File]+  assemblies_fasta     # one per genome

    call interhost.multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    call tasks_intrahost.isnvs_vcf {
        input:
            perSegmentMultiAlignments = mafft.alignments_by_chr,
            reference_fasta = reference_fasta
    }

}

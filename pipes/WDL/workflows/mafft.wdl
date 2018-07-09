import "interhost.wdl" as interhost

workflow mafft {
    Array[File]+  assemblies_fasta     # one per genome
    call interhost.multi_align_mafft {
        input:
            assemblies_fasta = assemblies_fasta
    }
}

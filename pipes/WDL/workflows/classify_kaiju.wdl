import "tasks_metagenomics.wdl" as metagenomics

workflow classify_kaiju {
    call metagenomics.kaiju
}

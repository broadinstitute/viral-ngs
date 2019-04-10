import "tasks_metagenomics.wdl" as metagenomics

workflow classify_krakenuniq {
    call metagenomics.krakenuniq
}

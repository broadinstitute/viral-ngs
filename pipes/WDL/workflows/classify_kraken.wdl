import "tasks_metagenomics.wdl" as metagenomics

workflow classify_kraken {
    call metagenomics.krakenuniq
}

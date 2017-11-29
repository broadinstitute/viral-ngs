import "tasks/metagenomics.wdl" as metagenomics

workflow classify_kraken {
  call metagenomics.kraken
}

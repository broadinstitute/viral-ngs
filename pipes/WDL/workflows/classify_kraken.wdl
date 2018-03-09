import "metagenomics.wdl" as metagenomics

workflow classify_kraken {
  call metagenomics.kraken
}

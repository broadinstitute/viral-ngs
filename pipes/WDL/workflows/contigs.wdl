import "metagenomics.wdl" as metagenomics
import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly

workflow contigs {

  call taxon_filter.deplete_taxa as deplete

  call assembly.assemble as spades {
    input:
      assembler = "spades",
      reads_unmapped_bam = deplete.cleaned_bam
  }

  # TO DO: taxonomic classification of contigs

}

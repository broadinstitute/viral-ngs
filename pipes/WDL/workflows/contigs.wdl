import "tasks/metagenomics.wdl" as metagenomics
import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/assembly.wdl" as assembly

workflow contigs {

  call taxon_filter.deplete_taxa as deplete

  call assembly.assemble as spades {
    input:
      assembler = "spades",
      reads_unmapped_bam = deplete.cleaned_bam
  }

  # TO DO: taxonomic classification of contigs

}

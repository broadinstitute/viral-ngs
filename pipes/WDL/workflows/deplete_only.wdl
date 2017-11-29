import "tasks/taxon_filter.wdl" as taxon_filter

workflow deplete_only {
  call taxon_filter.deplete_taxa
}
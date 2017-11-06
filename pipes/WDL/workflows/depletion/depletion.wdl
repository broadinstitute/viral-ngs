import "../../tasks/depletion.wdl" as depletion

workflow all_deplete {
  Array[File] bam_input_files

  scatter(i in bam_input_files) {
    call depletion.deplete as deplete{input: inputBam=i} 
  }

  scatter(i in deplete.cleaned) {
    call depletion.filterToTaxon as filter{input: inputBam=i}
  }

  output {
    filter.taxfiltBam
  }
}
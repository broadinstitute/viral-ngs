import "../../tasks/demux.wdl" as demux

workflow all_demux {
  call demux.illumina_demux as illumina_demux

  output {
    illumina_demux.outputBams
  }
}
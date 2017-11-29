import "tasks/demux.wdl" as tasks_demux

workflow demux_only {
  call tasks_demux.illumina_demux
}
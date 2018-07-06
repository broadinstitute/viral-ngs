import "demux.wdl" as tasks_demux

workflow demux_only {
    File flowcell_tgz
    call tasks_demux.illumina_demux {
        input:
            flowcell_tgz = flowcell_tgz
  }
}

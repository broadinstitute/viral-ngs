import "tasks/demux.wdl" as tasks_demux

workflow demux_only {
  call tasks_demux.illumina_demux as illumina_demux

  output {
    Array[File] raw_reads_unaligned_bams = illumina_demux.raw_reads_unaligned_bams
  }
}
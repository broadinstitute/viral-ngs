import "demux.wdl" as demux

workflow merge_bams {
    Array[File]+  in_bams
    String        out_basename
    call demux.merge_and_reheader_bams {
        input:
            in_bams = in_bams,
            out_basename = out_basename
    }
}

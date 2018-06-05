import "read_utils.wdl" as read_utils

workflow downsample_bams {
  call read_utils.downsample_bams as downsample_bams

}

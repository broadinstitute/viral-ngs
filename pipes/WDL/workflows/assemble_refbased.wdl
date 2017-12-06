import "tasks/assembly.wdl" as assembly

workflow assemble_refbased {
  call assembly.refine_2x_and_plot
}
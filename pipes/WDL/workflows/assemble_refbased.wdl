import "tasks_assembly.wdl" as assembly
import "tasks_reports.wdl" as reports

workflow assemble_refbased {
  call assembly.refine_2x_and_plot
  call reports.software_version
}

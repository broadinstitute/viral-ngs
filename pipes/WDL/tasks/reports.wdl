task fastqc_report {
  String sample

  command {

  }

  output {

  }
  runtime {
    memory: "3GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task consolidate_fastqc_on_all_runs {
  String sample

  command {

  }

  output {

  }
  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}

task spikein_report {
  String sample

  command {

  }

  output {

  }
  runtime {
    memory: "3GB"
    docker: "broadinstitute/viral-ngs"
  }
}

task consolidate_spike_count {
  String sample

  command {

  }

  output {

  }
  runtime {
    docker: "broadinstitute/viral-ngs"
  }
}
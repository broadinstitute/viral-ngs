#DX_SKIP_WORKFLOW

import "tasks_demux.wdl" as demux
import "tasks_metagenomics.wdl" as metagenomics
import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_assembly.wdl" as assembly
import "tasks_reports.wdl" as reports

workflow demux_metag {
  call demux.illumina_demux as illumina_demux

  scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
    call reports.spikein_report as spikein {
      input:
        reads_bam = raw_reads
    }
    call taxon_filter.deplete_taxa as deplete {
      input:
        raw_reads_unmapped_bam = raw_reads
    }
    call assembly.assemble as spades {
      input:
        assembler = "spades",
        reads_unmapped_bam = deplete.cleaned_bam
    }
  }

  call reports.spikein_summary as spike_summary {
      input:
          spikein_count_txt = spikein.report
  }

  call metagenomics.krakenuniq as kraken {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams,
  }
  call reports.aggregate_metagenomics_reports as metag_summary_report {
      input:
          kraken_summary_reports = kraken.krakenuniq_summary_reports
  }
  call metagenomics.kaiju as kaiju {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams,
  }
}

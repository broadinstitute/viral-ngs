import "tasks/demux.wdl" as demux
import "tasks/metagenomics.wdl" as metagenomics
import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/assembly.wdl" as assembly
import "tasks/reports.wdl" as reports

workflow demux_plus {

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

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams
  }

}

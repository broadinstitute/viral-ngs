import "tasks/demux.wdl" as demux
import "tasks/metagenomics.wdl" as metagenomics
import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/reports.wdl" as reports

workflow demux_plus {

  call demux.illumina_demux as illumina_demux

  scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
    call reports.fastqc as fastqc {
      input:
        reads_bam = raw_reads
    }
    call reports.spikein_report as spikein {
      input:
        reads_bam = raw_reads
    }
    call taxon_filter.deplete_taxa as deplete {
      input:
        raw_reads_unmapped_bam = raw_reads
    }
  }

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams
  }

  scatter(classified_reads in kraken.kraken_classified_reads) {
    call metagenomics.krona as krona {
      input:
        classified_reads_txt_gz = classified_reads,
    }
  }

}

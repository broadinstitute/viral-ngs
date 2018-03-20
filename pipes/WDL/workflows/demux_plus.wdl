import "demux.wdl" as demux
import "metagenomics.wdl" as metagenomics
import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly
import "reports.wdl" as reports

workflow demux_plus {

	File? metadata_path

  call demux.illumina_demux as illumina_demux {
	   input:
		 		metadata_path = metadata_path
	}

  scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
    call reports.spikein_report as spikein {
      input:
        reads_bam = raw_reads
    }
    call taxon_filter.deplete_taxa as deplete {
      input:
        raw_reads_unmapped_bam = raw_reads,
				metadata_path = metadata_path
    }
    call assembly.assemble as spades {
      input:
        assembler = "spades",
        reads_unmapped_bam = deplete.cleaned_bam,
				metadata_path = metadata_path
    }
  }

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams,
			metadata_path = metadata_path
  }

}

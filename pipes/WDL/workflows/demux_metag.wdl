import "tasks/demux.wdl" as demux
import "tasks/metagenomics.wdl" as metagenomics
import "tasks/taxon_filter.wdl" as taxon_filter
import "tasks/assembly.wdl" as assembly
import "tasks/reports.wdl" as reports

workflow demux_metag {
  File krona_taxonomy_db_tgz

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
    call metagenomics.diamond_contigs as diamond {
      input:
        contigs_fasta = spades.contigs_fasta,
        reads_unmapped_bam = deplete.cleaned_bam,
        krona_taxonomy_db_tar_lz4 = krona_taxonomy_db_tgz
    }
  }

  call metagenomics.kraken as kraken {
    input:
      reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams,
      krona_taxonomy_db_tgz = krona_taxonomy_db_tgz
  }

}

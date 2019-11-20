import "tasks_demux.wdl" as demux
import "tasks_metagenomics.wdl" as metagenomics
import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_assembly.wdl" as assembly
import "tasks_reports.wdl" as reports
import "tasks_read_utils.wdl" as reads

workflow demux_plus {

    call demux.illumina_demux as illumina_demux

    File spikein_db
    File trim_clip_db
    Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? bwaDbs

    scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
        # de-duplicate raw reads
        call reads.dedup_bam as dedup {
            input:
                in_bam = raw_reads        
        }
 
        # count spike-ins in the sample
        #   NB: the spike-in report is created from raw reads 
        #       that have NOT been de-duplicated
        call reports.spikein_report as spikein {
            input:
                reads_bam = raw_reads,
                spikein_db = spikein_db
        }

        # deplete human/host genomic reads
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = dedup.dedup_bam,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }

        # create de novo contigs from depleted reads via spaces
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = deplete.cleaned_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }

        # classify de-duplicated reads to taxa via krakenuniq
        call metagenomics.krakenuniq as krakenuniq {
            input:
                reads_unmapped_bam = dedup.dedup_bam
        }
    }

    # summarize spike-in reports from all samples
    call reports.spikein_summary as spike_summary {
        input:
            spikein_count_txt = spikein.report
    }

    # summarize kraken reports from all samples
    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}

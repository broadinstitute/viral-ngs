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
        call reads.dedup_bam as dedup {
            input:
                in_bam = raw_reads        
        }
    }
    
    scatter(reads_bam in dedup.dedup_bam) {
        call reports.spikein_report as spikein {
            input:
                reads_bam = reads_bam,
                spikein_db = spikein_db
        }
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = reads_bam,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = deplete.cleaned_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }
    }

    call metagenomics.krakenuniq as krakenuniq {
        input:
            reads_unmapped_bam = dedup.dedup_bam
    }

    call reports.spikein_summary as spike_summary {
        input:
            spikein_count_txt = spikein.report
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}

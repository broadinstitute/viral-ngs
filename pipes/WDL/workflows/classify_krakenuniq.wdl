import "tasks_metagenomics.wdl" as metagenomics
import "tasks_reports.wdl" as reports
import "tasks_read_utils.wdl" as reads

workflow classify_krakenuniq {
    Array[File] unclassified_bams
    scatter(reads_bam in unclassified_bams) {
        call reads.dedup_bam as dedup {
            input:
                in_bam = reads_bam        
        }
    }
    call metagenomics.krakenuniq {
        input:
            reads_unmapped_bam = dedup.dedup_bam
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}

import "tasks_metagenomics.wdl" as metagenomics
import "tasks_reports.wdl" as reports

workflow classify_krakenuniq {
    call metagenomics.krakenuniq

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}

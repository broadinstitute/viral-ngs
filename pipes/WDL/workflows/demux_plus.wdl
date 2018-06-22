import "demux.wdl" as demux
import "metagenomics.wdl" as metagenomics
import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly
import "reports.wdl" as reports

workflow demux_plus {

    File flowcell_tgz
    call demux.illumina_demux as illumina_demux {
        input:
            flowcell_tgz = flowcell_tgz
    }

    File spikein_db
    File trim_clip_db
    Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? bwaDbs
    scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
        call reports.spikein_report as spikein {
            input:
                reads_bam = raw_reads,
                spikein_db = spikein_db
        }
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = raw_reads,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = deplete.cleaned_bam,
                trim_clip_db = trim_clip_db
        }
    }

    File kraken_db_tar_lz4
    File krona_taxonomy_db_tgz
    call metagenomics.kraken as kraken {
        input:
            reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams,
            kraken_db_tar_lz4 = kraken_db_tar_lz4,
            krona_taxonomy_db_tgz = krona_taxonomy_db_tgz
    }
}

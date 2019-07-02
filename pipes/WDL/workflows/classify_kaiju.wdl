import "tasks_metagenomics.wdl" as metagenomics
import "tasks_read_utils.wdl" as reads

workflow classify_kaiju {
    Array[File] unclassified_bams
    scatter(reads_bam in unclassified_bams) {
        call reads.dedup_bam as dedup {
            input:
                in_bam = reads_bam        
        }
    }
    call metagenomics.kaiju {
        input:
            reads_unmapped_bam = dedup.dedup_bam
    }
}

import "metagenomics.wdl" as metagenomics
import "taxon_filter.wdl" as taxon_filter
import "assembly.wdl" as assembly

workflow contigs {
    File raw_reads_unmapped_bam
    Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
    Array[File]? bwaDbs
    call taxon_filter.deplete_taxa as deplete {
        input:
            raw_reads_unmapped_bam = raw_reads_unmapped_bam,
            bmtaggerDbs = bmtaggerDbs,
            blastDbs = blastDbs,
            bwaDbs = bwaDbs
    }

    File trim_clip_db
    call assembly.assemble as spades {
        input:
            assembler = "spades",
            reads_unmapped_bam = deplete.cleaned_bam,
            trim_clip_db = trim_clip_db
    }

  # TO DO: taxonomic classification of contigs

}

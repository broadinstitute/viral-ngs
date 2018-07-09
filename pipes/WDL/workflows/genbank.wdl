import "interhost.wdl" as interhost
import "ncbi.wdl" as ncbi

workflow genbank {

    File          reference_fasta
    Array[File]+  assemblies_fasta     # one per genome
    call interhost.multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    Array[File]+  ref_annotations_tbl  # one per chromosome
    call ncbi.annot_transfer as annot {
        input:
            multi_aln_fasta = mafft.alignments_by_chr,
            reference_fasta = reference_fasta,
            reference_feature_table = ref_annotations_tbl
    }
 
    File         authors_sbt
    File         biosampleMap
    File         genbankSourceTable
    File?        coverage_table # summary.assembly.txt (from Snakemake) -- change this to accept a list of mapped bam files and we can create this table ourselves
    String       sequencingTech
    String       comment # TO DO: make this optional
    String       organism
    String       molType = "cRNA"
    call ncbi.prepare_genbank as prep_genbank {
        input:
            assemblies_fasta = assemblies_fasta,
            annotations_tbl = annot.transferred_feature_tables,
            authors_sbt = authors_sbt,
            biosampleMap = biosampleMap,
            genbankSourceTable = genbankSourceTable,
            coverage_table = coverage_table, # summary.assembly.txt (from Snakemake) -- change this to accept a list of mapped bam files and we can create this table ourselves
            sequencingTech = sequencingTech,
            comment = comment, # TO DO: make this optional
            organism = organism,
            molType = molType
    }
}

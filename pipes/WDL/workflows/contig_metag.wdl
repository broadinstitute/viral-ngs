import "tasks/metagenomics.wdl" as metagenomics
import "tasks/reports.wdl" as reports

workflow contig_metagenomics {
  File contigs_fasta
  File? blastn_nt_tgz
  File? blastx_nr_tgz
  File? blast_tax_db_tgz
  File? tax_db_tgz
  File? cdd_db_tgz
  File? rfam_db_tgz

  call metagenomics.download_cdd_db as download_cdd {
    input: cdd_db_tgz=cdd_db_tgz
  }

  call metagenomics.download_rfam_db as download_rfam {
    input: rfam_db_tgz=rfam_db_tgz
  }

  call metagenomics.download_blastn_nt_db as download_nt {
    input: nt_tgz=blastn_nt_tgz,
           tax_db_tgz=tax_db_tgz,
           blast_tax_db_tgz=blast_tax_db_tgz

  }

  call metagenomics.download_blastx_nr_db as download_nr {
    input: nr_tgz=blastx_nr_tgz,
           tax_db_tgz=tax_db_tgz,
           blast_tax_db_tgz=blast_tax_db_tgz
  }

  call metagenomics.rpsblast_models as rpsblast {
    input: contigs_fasta=contigs_fasta,
           cdd_db=download_cdd.db
  }

  call metagenomics.infernal_rfam as infernal_rfam {
    input: contigs_fasta=contigs_fasta,
           rfam_db_files=download_rfam.db_files
  }

  call metagenomics.blastn_contigs as blastn_contigs {
    input: contigs_fasta=contigs_fasta,
           blast_db_files=download_nt.db,
           tax_db_files=download_nt.tax_db
  }

  call metagenomics.blastx_contigs as blastx_contigs {
    input: contigs_fasta=contigs_fasta,
           blast_db_files=download_nr.db,
           tax_db_files=download_nr.tax_db
  }

  output {
    File infernal_tbl = infernal_rfam.infernal_tbl
    File cdd_report = rpsblast.cdd_report
    File blastn = blastn_contigs.blastn
    File blastn_lca = blastn_contigs.blastn_lca
    File blastn_report = blastn_contigs.blastn_report
    File blastx = blastx_contigs.blastx
    File blastx_lca = blastx_contigs.blastx_lca
    File blastx_report = blastx_contigs.blastx_report
  }
}

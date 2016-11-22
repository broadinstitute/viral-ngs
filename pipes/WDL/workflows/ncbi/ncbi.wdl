import "../../tasks/ncbi.wdl" as ncbi

workflow build_lastal_db {
  String emailAddress
  Array[String] accessions
  String referenceName

  call ncbi.download_lastal_sources as download_lastal_sources{input: accessions=accessions, emailAddress=emailAddress, referenceName=referenceName} 
  call ncbi.build_lastal_db as build_lastal_db{input: sequences=download_lastal_sources.referenceFasta}

  output {
    build_lastal_db.lastalDbFiles
  }
}
#!/usr/bin/python

# built-ins
import sys, os, logging

# third-party
from Bio import Entrez, SeqIO

log = logging.getLogger(__name__)

def _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite=False, rettype="fasta"):
    """ 
        This function downloads and saves files from NCBI nuccore.
    """
    db           = "nuccore"
    Entrez.email = emailAddress

    outEx = {
        "fasta": ".fasta",
        "ft":".tbl"
    }

    assert rettype in outEx.keys(), "The return type requested, %s, is not compatible with the nuccore fetch." % rettype

    outputDirectory = os.path.abspath(os.path.expanduser(destinationDir))

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    log.info( "Fetching %s entries from GenBank: %s\n" % (len(accessionList), ", ".join(accessionList[:10])))
    outputFiles = []
    for i,acc in enumerate(accessionList):
        outputFilePath = os.path.join(outputDirectory, acc+outEx[rettype])

        if not forceOverwrite:
            log.info("not overwriting, checking for existence")
            assert not os.path.exists(outputFilePath), """File %s already exists. Consider removing 
                this file or specifying a different output directory. The files for the accessions specified 
                can be overwritten if you add --forceOverwrite flag. Processing aborted.""" % outputFilePath

        #try:
        log.info("Fetching file %s: %s" % (i+1, acc))
        handle = Entrez.efetch(db=db, rettype=rettype, id=acc)

        with open(outputFilePath, "w") as outf:
            for line in handle:
                outf.write(line)
        outputFiles.append(outputFilePath)
        #except Exception as e:
        #    raise IOError( "Error! Could not fetch: %s\n %s" % (acc, e.message))  

    # return list of files
    return outputFiles

def fetch_fastas_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, rettype="fasta"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype)

def fetch_feature_tables_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, rettype="ft"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype)

if __name__ == "__main__":
    fastaFilePaths = fetch_fastas_from_genbank(["NC_004296.1", "NC_004297.1"], "~/Desktop/")
    for fastaFilePath in fastaFilePaths:
        print fastaFilePath
        fastaFile = SeqIO.parse(fastaFilePath, "fasta")
        for seq in fastaFile:
            print seq.id

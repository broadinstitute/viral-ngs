#!/usr/bin/python

# built-ins
import sys, os, logging

# third-party
from Bio import Entrez, SeqIO

log = logging.getLogger(__name__)

def get_feature_table_id(featureTableFile):
    seqid = ""
    with open(featureTableFile, 'rt') as inf:
        for line in inf:
            line = line.rstrip('\r\n')
            if not line:
                pass
            elif line.startswith('>'):
                # sequence identifier
                if not line.startswith('>Feature '):
                    raise Exception("not sure how to handle a non-Feature record")
                seqid = line[len('>Feature '):].strip()
                if not ((seqid.startswith('gb|') or seqid.startswith('ref|')) and seqid.endswith('|') and len(seqid)>4):
                    raise Exception("reference annotation does not refer to a GenBank or RefSeq accession")
                seqid = seqid[seqid.find("|")+1:-1]
    if len(seqid) > 0:
        return seqid
        
def _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite=False, rettype="fasta", combinedGenomeFilePrefix=None, removeSeparateFastas=False):
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

    if rettype == "fasta":
        # assert that we are not trying to remove the intermediate files without writing a combined file
        if removeSeparateFastas:
            assert combinedGenomeFilePrefix, """The intermediate FASTA files 
                can only be removed if a combined file is written via --combinedGenomeFilePrefix"""

        # build a path to the combined genome file
        if combinedGenomeFilePrefix:
            concatenatedGenomeFilepath = os.path.join(destinationDir, combinedGenomeFilePrefix+outEx[rettype])

            if not forceOverwrite:
                assert not os.path.exists(concatenatedGenomeFilepath), """File %s already exists. Consider removing 
                    this file or specifying a different output directory. The files for the accessions specified 
                    can be overwritten if you add --forceOverwrite flag. Processing aborted.""" % outputFilePath

            # concatenate the files together into one genome file
            with open(concatenatedGenomeFilepath, 'w') as outfile:
                for filePath in outputFiles:
                    with open(filePath) as infile:
                        for line in infile:
                            outfile.write(line)

            # if the option is specified, remove the intermediate fasta files
            if removeSeparateFastas:
                while len(outputFiles) > 0:
                    os.unlink(outputFiles.pop())

            # add the combined file to the list of files returned
            outputFiles.append(concatenatedGenomeFilepath)


    # return list of files
    return outputFiles

def fetch_fastas_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, combinedGenomeFilePrefix, removeSeparateFastas, rettype="fasta"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype, combinedGenomeFilePrefix, removeSeparateFastas)

def fetch_feature_tables_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, combinedGenomeFilePrefix, removeSeparateFastas, rettype="ft"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype, combinedGenomeFilePrefix, removeSeparateFastas)

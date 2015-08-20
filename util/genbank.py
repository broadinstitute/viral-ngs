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
        
def _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite=False, rettype="fasta", retmode="text", fileExt=None, combinedFilePrefix=None, removeSeparateFiles=False):
    """ 
        This function downloads and saves files from NCBI nuccore.
    """
    db           = "nuccore"
    Entrez.email = emailAddress

    outEx = {
        "fasta": "fasta",
        "ft":"tbl",
        "gb":"gbk"
    }

    assert rettype in outEx.keys(), "The return type requested, %s, is not compatible with the nuccore fetch." % rettype

    outputDirectory = os.path.abspath(os.path.expanduser(destinationDir))

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    # if the file extension to use is specified as the fileExt arg, use it
    # otherwise, use the default for the rettype according to the outEx dict, 
    # falling back to the retmode if there is no match
    if not fileExt:
        outputExtension = outEx.get(rettype, retmode) 
    else:
        outputExtension = str(fileExt)

    # ensure the extension starts with a ".", also allowing for passed-in 
    # extensions that already have it
    if outputExtension[:1] != ".":
        outputExtension = "." + outputExtension

    log.info( "Fetching %s entries from GenBank: %s\n" % (len(accessionList), ", ".join(accessionList[:10])))
    outputFiles = []
    for i,acc in enumerate(accessionList):

        outputFilePath = os.path.join(outputDirectory, acc+outputExtension)

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

    #if rettype == "fasta":
    # assert that we are not trying to remove the intermediate files without writing a combined file
    if removeSeparateFiles:
        assert combinedFilePrefix, """The intermediate files 
            can only be removed if a combined file is written via --combinedFilePrefix"""

    # build a path to the combined genome file
    if combinedFilePrefix:
        concatenatedGenomeFilepath = os.path.join(outputDirectory, combinedFilePrefix+outputExtension)

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
        if removeSeparateFiles:
            while len(outputFiles) > 0:
                os.unlink(outputFiles.pop())

        # add the combined file to the list of files returned
        outputFiles.append(concatenatedGenomeFilepath)


    # return list of files
    return outputFiles

def fetch_fastas_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, removeSeparateFiles, fileExt=None, rettype="fasta", retmode="text"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix, removeSeparateFiles)

def fetch_feature_tables_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, removeSeparateFiles, fileExt=None, rettype="ft", retmode="text"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix, removeSeparateFiles)

def fetch_full_records_from_genbank(accessionList, destinationDir, emailAddress, forceOverwrite, combinedFilePrefix, removeSeparateFiles, fileExt=None, rettype="gb", retmode="text"):
    return _fetch_from_nuccore(accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix, removeSeparateFiles)

if __name__ == "__main__":
    fetch_full_records_from_genbank(["JQ610675.1","JQ610676.1","JQ610677.1","JQ610678.1","JQ610679.1","JQ610680.1","JQ610681.1","JQ610682.1","JQ610683.1","JQ610684.1"], "~/Desktop/", "tomkinsc@broadinstitute.org", forceOverwrite=True, combinedFilePrefix="orungo", removeSeparateFiles=False)



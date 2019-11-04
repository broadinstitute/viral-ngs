#!/usr/bin/python

# built-ins
import time
import os
import re
import logging

# third-party
from Bio import Entrez

log = logging.getLogger(__name__)


def parse_accession_str(chr_ref):
    '''
        This tries to match an NCBI accession as defined here:
            http://www.ncbi.nlm.nih.gov/Sequin/acc.html
    '''
    m = re.search(r"(?P<accession>(?:[a-zA-Z]{1,6}|NC_)\d{1,10})(?:\.(?P<version>\d+))?.*", chr_ref)
    if m:
        chr_ref = m.group("accession")
    return chr_ref

def get_feature_table_id(featureTableFile):
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
                if not (
                    (seqid.startswith('gb|') or seqid.startswith('ref|'))):
                    raise Exception("reference annotation ID does not appear to refer to a GenBank or RefSeq accession: %s" % seqid)
                m = re.search(r"(?P<db>(?:gb|ref|dbj))\|(?:(?P<accession>[a-zA-Z0-9\._]+))+.*", seqid)
                if m:
                   seqid = m.group("accession")
                else:
                    raise Exception("Error parsing accession: %s" % seqid)
                return seqid

def _seq_chunks(seq, n):
    # http://stackoverflow.com/a/312464/190597 (Ned Batchelder)
    """ Yield successive n-sized chunks from seq."""
    for i in range(0, len(seq), n):
        yield seq[i:i + n]

def _fetch_from_nuccore(accessionList, destinationDir, emailAddress,
                        forceOverwrite=False, rettype="fasta", retmode="text",
                        fileExt=None, combinedFilePrefix=None, removeSeparateFiles=False,
                        chunkSize=1, api_key=None):
    """
        This function downloads and saves files from NCBI nuccore.
    """
    db = "nuccore"
    Entrez.email = emailAddress
    if api_key is not None:
        Entrez.api_key = api_key
    Entrez.tool = "https://github.com/broadinstitute/viral-ngs"

    maxChunkSize = 500

    # Conform to NCBI retreival guidelines by chunking into 500-accession chunks if
    # >500 accessions are specified and chunkSize is set to 1
    # Also clamp chunk size to 500 if the user specified a larger value.
    if chunkSize > maxChunkSize or (len(accessionList) > maxChunkSize and chunkSize == 1):
        chunkSize = maxChunkSize

    outEx = {"fasta": "fasta", "ft": "tbl", "gb": "gbk"}

    assert rettype in outEx.keys(
    ), "The return type requested, %s, is not compatible with the nuccore fetch." % rettype

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

    log.info("Fetching %s entries from GenBank: %s\n", str(len(accessionList)), ", ".join(accessionList[:10]))
    outputFiles = []

    for chunkNum, chunk in enumerate(_seq_chunks(accessionList, chunkSize)):

        accString = ",".join(chunk)

        # if the filename would be longer than Linux allows, simply say "chunk-chunkNum"
        if len(accString) + len(outputExtension) <= 254:
            outputFilePath = os.path.join(outputDirectory, accString + outputExtension)
        else:
            outputFilePath = os.path.join(outputDirectory, "chunk-{}".format(chunkNum) + outputExtension)

        if not forceOverwrite:
            log.info("not overwriting, checking for existence")
            assert not os.path.exists(outputFilePath), """File %s already exists. Consider removing
                this file or specifying a different output directory. The files for the accessions specified
                can be overwritten if you add --forceOverwrite flag. Processing aborted.""" % outputFilePath

        tryCount = 1
        while True:
            try:
                log.info("Fetching file %s: %s, try #%s", chunkNum + 1, accString, tryCount)
                handle = Entrez.efetch(db=db, rettype=rettype, id=accString)

                with open(outputFilePath, "w") as outf:
                    for line in handle:
                        outf.write(line)
                outputFiles.append(outputFilePath)
            except IOError:

                log.warning(
                    "Error fetching file %s: %s, try #%s probably because NCBI is too busy.", chunkNum + 1, accString,
                    tryCount)

                tryCount += 1
                if tryCount > 4:
                    log.warning("Tried too many times. Aborting.")
                    raise

                # if the fetch failed, wait a few seconds and try again.
                log.info("Waiting and retrying...")
                time.sleep(2)

                continue
            break

    # assert that we are not trying to remove the intermediate files without writing a combined file
    if removeSeparateFiles:
        assert combinedFilePrefix, """The intermediate files
            can only be removed if a combined file is written via --combinedFilePrefix"""

    # build a path to the combined genome file
    if combinedFilePrefix:
        concatenatedGenomeFilepath = os.path.join(outputDirectory, combinedFilePrefix + outputExtension)

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


def fetch_fastas_from_genbank(
        accessionList,
        destinationDir,
        emailAddress,
        forceOverwrite,
        combinedFilePrefix,
        removeSeparateFiles,
        fileExt=None,
        rettype="fasta",
        retmode="text",
        chunkSize=1, 
        api_key=None):
    return _fetch_from_nuccore(
        accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix,
        removeSeparateFiles, chunkSize, api_key=api_key)


def fetch_feature_tables_from_genbank(
        accessionList,
        destinationDir,
        emailAddress,
        forceOverwrite,
        combinedFilePrefix,
        removeSeparateFiles,
        fileExt=None,
        rettype="ft",
        retmode="text",
        chunkSize=1,
        api_key=None):
    return _fetch_from_nuccore(
        accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix,
        removeSeparateFiles, chunkSize, api_key=api_key)


def fetch_full_records_from_genbank(
        accessionList,
        destinationDir,
        emailAddress,
        forceOverwrite,
        combinedFilePrefix,
        removeSeparateFiles,
        fileExt=None,
        rettype="gb",
        retmode="text",
        chunkSize=1,
        api_key=None):
    return _fetch_from_nuccore(
        accessionList, destinationDir, emailAddress, forceOverwrite, rettype, retmode, fileExt, combinedFilePrefix,
        removeSeparateFiles, chunkSize, api_key=api_key)

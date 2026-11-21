import gzip
import os
from os.path import join
import shutil

import pysam

from viral_ngs.classify import krona as krona_module
from viral_ngs.core import file

CONTIG_REF_IDS = {
    "NODE_1": 0,
    "NODE_2": 1,
    "NODE_3": 2,
}

CONTIG_BAM_HEADER = {
    "HD": {"VN": "1.6"},
    "SQ": [
        {"SN": "NODE_1", "LN": 1000},
        {"SN": "NODE_2", "LN": 1000},
        {"SN": "NODE_3", "LN": 1000},
    ],
}

def write_contig_bam(path, records):
    with pysam.AlignmentFile(str(path), "wb", header=CONTIG_BAM_HEADER) as outf:
        for record in records:
            outf.write(record)


def contig_segment(
    query_name,
    reference_name="NODE_1",
    reference_start=0,
    cigar=None,
    mapq=60,
    nm=0,
    flag=0,
    query_length=None,
):
    cigar = cigar or [(pysam.CMATCH, 100)]
    if query_length is None:
        query_length = sum(
            length
            for op, length in cigar
            if op in (
                pysam.CMATCH,
                pysam.CINS,
                pysam.CSOFT_CLIP,
                pysam.CEQUAL,
                pysam.CDIFF,
            )
        )
    if query_length == 0:
        query_length = 100

    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.query_sequence = "A" * query_length
    record.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    record.flag = flag

    if flag & 0x4:
        record.reference_id = -1
        record.reference_start = -1
        record.mapping_quality = 0
    else:
        record.reference_id = CONTIG_REF_IDS[reference_name]
        record.reference_start = reference_start
        record.mapping_quality = mapq
        record.cigartuples = cigar
        if nm is not None:
            record.set_tag("NM", nm)

    return record


# @pytest.fixture(scope='module')
def krona():
    krona_tool = krona_module.Krona()
    krona_tool.install()
    return krona_tool


# @pytest.fixture(scope='module', params=['TestMetagenomicsSimple', 'TestMetagenomicsViralMix'])
def db_type(request):
    return request.param


# @pytest.fixture(scope='module')
def krona_db(request, tmpdir_module, krona, db_type):
    data_dir = join(file.get_test_input_path(), db_type)
    db_dir = join(data_dir, 'db')

    db = join(tmpdir_module, 'krona_db_{}'.format(db_type))
    os.mkdir(db)
    for d in ('names.dmp', 'nodes.dmp'):
        src = join(db_dir, 'taxonomy', d)
        dest = join(db, d)
        if os.path.isfile(src):
            os.symlink(src, dest)
    krona.create_db(db)
    return db


def taxonomy_db(request, tmpdir_module, db_type):
    taxonomy = join(tmpdir_module, db_type, 'taxonomy')
    shutil.copytree(join(file.get_test_input_path(), db_type, 'db', 'taxonomy'),
                    taxonomy)
    prot = join(taxonomy, 'accession2taxid', 'prot.accession2taxid')
    prot_gz = prot + '.gz'

    with open(prot, 'rb') as f_in:
        with gzip.open(prot_gz, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return taxonomy

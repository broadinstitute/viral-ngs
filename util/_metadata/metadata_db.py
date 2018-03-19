"""Management of the metadata database"""

import json
import collections
import builtins
import os
import functools
import operator
import io
import gzip
import contextlib
import warnings
import shlex

import util.misc
import util.file
from . import md_utils

import fs
import fs.multifs
import fs_s3fs

def metadata_dir():
    """Returns the string describing directories to which metadata was recorded, 
    as specified by the environment variable VIRAL_NGS_METADATA_PATH.
    Raises an error if the environment variable is not defined.
    """
    return os.environ['VIRAL_NGS_METADATA_PATH']

def metadata_dirs():
    """Return list of dirs to which metadata was recorded"""

    def read_from_file(x):
        if x.startswith('@'): return util.file.slurp_file(x[1:]).strip()
        return x

    return list(map(read_from_file, shlex.split(metadata_dir().strip())))

def _mask_secret_info(fs_url):
    """Mask any secret info, such as AWS keys, from fs_url. This is to keep such info from any printed logs."""
    if fs_url.startswith('s3://') and '@' in fs_url:
        fs_url = 's3://' + fs_url[fs_url.index('@')+1:]
    return fs_url

def metadata_dir_sanitized():
    """Return version `metadata_dir()` suitable for display in log messsages.  Any sensitive information like AWS keys will be scrubbed."""
    return ' '.join([_mask_secret_info(p) for p in metadata_dirs()])

def is_metadata_tracking_enabled():
    return bool(metadata_dirs())

def is_valid_step_record(d):
    """Test whether `d` is a dictionary containing all the expected elements of a step, as recorded by the code above"""
    return md_utils.dict_has_keys(d, 'format step') and \
        md_utils.dict_has_keys(d['step'], 'args step_id cmd_module')

def load_all_records():
    """Load all step records from the database."""

    records = []

    with contextlib.closing(fs.multifs.MultiFS()) as metadata_fs:
        for i, fsys in enumerate(metadata_dirs()):
            metadata_fs.add_fs('fs_{}'.format(i), fs.open_fs(fsys))

        fnames = sorted(set(metadata_fs.listdir(u'/')))

        for i, f in enumerate(fnames):
            if f.endswith('.json.gz'):
                step_record = util.file.from_json_gz(metadata_fs.getbytes(f))
                if is_valid_step_record(step_record):
                    records.append(step_record)
                else:
                    warnings.warn('Invalid step record ignored from {}'.format(f))
    return records

def store_step_record(step_data, write_obj=None):
    """Store step record to metadata database(s)"""
    json_fname = u'{}.json.gz'.format(step_data['step']['step_id'])
    json_data_gzipped = util.file.to_json_gz(step_data, write_obj=write_obj, filename=json_fname)
    
    for mdir in metadata_dirs():
        with fs.open_fs(mdir) as metadata_fs:
            metadata_fs.setbytes(json_fname, json_data_gzipped)

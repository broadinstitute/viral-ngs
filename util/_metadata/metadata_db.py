"""Management of the metadata database"""

import json
import collections
import builtins
import os

from . import md_utils

import fs
import fs_s3fs

def metadata_dir():
    """Returns the directory to which metadata was recorded, as specified by the environment variable VIRAL_NGS_METADATA_PATH.
    Raises an error if the environment variable is not defined.
    """
    return os.environ['VIRAL_NGS_METADATA_PATH']

def metadata_dir_sanitized():
    """Return version `metadata_dir()` suitable for display in log messsages.  Any sensitive information like AWS keys will be scrubbed."""
    return md_utils._mask_secret_info(metadata_dir())

def is_metadata_tracking_enabled():
    return 'VIRAL_NGS_METADATA_PATH' in os.environ
    # todo: check also that the only VIRAL_NGS_METADATA env vars are known ones

def is_valid_step_record(d):
    """Test whether `d` is a dictionary containing all the expected elements of a step, as recorded by the code above"""
    return md_utils.dict_has_keys(d, 'format step') and \
        md_utils.dict_has_keys(d['step'], 'args step_id cmd_module')

def load_all_records():
    """Load all step records from the database"""

    def byteify(input):
        if not hasattr(builtins, 'unicode'): return input
        if isinstance(input, dict):
            return {byteify(key): byteify(value)
                    for key, value in input.items()}
        elif isinstance(input, list):
            return [byteify(element) for element in input]
        elif isinstance(input, unicode):
            return input.encode('utf-8')
        else:
            return input

    records = []
    with fs.open_fs(metadata_dir()) as metadata_fs:
        for f in metadata_fs.listdir(u'/'):
            if f.endswith('.json'):
                json_str = metadata_fs.gettext(f)
                step_record = byteify(json.loads(json_str))
                if is_valid_step_record(step_record):
                    records.append(step_record)
    return records

def store_step_record(step_data, write_obj):
    json_fname = u'{}.json'.format(step_data['step']['step_id'])
    json_str = json.dumps(step_data, sort_keys=True, indent=4, default=write_obj)
    if hasattr(builtins, 'unicode'): json_str = unicode(json_str)
    with fs.open_fs(metadata_dir()) as metadata_fs:
        metadata_fs.settext(json_fname, json_str)
        
    
    

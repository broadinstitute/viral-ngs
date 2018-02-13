"""Management of the metadata database"""

from .md_utils import dict_has_keys

def is_valid_step_record(d):
    """Test whether `d` is a dictionary containing all the expected elements of a step, as recorded by the code above"""
    return dict_has_keys(d, 'format step') and \
        dict_has_keys(d['step'], 'args step_id cmd_module')

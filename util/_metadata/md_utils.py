"""Misc utils used by the metadata package"""

import collections
import util.misc

def _make_list(*x): return x

def _shell_cmd(cmd, *args, **kwargs):
    """Run a command and return its output; if command fails and `check` keyword arg is False, return the empty string."""
    out = ''
    result = util.misc.run_and_print(cmd.strip().split(), silent=True, *args, **kwargs)
    if result.returncode == 0:
        out = result.stdout
        if not isinstance(out, str):
            out = out.decode('utf-8')
        out = out.strip()
    return out

def _mask_secret_info(fs_url):
    """Mask any secret info, such as AWS keys, from fs_url. This is to keep such info from any printed logs."""
    if fs_url.startswith('s3://') and '@' in fs_url:
        fs_url = 's3://' + fs_url[fs_url.index('@')+1:]
    return fs_url

def dict_has_keys(d, keys_str):
    """Test whether a `d` is a dict containing all the given keys (given as tokens of `keys_str`)"""
    return isinstance(d, collections.Mapping) and d.keys() >= set(keys_str.split())


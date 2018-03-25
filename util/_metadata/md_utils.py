"""Misc utils used by the metadata package"""

import builtins
import os
import contextlib
import collections
import warnings
import traceback

import util.misc

@contextlib.contextmanager
 def errors_as_warnings(msg='handling metadata', mask_errors=True):
     """Context manager for wrapping non-essential functionality; if mask_errors=True (default), turns errors into
     warnings, so that errors in non-essential code do not stop primary functionality from being carried out.
     If runnning under pytest, exceptions are left as errors.
     
     Context manager returns a new empty list, to which exception string is appended if exception occurs; this lets later
     code see if an exception happened.
     """

    exc = []
    try:
        yield exc
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        # metadata recording is not an essential operation, so if anything goes wrong we just print a warning
        e_str = traceback.format_exc()
        exc.append(e_str)
        if mask_errors: 
            warnings.warn('Error in {} ({})'.format(msg, e_str))
        else:
            raise

def _shell_cmd(cmd, *args, **kwargs):
    """Run a command and return its output; if command fails and `check` keyword arg is False, return the empty string."""
    out = ''
    with errors_as_warnings(mask_errors=not kwargs.get('check')):
        result = util.misc.run_and_print(cmd.strip().split(), *args, **kwargs)
        if result.returncode == 0:
            out = result.stdout
            if not isinstance(out, str):
                out = out.decode('utf-8')
            out = out.strip()
        else:
            warnings.warn('Error running command: {}'.format(cmd))
    return out

def dict_has_keys(d, keys_str):
    """Test whether a `d` is a dict containing all the given keys (given as tokens of `keys_str`)"""
    return isinstance(d, collections.Mapping) and set(d.keys()) >= set(keys_str.split())

def tuple_key_matches(k, patterns):
    """Test whether a nested-dict key `k` (tuple of keys) matches one of the `patterns`.
    Each pattern is a tuple of strings giving possible values for each position in the key."""
    return any(all(str(k_elt) in p_elt.split() for k_elt, p_elt in zip(k, p) if p_elt) for p in patterns if len(k)>=len(p))

    

import os
import os.path
import collections
import stat
import shutil
import pwd
import traceback
import warnings
import tempfile
import multiprocessing

import util.file
from util.argparse_arg_types import OptionalFile
import util._metadata.hashing as hashing
from . import _log


# ** class FileArg
class FileArg(object):

    '''The value of an argparse parser argument denoting input or output file(s).  In addition to the string representing the
    argument value, keeps track of any filename(s) derived from the argument, and has methods for capturing metadata about the
    files they denote.'''
    
    def __init__(self, val, mode, compute_fnames):
        """Construct a FileArg.

        Args:
           val: the value of the command-line argument.  Most commonly this is just the filename of an input or output file of a command,
              but can also be e.g. the prefix for a group of files.
           mode: 'r' if `val` denotes input file(s), 'w' if to output files
           compute_fnames: function that will compute, from `val`, the list of actual filenames of the file(s) denoted by this 
             command-line argument.  See util.argparse_arg_types.
        """
        self.val, self.mode, self.compute_fnames = val, mode, compute_fnames
        self._set_up_hash_computation_for_pipes()


    def _set_up_hash_computation_for_pipes(self):
        """Set things up to let us compute file hash and size when the file is a pipe.

        If an input or output file is a pipe, this creates a problem for determining its hash and size.  A pipe can only be
        read or written once, but we need to run the data through both the command implementation and our hash-computing code.
        So we implement tee-like functionality to let us compute file hash/size as the command reads or writes the data.
        Even for a regular file this may be a more efficient way to compute the hash, as the file then does not need to be
        read twice; but that's only possible if the original command can handle piped input/output.

        Currently, this is only supported when the arg denotes a single file, the pipe.  Can be extended to other cases if such 
        use cases arise.
        """
        if self.fnames == [self.val] and util.file.ispipe(self.val):
            self.pipe_hasher = hashing.PipeHasher(old_pipe=self.val, mode=self.mode)
            self.val = self.pipe_hasher.new_pipe
            assert self.fnames == [self.val]
            
    @property
    def required_fnames(self):
        """List of required filename(s) denoted by this command-line argument."""
        return [str(f) for f in self.compute_fnames(self.val) if not isinstance(f, OptionalFile)]

    @property
    def optional_fnames(self):
        """List of optional filename(s) denoted by this command-line argument."""
        return [str(f) for f in self.compute_fnames(self.val) if isinstance(f, OptionalFile)]

    @property
    def fnames(self):
        """List all required and optional filename(s) denoted by this command-line argument"""
        return list(map(str, self.compute_fnames(self.val)))
        #return self.required_fnames + self.optional_fnames

    def gather_file_info(self, hasher, out_files_exist):
        """Return a dict representing metadata about the file(s) denoted by this argument.

        Args:
            hasher: callable for computing the hash value of a file
            out_files_exist: if False, don't expect output files to exist (because the command raised an exception)
        """

        def file2info(fname):
            """Compute a dictionary of info about one file"""
            file_info = dict(fname=fname, realpath=os.path.realpath(fname), abspath=os.path.abspath(fname))
            if self.mode=='r' or out_files_exist:
                try:
                    hash_val, size = None, None
                    if hasattr(self, 'pipe_hasher'):
                        hash_val, size = self.pipe_hasher.get_results()
                    else:
                        hash_val, size = hasher(fname)

                    file_info.update(hash=hash_val)

                    file_stat = os.stat(fname)
                    file_info.update(size=size,
                                     mtime=file_stat[stat.ST_MTIME], ctime=file_stat[stat.ST_CTIME])
                    file_info.update(owner=pwd.getpwuid(file_stat[stat.ST_UID]).pw_name)
                    file_info.update(inode=file_stat[stat.ST_INO], device=file_stat[stat.ST_DEV])
                except Exception:
                    if fname not in self.optional_fnames:
                        warnings.warn('Error getting file info for {} ({})'.format(fname, traceback.format_exc()))
                finally:
                    if hasattr(self, 'pipe_hasher'):
                        self.pipe_hasher.close()
            return file_info
        # end: def file2info(fname):

        return dict(__FileArg__=True, val=self.val, mode=self.mode, files=list(map(file2info, self.fnames)))
    # end: def gather_file_info(self, hasher, out_files_exist):

    @staticmethod
    def is_from_dict(val):
        """Tests whether `val` is a valid dict representation of a FileArg object (as constructed by gather_file_info() method)."""
        return isinstance(val, collections.Mapping) and '__FileArg__' in val and isinstance(val['files'], list)

    def __str__(self):
        return '{}({})'.format('InFile' if self.mode=='r' else 'OutFile', self.val)

    def __repr__(self): return str(self)

# end: class FileArg(object)
    

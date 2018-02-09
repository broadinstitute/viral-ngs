import os
import os.path
import collections
import stat
import pwd

from .md_utils import _make_list


# ** class FileArg
class FileArg(object):

    '''The value of an argparse parser argument denoting input or output file(s).  In addition to the string representing the
    argument value, keeps track of any filename(s) derived from the argument, and has methods for capturing metadata about the
    files they denote.'''
    
    def __init__(self, val, mode, compute_fnames=_make_list):
        """Construct a FileArg.

        Args:
           val: the value of the command-line argument.  Most commonly this is just the filename of an input or output file of a command,
              but can also be e.g. the prefix for a group of files.
           mode: 'r' if `val` denotes input file(s), 'w' if to output files
           compute_fnames: function that will compute, from `val`, the list of actual filenames of the file(s) denoted by this 
             command-line argument.  By default, this is just one file and `val` contains its full name.  But `val` can be a 
             common prefix for a set of files with a given list of suffixes, or `val` can be a directory denoting all the files
             in the directory or just those matching a wildcard; and in those cases, compute_fnames will compute the actual file names
             by some non-trivial operation.
        """
        self.val, self.mode, self.compute_fnames = val, mode, compute_fnames

    @property
    def fnames(self):
        """List of filename(s) denoted by this command-line argument."""
        return self.compute_fnames(self.val)

    def gather_file_info(self, hasher, out_files_exist):
        """Return a dict representing metadata about the file(s) denoted by this argument.

        Args:
            hasher: callable for computing the hash value of a file
            out_files_exist: if False, don't expect output files to exist (because the command raised an exception)
        """

        def file2dict(file_arg):
            """Compute a dictionary of info about one file"""
            file_info = dict(fname=file_arg, realpath=os.path.realpath(file_arg), abspath=os.path.abspath(file_arg))
            if self.mode=='r' or out_files_exist:
                file_info.update(hash=hasher(file_arg))

                try:
                    file_stat = os.stat(file_arg)
                    file_info.update(size=file_stat[stat.ST_SIZE],
                                     mtime=file_stat[stat.ST_MTIME], ctime=file_stat[stat.ST_CTIME])
                    file_info.update(owner=pwd.getpwuid(file_stat[stat.ST_UID]).pw_name)
                    file_info.update(inode=file_stat[stat.ST_INO], device=file_stat[stat.ST_DEV])
                except Exception:
                    _log.warning('Error getting file info for {} ({})'.format(file_arg, traceback.format_exc()))
            return file_info
        # end: def file2dict(file_arg):

        return dict(__FileArg__=True, val=self.val, mode=self.mode, files=list(map(file2dict, self.fnames)))
    # end: def gather_file_info(self, hasher, out_files_exist):

    @staticmethod
    def is_from_dict(val):
        """Tests whether `val` is a valid dict representation of a FileArg object (as constructed by gather_file_info() method)."""
        return isinstance(val, collections.Mapping) and '__FileArg__' in val and isinstance(val['files'], list)

    def __str__(self):
        return '{}({})'.format('InFile' if self.mode=='r' else 'OutFile', self.val)

    def __repr__(self): return str(self)

# end: class FileArg(object)
    

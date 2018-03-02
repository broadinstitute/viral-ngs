"""Hashing of files"""

import os
import os.path
import io
import tempfile
import stat
import traceback
import warnings
import shutil
import hashlib
import errno
import multiprocessing
import time

import util.file

from util._metadata import _log

class Hasher(object):
    """Manages computation of file hashes.

    It might also cache the actual contents of files, depending on size.
    """

    def __init__(self, hash_algorithm='sha1'):
        self.hash_algorithm = hash_algorithm

    def __call__(self, file):
        file_hash = ''
        try:
            if os.path.isfile(file) and not stat.S_ISFIFO(os.stat(file).st_mode):
                file_hash = self.hash_algorithm + '_' + util.file.hash_file(file, hash_algorithm=self.hash_algorithm)
        except Exception:
            warnings.warn('Cannot compute hash for {}: {}'.format(file, traceback.format_exc()))
        return file_hash

# end: class Hasher(object)

def compute_hash_and_size(fname, copy_data_to=None, copy_pipe_may_break=False,
                          write_result_to=None, done_file=None, buf_size=io.DEFAULT_BUFFER_SIZE, hash_algorithm = 'sha1'):
    '''Compute the hashsum and of a file, or of data in a pipe.

    Args:
        fname: path to file or named pipe
        copy_data_to: if not None, data is copied to this pipe
        copy_pipe_may_break: if True, if pipe errors happen writing to `copy_data_to`, we just stop writing there (but still
           finish reading the input and computing its hash and size)
        write_result_to: write results to this file
        done_file: touch this file after writing results
        buf_size: process data in chunks of this size
        hash_algorithm: hash algorithm to use
    '''
    
    assert not hasattr(hashlib, 'algorithms_available') or hash_algorithm in hashlib.algorithms_available
    hasher = getattr(hashlib, hash_algorithm)()
    data_size = 0
    _log.debug('HASHING: orig={} new={}'.format(fname, copy_data_to))
    with open(fname, 'rb') as f, open(copy_data_to or os.devnull, 'wb') as copy_data_out:
        _log.debug('OPENED OK')
        while True:
            data = f.read(buf_size)
            _log.debug('READ OK')
            data_size += len(data)
            _log.debug('GOT DATA: {} {} {}'.format(len(data), data_size, data))
            if data:
                hasher.update(data)
                if copy_data_to: 
                    try:
                        copy_data_out.write(data)
                    except IOError as e:
                        if e.errno == errno.EPIPE:
                            copy_data_to = None
                        else:
                            raise
                _log.debug('PROCESSED DATA: {} {}'.format(len(data), data))
            else:
                result = (hasher.hexdigest(), data_size)
                if write_result_to:
                    _log.debug('SENDING DATA BACK')
                    util.file.dump_file(write_result_to, '\n'.join(map(str, result)))
                    if done_file: util.file.make_empty(done_file)
                return result
            
class PipeHasher(object):
    """Manages the hashing of data in a pipe"""

    def __init__(self, old_pipe, mode):
        self.old_pipe = old_pipe
        self.pipe_dir = tempfile.mkdtemp()
        self.new_pipe = os.path.join(self.pipe_dir, '0.pipe')
        os.mkfifo(self.new_pipe)
        hasher_reads, hasher_writes = (old_pipe, self.new_pipe) if mode=='r' else (self.new_pipe, old_pipe)
        self.pipe_result_fname = os.path.join(self.pipe_dir, 'hash_result')
        self.pipe_done_fname = os.path.join(self.pipe_dir, 'hash_result_done')
        self.pipe_hasher_proc = multiprocessing.Process(target=compute_hash_and_size, 
                                                       kwargs=dict(fname=hasher_reads, copy_data_to=hasher_writes,
                                                                   write_result_to=self.pipe_result_fname,
                                                                   done_file=self.pipe_done_fname))
        self.pipe_hasher_proc.start()

    def get_results(self):
        """Return the hash and size of the data (or '' and -1 respectivevly in case of failure), and clean up."""

        hash_val, size = '', -1
        self.pipe_hasher_proc.join(5)
        if os.path.isfile(self.pipe_done_fname):
            hash_val, size = util.file.slurp_file(self.pipe_result_fname).strip().split()
            size = int(size)

        return hash_val, size

    def close(self):
        """Remove the temp dir"""
        if self.pipe_hasher_proc.is_alive():
            self.pipe_hasher_proc.terminate()
        shutil.rmtree(self.pipe_dir, ignore_errors=True)

"""Hashing of files"""

import os
import os.path
import sys
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
import tools.samtools

from util._metadata import _log

class Hasher(object):
    """Manages computation of file hashes.
    """

    def __init__(self, hash_algorithm='sha1'):
        self.hash_algorithm = hash_algorithm

    def __call__(self, fname):
        file_hash = ''
        file_size = 0
        try:
            if os.path.isfile(fname) and not stat.S_ISFIFO(os.stat(fname).st_mode):
                if False and fname.endswith('.bam'):
                    file_hash, file_size = canonicalize_bam(fname, self.hash_algorithm)
                else:
                    file_hash = self.hash_algorithm + '_' + util.file.hash_file(fname, hash_algorithm=self.hash_algorithm)
                    file_size = os.path.getsize(fname)
        except Exception:
            warnings.warn('Cannot compute hash for {}: {}'.format(fname, traceback.format_exc()))
        return file_hash, file_size

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
    with open(fname, 'rb') as f, open(copy_data_to or os.devnull, 'wb') as copy_data_out:
        while True:
            data = f.read(buf_size)
            data_size += len(data)
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
            else:
                result = (hasher.hexdigest(), data_size)
                if write_result_to:
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


def canonicalize_bam(fname, hash_algorithm):  # pragma: no cover
    """Computed hash and size for a canonicalized bam, where details such as exact command lines used to produce the file
    are ignored."""
    one-time temp files) 
    with util.file.tempfnames(suffixes=('.oldhr.txt','.newhdr.txt', 'canon.bam')) as (old_header_fname, 
                                                                                      new_header_fname, canon_bam):
        # make this a plugin call
        tools.samtools.SamtoolsTool().dumpHeader(fname, old_header_fname)
        with open(new_header_fname, 'wb') as out_h, open(old_header_fname, 'rb') as in_h:
            for line in in_h:
                line = line.decode("latin-1")
                if line.startswith('@PG'):
                    line = '\t'.join(filter(lambda s: not s.startswith('CL:'), line.rstrip('\n').split('\t')))+'\n'
                if line.startswith('@SQ'):
                    line = '\t'.join(filter(lambda s: not s.startswith('AS:'), line.rstrip('\n').split('\t')))+'\n'
                out_h.write(line.encode("latin-1"))

        tools.samtools.SamtoolsTool().reheader_no_PG(fname, new_header_fname, canon_bam)
        file_hash = hash_algorithm + '_bamcanon_' + util.file.hash_file(canon_bam, hash_algorithm=hash_algorithm)
        file_size = os.path.getsize(canon_bam)
        return file_hash, file_size

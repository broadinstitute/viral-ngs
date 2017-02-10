'''A few miscellaneous tools. '''
from __future__ import print_function, division  # Division of integers with / should never round!
import collections
import itertools
import logging
import os
import re
import subprocess
import multiprocessing
import sys

import util.file

log = logging.getLogger(__name__)

__author__ = "dpark@broadinstitute.org"


def unique(items):
    ''' Return unique items in the same order as seen in the input. '''
    seen = set()
    for i in items:
        if i not in seen:
            seen.add(i)
            yield i


def histogram(items):
    ''' I count the number of times I see stuff and return a dict of counts. '''
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
    return out


def freqs(items, zero_checks=None):
    ''' Given a list of comparable, non-unique items, produce an iterator of
            (item, count, freq) tuples.
            item is a unique instance of one of the items seen on input
            count is a positive integer describing the number of times this item was observed
            freq is count / sum(count) for all outputs.
        If zero_checks is specified, then the output iterator will emit tuples for the
        items in zero_checks even if they do not appear in the input.  If they are not in
        the input, they will be emitted with a zero count and freq.
        See histogram(items)
    '''
    zero_checks = zero_checks or set()

    tot = 0
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
        tot += 1
    for k, v in out.items():
        yield (k, v, float(v) / tot)
    for i in zero_checks:
        if i not in out:
            yield (i, 0, 0.0)


def intervals(i, n, l):
    ''' Divide something of length l into n equally sized parts and return the
        start-stop values of the i'th part.  Values are 1-based.  Each part
        will be adjacent and non-overlapping with the next part. i must be a
        number from 1 to n.
    '''
    assert 1 <= i <= n and l >= n
    part_size = l // n
    start = 1 + part_size * (i - 1)
    stop = part_size * i
    if i == n:
        stop = l
    return (start, stop)

# from http://stackoverflow.com/a/312467


def pairwise(iterable):
    """ from itertools recipes
        s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    if hasattr(itertools, 'izip'):
        # Python 2
        return itertools.izip(a, b)
    else:
        # Python 3
        return zip(a, b)


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    it = iter(iterator)
    item = list(itertools.islice(it, batch_size))
    while item:
        yield item
        item = list(itertools.islice(it, batch_size))


def list_contains(sublist, list_):
    """Tests whether sublist is contained in full list_."""

    for i in range(len(list_)-len(sublist)+1):
        if sublist == list_[i:i+len(sublist)]:
            return True
    return False


try:
    from subprocess import run
except ImportError:
    CompletedProcess = collections.namedtuple(
        'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])

    def run(args, stdin=None, stdout=None, stderr=None, shell=False,
            env=None, cwd=None, timeout=None, check=False):
        '''A poor man's substitute of python 3.5's subprocess.run().

        Should only be used for capturing stdout. If stdout is unneeded, a
        simple subprocess.call should suffice.
        '''
        assert stdout is not None, (
            'Why are you using this util function if not capturing stdout?')

        stdout_pipe = stdout == subprocess.PIPE
        stderr_pipe = stderr == subprocess.PIPE
        # A little optimization when we don't need temporary files.
        if stdout_pipe and (
                stderr == subprocess.STDOUT or stderr is None):
            try:
                output = subprocess.check_output(
                    args, stdin=stdin, stderr=stderr, shell=shell,
                    env=env, cwd=cwd)
                return CompletedProcess(args, 0, output, b'')
            # Py3.4 doesn't have stderr attribute
            except subprocess.CalledProcessError as e:
                if check:
                    raise
                returncode = e.returncode
                stderr_text = getattr(e, 'stderr', b'')
                return CompletedProcess(args, e.returncode, e.output, stderr_text)

        # Otherwise use temporary files as buffers, since subprocess.call
        # cannot use PIPE.
        if stdout_pipe:
            stdout_fn = util.file.mkstempfname('.stdout')
            stdout = open(stdout_fn, 'wb')
        if stderr_pipe:
            stderr_fn = util.file.mkstempfname('.stderr')
            stderr = open(stderr_fn, 'wb')
        try:
            returncode = subprocess.call(
                args, stdin=stdin, stdout=stdout,
                stderr=stderr, shell=shell, env=env, cwd=cwd)
            if stdout_pipe:
                stdout.close()
                with open(stdout_fn, 'rb') as f:
                    output = f.read()
            else:
                output = ''
            if stderr_pipe:
                stderr.close()
                with open(stderr_fn, 'rb') as f:
                    error = f.read()
            else:
                error = ''
            if check and returncode != 0:
                print(output.decode("utf-8"))
                print(error.decode("utf-8"))
                try:
                    raise subprocess.CalledProcessError(
                        returncode, args, output, error) #pylint: disable-msg=E1121
                except TypeError: # py2 CalledProcessError does not accept error
                    raise subprocess.CalledProcessError(
                        returncode, args, output)
            return CompletedProcess(args, returncode, output, error)
        finally:
            if stdout_pipe:
                stdout.close()
                os.remove(stdout_fn)
            if stderr_pipe:
                stderr.close()
                os.remove(stderr_fn)


def run_and_print(args, stdout=None, stderr=None,
                  stdin=None, shell=False, env=None, cwd=None,
                  timeout=None, silent=False, buffered=False, check=False,
                  loglevel=None):
    '''Capture stdout+stderr and print.

    This is useful for nose, which has difficulty capturing stdout of
    subprocess invocations.
    '''
    if loglevel:
        silent = False
    if not buffered:
        if check and not silent:
            try:
                result = run(
                    args,
                    stdin=stdin,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    env=env,
                    cwd=cwd,
                    timeout=timeout,
                    check=check
                )
                print(result.stdout.decode('utf-8'))
                try:
                    print(result.stderr.decode('utf-8'))
                except AttributeError:
                    pass
            except subprocess.CalledProcessError as e:
                if loglevel:
                    try:
                        log.log(loglevel, result.stdout.decode('utf-8'))
                    except NameError:
                        # in some situations, result does not get assigned anything
                        pass
                    except AttributeError:
                        log.log(loglevel, result.output.decode('utf-8'))
                else:
                    print(e.output.decode('utf-8'))
                    try:
                        print(e.stderr.decode('utf-8'))
                    except AttributeError:
                        pass
                    sys.stdout.flush()
                raise(e)
        else:
            result = run(args, stdin=stdin, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, env=env, cwd=cwd,
                         timeout=timeout, check=check)
            if not silent and not loglevel:
                print(result.stdout.decode('utf-8'))
                sys.stdout.flush()
            elif loglevel:
                log.log(loglevel, result.stdout.decode('utf-8'))

    else:
        CompletedProcess = collections.namedtuple(
        'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])

        process = subprocess.Popen(args, stdin=stdin, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, env=env,
                                    cwd=cwd)
        output = []
        while process.poll() is None:
            for c in iter(process.stdout.read, b""):
                output.append(c)
                if not silent:
                   print(c.decode('utf-8'), end="") # print for py3 / p2 from __future__
                sys.stdout.flush() # flush buffer to stdout

        # in case there are still chars in the pipe buffer
        for c in iter(lambda: process.stdout.read(1), b""):
           if not silent:
               print(c, end="")
           sys.stdout.flush() # flush buffer to stdout

        if check and process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, args,
                                                b''.join(output))
        result = CompletedProcess(
            args, process.returncode, b''.join(output), b'')

    return result


def run_and_save(args, stdout=None, stdin=None,
                 outf=None, stderr=None, preexec_fn=None,
                 close_fds=False, shell=False, cwd=None, env=None, check=True):
    assert outf is not None

    sp = subprocess.Popen(args,
                          stdin=stdin,
                          stdout=outf,
                          stderr=subprocess.PIPE,
                          preexec_fn=preexec_fn,
                          close_fds=close_fds,
                          shell=shell,
                          cwd=cwd,
                          env=env)
    out, err = sp.communicate()

    # redirect stderror to stdout so it can be captured by nose
    if err:
        sys.stdout.write(err.decode("UTF-8"))

    if sp.returncode != 0 and check:
        raise subprocess.CalledProcessError(sp.returncode, str(args[0]))

    return sp

class FeatureSorter(object):
    ''' This class helps sort genomic features. It's not terribly optimized
        for speed or anything. Slightly inspired by calhoun's MultiSequenceRangeMap.
    '''
    def __init__(self, collection=None):
        self.seqids = []
        self.seq_to_features = {}
        self.seq_to_breakpoints = {}
        self.dirty = False
        if collection is not None:
            for args in collection:
                self.add(*args)

    def add(self, c, start, stop, strand='+', other=None):
        ''' Add a "feature", which is a chrom,start,stop,strand tuple (with
            optional other info attached)
        '''
        if c not in self.seq_to_features:
            self.seqids.append(c)
            self.seq_to_features[c] = []
            self.seq_to_breakpoints[c] = set()
            #self.seq_to_breakpoints[c].add(1) # do we want the empty interval in front?
        self.seq_to_features[c].append((int(start), int(stop), strand, other))
        self.seq_to_breakpoints[c].add(start)
        self.seq_to_breakpoints[c].add(stop+1)
        self.dirty = True

    def _cleanup(self):
        if self.dirty:
            self.dirty = False
            for c in self.seqids:
                self.seq_to_features[c].sort()

    def get_seqids(self):
        return tuple(self.seqids)

    def get_features(self, c=None, left=0, right=float('inf')):
        ''' Get all features on all chromosomes in sorted order. Chromosomes
            are emitted in order of first appearance (via add). Features on
            each chromosome are emitted in start, then stop order.  If
            boundaries are specified, we restrict to features that contain
            the specified interval.
        '''
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for start, stop, strand, other in self.seq_to_features[c]:
                if stop>=left and start<=right:
                    yield (c, start, stop, strand, other)

    def get_intervals(self, c=None):
        ''' Get all intervals on the reference where the overlapping feature
            set remains the same. Output will be sorted, adjacent intervals
            and will describe how many and which features overlap it.
        '''
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for left, right in pairwise(sorted(self.seq_to_breakpoints[c])):
                right = right - 1
                features = list(self.get_features(c, left, right))
                yield (c, left, right, len(features), features)


def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """
    try:
        with open('/proc/self/status') as f:
            status = f.read()
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return min(res, multiprocessing.cpu_count())
    except IOError:
        pass

    return multiprocessing.cpu_count()

def which(application_binary_name):
    """
        Similar to the *nix "which" command, 
        this function finds the first executable binary present
        in the system PATH for the binary specified. 
        It differs in that it resolves symlinks.
    """
    path=os.getenv('PATH')
    for path in path.split(os.path.pathsep):
        full_path=os.path.join(path, application_binary_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            link_resolved_path = os.path.realpath(full_path)
            return link_resolved_path
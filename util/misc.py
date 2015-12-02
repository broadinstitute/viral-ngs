'''A few miscellaneous tools. '''
from __future__ import division  # Division of integers with / should never round!
import collections
import itertools
import subprocess

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


try:
    from subprocess import run
except ImportError:
    CompletedProcess = collections.namedtuple(
        'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])

    def run(args, stdin=None, stdout=None, stderr=None, shell=False,
            env=None, cwd=None, timeout=None):
        '''A poor man's substitute of python 3.5's subprocess.run().

        Definitely a poor man's substitute because stdout and stderr are
        forcibly merged into stdout and capturing always takes place even when
        they should require subprocess.PIPE assignments, but the interface is
        fairly similar.
        '''
        try:
            output = subprocess.check_output(
                args, stdin=stdin, stderr=subprocess.STDOUT, shell=shell,
                env=env, cwd=cwd)
            returncode = 0
        except subprocess.CalledProcessError as e:
            output = e.output
            returncode = e.returncode

        return CompletedProcess(args, returncode, output, '')


def run_and_print(args, stdin=None, shell=False, env=None, cwd=None,
                  timeout=None):
    '''Capture stdout+stderr and print.

    This is useful for nose, which has difficulty capturing stdout of
    subprocess invocations.
    '''
    result = run(args, stdin=stdin, stdout=subprocess.PIPE,
                 stderr=subprocess.STDOUT, env=env, cwd=cwd, timeout=timeout)
    print(result.stdout.decode('utf-8'))
    return result

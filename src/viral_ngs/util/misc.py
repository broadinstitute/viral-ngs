'''A few miscellaneous tools. '''

import itertools
# import scipy # Commented out until installation on Travis is working

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

def freqs(items, zero_checks = set()):
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
    tot = 0
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
        tot += 1
    for k,v in out.items():
        yield (k,v,float(v)/tot)
    for i in zero_checks:
        if i not in out:
            yield (i,0,0.0)

def intervals(i, n, l):
    ''' Divide something of length l into n equally sized parts and return the
        start-stop values of the i'th part.  Values are 1-based.  Each part
        will be adjacent and non-overlapping with the next part. i must be a
        number from 1 to n.
    '''
    assert 1 <= i <= n and l>=n
    part_size = l//n
    start = 1 + part_size * (i-1)
    stop = part_size * i
    if i==n:
        stop = l
    return (start,stop)


try:
    # Python 3.4
    from statistics import mean, median
except ImportError:
    # Python <3.4, avoid numpy if these two methods are all we really need
    def mean(l):
        if len(l)>0:
            return float(sum(l))/len(l)
        else:
            raise Exception("empty list for mean")
    def median(l):
        if len(l)>0:
            half = len(l) / 2
            l.sort()
            if len(l) % 2 == 0:
                return (l[half-1] + l[half]) / 2.0
            else:
                return l[half]
        else:
            raise Exception("empty list for median")

# from http://stackoverflow.com/a/312467
def batch_iterator(iterator, batch_size) :
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

def fisher_exact(contingencyTable2x2) :
    """ Return two-tailed p-value for a 2 x 2 contingency table using Fisher's
        exact test. Input is a sequence of 2 sequences of size 2.
    """
    # Commented out until scipy installation is working:
    # return scipy.stats.fisher_exact(contingencyTable)[1]

    return python_fisher_exact(contingencyTable2x2)

def chi2_contingency(contingencyTable) :
    """ Return two-tailed p-value for an n x m contingency table using chi-square
        distribution. Not recommended if any of the counts or expected counts are
        less than 5. Input is a sequence of n sequences of size m.
    """
    # Commented out scipy version until scipy installation is working:
    # return scipy.stats.chi2_contingency(contingencyTable)[1]
    
    return 1.0

from math import exp, log, pi

def log_stirling(n) :
    """Return Stirling's approximation for log(n!) use up to n^3 term."""
    return n * log(n) - n + 0.5 * log(2 * pi * n) +  \
           1 / 12 / n - 1 / 360 / n ** 3

def python_fisher_exact(contingencyTable2x2) :
    """ Return two-tailed p-value for a 2 x 2 contingency table using Fisher's
            exact test. Input is a sequence of 2 sequences of size 2.
    """
    a, b = contingencyTable2x2[0][0], contingencyTable2x2[0][1]
    c, d = contingencyTable2x2[1][0], contingencyTable2x2[1][1]

    def log_choose(n, k) :
        k = min(k, n - k)
        if k <= 10 :
            result = 0.0
            for ii in range(1, k + 1) :
                result += log(n - ii + 1) - log(ii)
        else :
            result = log_stirling(n) - log_stirling(k) - log_stirling(n-k)
        return result
    
    def mydhyper(a, s, t, n) :
        # Probability that intersection of subsets of size s and t of a set with
        # n elements has size exactly a.
        return exp(log_choose(s, a) + log_choose(n - s, t - a) - log_choose(n, t))

    s = a + b
    t = a + c
    n = a + b + c + d
    result = p0 = mydhyper(a, s, t, n)

    # Count up from 0 and down from min(s, t) adding all smaller p-vals
    for x in range(a) :
        prob = mydhyper(x, s, t, n)
        if prob <= p0 :
            result += prob
        else :
            break
    for y in range(min(s, t), a, -1) :
        prob = mydhyper(y, s, t, n)
        if prob <= p0 :
            result += prob
        else :
            break
    return result


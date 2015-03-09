'''A few pure-python statistical tools to avoid the need to install scipy. '''
from __future__ import division # Division of integers with / should never round!
from math import exp, log, pi
import itertools
# import scipy # Commented out until installation on Travis is working

__author__ = "dpark@broadinstitute.org, irwin@broadinstitute.org"

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
            half = len(l) // 2
            l.sort()
            if len(l) % 2 == 0:
                return (l[half-1] + l[half]) / 2.0
            else:
                return l[half]
        else:
            raise Exception("empty list for median")

def chi2_contingency(contingencyTable) :
    """ Return two-tailed p-value for an n x m contingency table using chi-square
        distribution. Not recommended if any of the counts or expected counts are
        less than 5. Input is a sequence of n sequences of size m.
    """
    # Commented out scipy version until scipy installation is working:
    # return scipy.stats.chi2_contingency(contingencyTable)[1]
    
    return 1.0

def log_stirling(n) :
    """Return Stirling's approximation for log(n!) using up to n^7 term.
       Provides exact right answer (when rounded to int) for 16! and lower.
       Correct to 10 digits for 5! and above."""
    n2 = n * n
    n3 = n * n2
    n5 = n3 * n2
    n7 = n5 * n2
    return n * log(n) - n + 0.5 * log(2 * pi * n) +  \
           1 / 12.0 / n - 1 / 360.0 / n3 + 1 / 1260.0 / n5 - 1 / 1680.0 / n7

def log_choose(n, k) :
    k = min(k, n - k)
    if k <= 10 :
        result = 0.0
        for ii in range(1, k + 1) :
            result += log(n - ii + 1) - log(ii)
    else :
        result = log_stirling(n) - log_stirling(k) - log_stirling(n-k)
    return result

def fisher_exact(contingencyTable2x2) :
    """ Return two-tailed p-value for a 2 x 2 contingency table using Fisher's
            exact test. Input is a sequence of 2 sequences of size 2.
    """
    # Python code below is accurate to around 1e-11 and is much faster than
    # the scipy version. However, if for some reason we want to use the scipy
    # version instead, here's the call:
    # return scipy.stats.fisher_exact(contingencyTable)[1]
    
    a, b = contingencyTable2x2[0][0], contingencyTable2x2[0][1]
    c, d = contingencyTable2x2[1][0], contingencyTable2x2[1][1]

    s = a + b
    t = a + c
    n = a + b + c + d
    logChooseNT = log_choose(n, t)

    def mydhyper(a, s, t, n) :
        # Probability that intersection of subsets of size s and t of a set with
        # n elements has size exactly a.
        if max(0, s + t - n) <= a <= min(s, t) :
            return exp(log_choose(s, a) + log_choose(n - s, t - a) - logChooseNT)
        else :
            return 0

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


'''A few pure-python statistical tools to avoid the need to install scipy. '''
from __future__ import division # Division of integers with / should never round!
from math import exp, log, pi, sqrt
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

def product(iter) :
    prod = 1
    for x in iter :
        prod *= x
    return prod

def chi2_contingency(contingencyTable, correction = True) :
    """ contingencyTable is a sequence of m sequences, each of length n.
        Return an estimate using the chi-square distribution of the two-tailed 
            p-value for an m x n contingency table against the null hypothesis
            that the row and column criteria are independent. This is not as 
            accurate as fisher_exact, but is faster (and is implemented for all 
            m and n).
        If correction is True and there is 1 degree of freedom, apply Yates's
            correction for continuity, i.e., adjust each observed value
            by 0.5 towards the corresponding expected value, which brings
            the result closer to the Fisher exact test result.
        Not recommended if any of the counts or expected counts are
            less than 5.
    """
    # scipy equivalent:  scipy.stats.chi2_contingency(contingencyTable)[1]
    
    if len(contingencyTable) == 0 :
        return 1.0
    if len(set(map(len, contingencyTable))) != 1 :
        raise ValueError('Not all rows have the same length')
    
    # Eliminate rows and columns with 0 sum
    colSums = [sum(row[col] for row in contingencyTable)
               for col in range(len(contingencyTable[0]))]
    table = [[x for x, colSum in zip(row, colSums) if colSum != 0]
             for row in contingencyTable if sum(row) != 0]

    if len(table) < 2 or len(table[0]) < 2 :
        return 1.0

    m = len(table)
    n = len(table[0])
    rowSums = [sum(row) for row in table]
    colSums = [sum(row[col] for row in table) for col in range(n)]
    N = sum(rowSums)
    expect = [[rowSums[i] * colSums[j] / N for j in range(n)] for i in range(m)]
    if correction and m == n == 2 :
        def corr(i, j) :
            if expect[i][j] > table[i][j] :
                return min(table[i][j] + 0.5, expect[i][j])
            else :
                return max(table[i][j] - 0.5, expect[i][j])
        table = [[corr(i, j) for j in range(n)] for i in range(m)]
    chisq = sum((table[i][j] - expect[i][j]) ** 2 / expect[i][j]
                for j in range(n)
                for i in range(m))
    pval = 1 - pchisq(chisq, (m - 1) * (n - 1))
    return pval

def fisher_exact(contingencyTable2xn) :
    """ contingencyTable2xn is a sequence of 2 sequences, each of length n.
        Return the two-tailed p-value for a 2 x n contingency table against the
            null hypothesis that the row and column criteria are independent,
            using Fisher's exact test.
        For n larger than 2, this is very slow, O(S^(n-1)), where S is the
            smaller of the two row sums. Better to use chi2_contingency unless
            one of the row sums is small.
    """
    if len(contingencyTable2xn) != 2 or \
       len(contingencyTable2xn[0]) != len(contingencyTable2xn[1])  :
        raise NotImplementedError('Input must be 2 rows of equal length.')

    n = len(contingencyTable2xn[0])
    
    if n == 2 :
        return fisher_exact_2x2(contingencyTable2xn)

    # There are many optimizations possible for the following code, but it would
    #     still be O(S^(n-1)) so it would still be too slow for anything
    #     sizeable, and it's usable as it for small things.

    # Put row with smaller sum first. Makes the loop iterations simpler.
    table = list(contingencyTable2xn)
    if sum(table[0]) > sum(table[1]) :
        table.reverse()

    rowSum0 = sum(table[0])
    rowsum1 = sum(table[1])
    colSums = [table[0][i] + table[1][i] for i in range(n)]

    logChooseNrowSum = log_choose(rowSum0 + rowsum1, rowSum0)

    def prob_of_table(firstRow) :
        return exp(sum(log_choose(cs, a) for cs, a in zip(colSums, firstRow)) -
                   logChooseNrowSum)

    p0 = prob_of_table(table[0])

    result = 0
    for firstRowM1 in itertools.product(*[range(rowSum0 + 1) for i in range(n - 1)]) :
        lastElmt = rowSum0 - sum(firstRowM1)
        if lastElmt < 0 :
            continue
        firstRow = list(firstRowM1) + [lastElmt]
        if any(x > colSum for x, colSum in zip(firstRow, colSums)) :
            continue
        prob = prob_of_table(firstRow)
        if prob <= p0 + 1e-9 : # Handle floating point round off in equality check
            result += prob
    return result

def fisher_exact_2x2(contingencyTable2x2) :
    """ Return two-tailed p-value for a 2 x 2 contingency table using Fisher's
            exact test. Input is a sequence of 2 sequences of size 2.
    """
    # Python code below is accurate to around 1e-11 and is much faster than
    # the scipy version. However, if for some reason we want to use the scipy
    # version instead, here's the call:
    # scipy equivalent: scipy.stats.fisher_exact(contingencyTable2x2)[1]
    
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

def factorial(n) :
    """Return n factorial exactly up to n = 16, otherwise approximate to 14 digits.
       n must be a non-negative integer."""
    if n < 0 or n != int(n) :
        raise ValueError('%s is not a non-negative integer' % n)
    return int(round(exp(log_stirling(n)))) if n > 0 else 1

def gamma(s) :
    """ Gamma function = integral from 0 to infinity of t ** (s-1) exp(-t) dt.
        Implemented only for s >= 0,
    """
    # Accurate to better than 1 in 1e11
    # scipy equivalent: scipy.special.gamma(s)
    if s <= 0 :
        raise ValueError('%s is not positive' % s)
    if s == int(s) :
        return factorial(int(s - 1))
    elif 2 * s == int(2 * s) :
        return sqrt(pi) * factorial(int(2 * s - 1)) / factorial(int(s - 0.5)) /\
                   4 ** (s - 0.5)
    else :
        # stirling is more accurate for larger s, so call it for larger s and
        # use gamma recursion to get back to lower value.
        return exp(log_stirling(s + 9)) / product(s + i for i in range(10))

def gammainc(s, x) :
    """ Lower incomplete gamma function = 
            integral from 0 to x of t ** (s-1) exp(-t) dt divided by gamma(s),
        i.e., the fraction of gamma that you get if you integrate only until
            x instead of all the way to infinity.
        Implemented only for s >= 0.
    """
    # scipy equivalent: scipy.special.gammainc(s,x)
    if s <= 0 :
        raise ValueError('%s is not positive' % s)
    if x < 0 :
        raise ValueError('%s < 0' % x)
    # Need better stopping condition than k == 100, and for large x
    #     should instead compute as 1 - upper incomplete gamma...
    return sum((-1) ** k / factorial(k) * x ** (s + k) / (s + k)
               for k in range(100)) / gamma(s)

def pchisq(x, k) :
    "Cumulative distribution function of chi squared with k degrees of freedom."
    if k < 1 or k != int(k) :
        raise ValueError('%s is not a positive integer' % k)
    if x < 0 :
        raise ValueError('%s < 0' % x)
    return gammainc(k / 2, x / 2)


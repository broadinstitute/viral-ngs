'''A few pure-python statistical tools to avoid the need to install scipy. '''
from __future__ import division # Division of integers with / should never round!
from math import exp, log, pi, sqrt
import itertools

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

def fisher_exact(contingencyTable) :
    """ Fisher exact test for the 2 x n case.
        contingencyTable is a sequence of 2 length-n sequences of integers.
        Return the two-tailed p-value against the null hypothesis that the row
            and column criteria are independent, using Fisher's exact test.
        For n larger than 2, this is very slow, O(S^(n-1)), where S is the
            smaller of the two row sums. Better to use chi2_contingency unless
            one of the row sums is small.
        Handles m x n contingencyTable with m > 2 if it can be reduced to the
            2 x n case by transposing or by removing rows that are all 0s. Also
            handles degenerate cases of 0 or 1 row by returning 1.0.
    """
    if len(contingencyTable) == 0 :
        return 1.0
    if len(set(map(len, contingencyTable))) != 1 :
        raise ValueError('Not all rows have the same length')
    if any(x != int(x) for row in contingencyTable for x in row) :
        raise ValueError('Some table entry is not an integer')

    # Eliminate rows and columns with 0 sum
    colSums = [sum(row[col] for row in contingencyTable)
               for col in range(len(contingencyTable[0]))]
    table = [[x for x, colSum in zip(row, colSums) if colSum != 0]
             for row in contingencyTable if sum(row) != 0]

    if len(table) < 2 or len(table[0]) < 2 :
        return 1.0

    if len(table) > len(table[0]) :
        table = list(zip(*table)) # Transpose

    m = len(table)
    n = len(table[0])
    if m != 2 :
        raise NotImplementedError('More than 2 non-zero rows and columns.')

    # Put row with smaller sum first. Makes the loop iterations simpler.
    table.sort(key = sum)
    # Put column with largest sum last. Makes loop quick rejection faster.
    table = list(zip(*table)) # Transpose
    table.sort(key = sum)
    table = list(zip(*table)) # Transpose back

    # There are many optimizations possible for the following code, but it would
    #     still be O(S^(n-1)) so it would still be too slow for anything
    #     sizeable, and it's usable as it for small things.

    rowSums = [sum(row) for row in table]
    colSums = [sum(row[col] for row in table) for col in range(n)]

    logChooseNrowSum = log_choose(sum(rowSums), rowSums[0])

    def prob_of_table(firstRow) :
        return exp(sum(log_choose(cs, a) for cs, a in zip(colSums, firstRow)) -
                   logChooseNrowSum)

    p0 = prob_of_table(table[0])
    result = 0
    for firstRowM1 in itertools.product(*[range(min(rowSums[0], colSums[i]) + 1)
                                         for i in range(n - 1)]) :
        lastElmt = rowSums[0] - sum(firstRowM1)
        if lastElmt < 0 or lastElmt > colSums[-1] :
            continue
        prob = prob_of_table(firstRowM1 + (lastElmt,))
        if prob <= p0 + 1e-9 : # (1e-9 handles floating point round off)
            result += prob

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
        # stirling is more accurate for larger values, so call it for a
        # larger value of s and use gamma recursion to get back to lower value.
        return exp(log_stirling(s + 9)) / product(s + i for i in range(10))

def gammainc(s, x) :
    """ Lower incomplete gamma function = 
            integral from 0 to x of t ** (s-1) exp(-t) dt divided by gamma(s),
        i.e., the fraction of gamma that you get if you integrate only until
            x instead of all the way to infinity.
        Implemented only for s > 0.
    """
    # scipy equivalent: scipy.special.gammainc(s,x)
    if s <= 0 :
        raise ValueError('%s is not positive' % s)
    if x < 0 :
        raise ValueError('%s < 0' % x)
    
    # Handle integers analytically
    if s == int(s) :
        term = 1
        total = 1
        for k in range(1, int(s)) :
            term *= x /  k
            total += term
        return 1 - exp(-x) * total

    # Otherwise use infinite series:
    # gammainc(s,x) = x ** s * exp(-x) / s / gamma(s) *
    #                 sum_k=0_to_infinity(x ** k / product_j=1_to_k(s + j))
    # which follows from the recursion formula:
    # gammainc(s, x) = gammainc(s - 1, x) - x ** (s - 1) * exp(-x) / gamma(s)
    absTol = 1e-9
    factor = x ** s * exp(-x) / s / gamma(s)
    term = 1
    total = 1
    for k in itertools.count(1) :
        term *= x / (s + k)
        if s + k + 1 > x :
            # If this will be the last term, we will sum a geometric progression
            # that approximates all remaining terms.
            ratio = x / (s + k + 1)
            potentialLastTerm = term / (1 - ratio)
            # Estimate error (1st two terms of difference between real series
            # and geometric progression are 0 so gain a factor of ratio ** 2)
            err = factor * potentialLastTerm * ratio ** 2
            if err < absTol :
                total += potentialLastTerm
                break
        total += term
    return factor * total

def pchisq(x, k) :
    "Cumulative distribution function of chi squared with k degrees of freedom."
    if k < 1 or k != int(k) :
        raise ValueError('%s is not a positive integer' % k)
    if x < 0 :
        raise ValueError('%s < 0' % x)
    return gammainc(k / 2, x / 2)


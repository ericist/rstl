import math
import numpy


class STL:
    """Python port of R's stl function.

    For more information see https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stl.html

    Args:
        ts: The time series (numpy array).
        freq: The frequency of the time series.
        s_window: Either the character string "periodic" or the span (in lags) of the loess window
        for seasonal extraction, which should be odd and at least 7, according to Cleveland et al.

    Optional args:
        s_degree, t_window, t_degree, l_window, l_degree, s_jump, t_jump, l_jump, robust, inner and outer.

        See https://stat.ethz.ch/R-manual/R-devel/library/stats/html/stl.html for their meanings and defaults.

    Attributes:
        trend: The trend component of the time series (numpy array).
        seasonal: The seasonal component of the time series (numpy array).
        remainder: The remainder of the time series not explained by trend and seasonal components (numpy array).
        weights: The final robust weights (all one if fitting is not done robustly) (numpy array).

        s_window, s_degree, t_window, t_degree, l_window, l_degree, s_jump, t_jump, l_jump, inner and outer.
        Note that these may have been altered by the program.
    """

    def __init__(self, ts, freq, s_window, s_degree=0, t_window=None,
                 t_degree=1, l_window=None, l_degree=None, s_jump=None,
                 t_jump=None, l_jump=None, robust=False, inner=None, outer=None):

        if ts.ndim != 1:
            raise ValueError("The time series must be 1-dimensional.")
        if numpy.isnan(numpy.sum(ts)):
            raise ValueError("The time series contains NaNs.")

        n = len(ts)

        if freq < 2:
            raise ValueError("The frequency must be greater than 1.")
        if n <= 2*freq:
            raise ValueError("The time series must contain more than 2 full periods of data.")

        if s_window == 'periodic':
            s_window = 10*n + 1
        if s_jump is None:
            s_jump = math.ceil(s_window/10)

        if t_window is None:
            t_window = nextodd(math.ceil(1.5*freq / (1 - 1.5/s_window)))
        if t_jump is None:
            t_jump = math.ceil(t_window/10)

        if l_window is None:
            l_window = nextodd(freq)
        if l_degree is None:
            l_degree = t_degree
        if l_jump is None:
            l_jump = math.ceil(l_window/10)

        if inner is None:
            inner = 1 if robust else 2
        if outer is None:
            outer = 15 if robust else 0

        weights = numpy.zeros(n)
        seasonal = numpy.zeros(n)
        trend = numpy.zeros(n)
        work = numpy.zeros((n+2*freq, 5))

        s_window = max(3, s_window)
        t_window = max(3, t_window)
        l_window = max(3, l_window)
        if s_window % 2 == 0:
            s_window += 1
        if t_window % 2 == 0:
            t_window += 1
        if l_window % 2 == 0:
            l_window += 1

        userw = False
        stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work)

        userw = True
        for _ in range(outer):
            work[:n, 0] = trend + seasonal
            stlrwt(ts, n, work[:n, 0], weights)
            stlstp(ts, n, freq, s_window, t_window, l_window, s_degree, t_degree, l_degree, s_jump, t_jump, l_jump, inner, userw, weights, seasonal, trend, work)

        if outer <= 0:
            weights.fill(1)

        self.seasonal = seasonal
        self.trend = trend
        self.remainder = ts - trend - seasonal
        self.weights = weights

        self.s_window = s_window
        self.t_window = t_window
        self.l_window = l_window
        self.s_degree = s_degree
        self.t_degree = t_degree
        self.l_degree = l_degree
        self.s_jump = s_jump
        self.t_jump = t_jump
        self.l_jump = l_jump
        self.inner = inner
        self.outer = outer


def nextodd(x):
    x = round(x)
    if x % 2 == 0:
        x += 1
    return x


def stless(y, n, length, ideg, njump, userw, rw, ys, res):
    if n < 2:
        ys[0] = y[0]
        return

    newnj = min(njump, n-1)
    if length >= n:
        nleft = 1
        nright = n
        for i in range(0, n, newnj):
            nys = stlest(y, n, length, ideg, i+1, ys[i], nleft, nright, res, userw, rw)
            if nys is not None:
                ys[i] = nys
            else:
                ys[i] = y[i]
    else:
        if newnj == 1:
            nsh = int((length+1)/2)
            nleft = 1
            nright = length
            for i in range(n):
                if (i+1) > nsh and nright != n:
                    nleft += 1
                    nright += 1
                nys = stlest(y, n, length, ideg, i+1, ys[i], nleft, nright, res, userw, rw)
                if nys is not None:
                    ys[i] = nys
                else:
                    ys[i] = y[i]
        else:
            nsh = int((length+1)/2)
            for i in range(1, n+1, newnj):
                if i < nsh:
                    nleft = 1
                    nright = length
                elif i >= (n-nsh+1):
                    nleft = n-length+1
                    nright = n
                else:
                    nleft = i-nsh+1
                    nright = length+i-nsh
                nys = stlest(y, n, length, ideg, i, ys[i-1], nleft, nright, res, userw, rw)
                if nys is not None:
                    ys[i-1] = nys
                else:
                    ys[i-1] = y[i-1]

    if newnj != 1:
        for i in range(0, n-newnj, newnj):
            delta = (ys[i+newnj] - ys[i]) / newnj
            ys[(i+1):(i+newnj)] = ys[i] + delta*numpy.arange(1, newnj)
        k = int(((n-1) // newnj)*newnj+1)

        if k != n:
            nys = stlest(y, n, length, ideg, n, ys[n-1], nleft, nright, res, userw, rw)
            if nys is not None:
                ys[n-1] = nys
            else:
                ys[n-1] = y[n-1]

            if k != (n-1):
                delta = (ys[n-1] - ys[k-1]) / (n-k)
                ys[k:(n-1)] = ys[k-1] + delta*numpy.arange(1, n-k)


def stlest(y, n, length, ideg, xs, ys, nleft, nright, w, userw, rw):
    nleft = int(nleft)
    nright = int(nright)

    h = max(xs-nleft, nright-xs)
    if length > n:
        h += (length-n) // 2

    r = numpy.abs(numpy.arange(nleft-xs, nright-xs+1))
    window = numpy.arange(nleft-1, nright)

    low_mask = r <= 0.001*h
    high_mask = r > 0.999*h
    mid_mask = numpy.logical_not(numpy.logical_or(low_mask, high_mask))
    lowmid_mask = numpy.logical_not(high_mask)

    low = window[low_mask]
    high = window[high_mask]
    mid = window[mid_mask]
    lowmid = window[lowmid_mask]

    w[low] = 1
    w[mid] = numpy.power(1 - numpy.power(r[mid_mask]/h, 3), 3)
    if userw:
        w[lowmid] *= rw[lowmid]
    a = numpy.sum(w[lowmid])

    w[high] = 0

    if a <= 0:
        ret = None
    else:
        w[(nleft-1):nright] /= a
        if h > 0 and ideg > 0:
            a = numpy.sum(w[(nleft-1):nright] * numpy.arange(nleft, nright+1))
            b = xs-a
            c = numpy.sum(w[(nleft-1):nright] * numpy.square(numpy.arange(nleft-a, nright-a+1)))
            if math.sqrt(c) > 0.001*(n-1):
                b /= c
                w[(nleft-1):nright] *= (b * numpy.arange(nleft-a, nright+1-a) + 1)
        ret = numpy.sum(w[(nleft-1):nright] * y[(nleft-1):nright])

    return ret


def stlfts(x, n, np, trend, work):
    stlma(x, n, np, trend)
    stlma(trend, n-np+1, np, work)
    stlma(work, n-2*np+2, 3, trend)


def stlma(x, n, length, ave):
    v = numpy.sum(x[:length])
    ave[0] = v / length

    newn = n - length + 1
    if newn > 1:
        k = length
        m = 0
        for j in range(1, newn):
            k += 1
            m += 1
            v = v - x[m-1] + x[k-1]
            ave[j] = v / length


def stlstp(y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, work):
    for _ in range(ni):
        work[:n, 0] = y - trend
        stlss(work[:, 0], n, np, ns, isdeg, nsjump, userw, rw, work[:, 1], work[:, 2], work[:, 3], work[:, 4], season)
        stlfts(work[:, 1], n+2*np, np, work[:, 2], work[:, 0])
        stless(work[:, 2], n, nl, ildeg, nljump, False, work[:, 3], work[:, 0], work[:, 4])
        season[:] = work[np:np+n, 1] - work[:n, 0]
        work[:n, 0] = y - season
        stless(work[:, 0], n, nt, itdeg, ntjump, userw, rw, trend, work[:, 2])


def stlrwt(y, n, fit, rw):
    r = numpy.abs(y - fit)

    med = 6*numpy.median(r)
    low = r <= 0.001*med
    high = r > 0.999*med
    mid = numpy.logical_not(numpy.logical_or(low, high))

    rw[low] = 1
    rw[mid] = numpy.square(1 - numpy.square(r[mid] / med))
    rw[high] = 0


def stlss(y, n, np, ns, isdeg, nsjump, userw, rw, season, work1, work2, work3, work4):
    for j in range(np):
        k = (n-j-1) // np + 1
        work1[:k] = y[numpy.arange(k)*np+j]

        if userw:
            work3[:k] = rw[numpy.arange(k)*np+j]

        stless(work1, k, ns, isdeg, nsjump, userw, work3, work2[1:], work4)
        nright = min(ns, k)

        nval = stlest(work1, k, ns, isdeg, 0, work2[0], 1, nright, work4, userw, work3)
        if nval is not None:
            work2[0] = nval
        else:
            work2[0] = work2[1]

        nleft = max(1, k-ns+1)

        nval = stlest(work1, k, ns, isdeg, k+1, work2[k+1], nleft, k, work4, userw, work3)
        if nval is not None:
            work2[k+1] = nval
        else:
            work2[k+1] = work2[k]

        for m in range(k+2):
            season[m*np+j] = work2[m]

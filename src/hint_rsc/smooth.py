"""Residual smoothing filters, applied per chromosome (identical to the paper's implementation)."""
import math

import numpy as np

METHODS = ('mean', 'gaussian', 'median', 'savgol', 'haar', 'sym4')


def mean_filter(x, w):
    """Mean over the w bins centred on each bin; the window is truncated at chromosome ends."""
    n = len(x)
    out = np.zeros(n)
    for i in range(n):
        s = max(0, i - math.floor(w / 2))
        e = min(i + math.floor(w / 2), n)
        out[i] = np.mean(x[s:e + 1])
    return out


def smooth(x, method, p):
    """p is the window width in bins (mean, gaussian, median, savgol) or the number of detail levels (haar, sym4).

    gaussian: sigma = p / sqrt(12), the standard deviation of a width-p box; savgol: local quadratic fit over p bins
    (p odd, >= 5); haar, sym4: soft universal threshold on detail levels 1..p. p = 0 returns x unchanged.
    """
    if not p:
        return x
    if method == 'mean':
        return mean_filter(x, p)
    if method == 'gaussian':
        from scipy.ndimage import gaussian_filter1d
        return gaussian_filter1d(x, sigma=p / math.sqrt(12), mode='nearest', truncate=4.0)
    if method == 'median':
        from scipy.ndimage import median_filter
        return median_filter(x, size=p, mode='nearest')
    if method == 'savgol':
        from scipy.signal import savgol_filter
        return savgol_filter(x, p, 2, mode='nearest')
    if method in ('haar', 'sym4'):
        import pywt
        coeffs = pywt.wavedec(x, method, level=p, mode='symmetric')
        sigma = np.median(np.abs(coeffs[-1])) / 0.6745
        thr = sigma * math.sqrt(2 * math.log(len(x)))
        if thr > 0:                                            # thr = 0 (no noise): nothing to remove
            coeffs = [coeffs[0]] + [pywt.threshold(c, thr, 'soft') for c in coeffs[1:]]
        return pywt.waverec(coeffs, method, mode='symmetric')[:len(x)]
    raise ValueError(f'unknown smoothing method {method!r}; choose from {", ".join(METHODS)}')


def smooth_dict(resid, method, p):
    return {c: smooth(v, method, p) for c, v in resid.items()}

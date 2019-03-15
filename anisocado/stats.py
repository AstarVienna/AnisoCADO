import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from anisocado import AnalyticalScaoPsf


def strehl_map(r=25, dr=3, **kwargs):
    psf = AnalyticalScaoPsf(**kwargs)
    x, y = np.mgrid[-r:r+dr:dr, -r:dr+r:dr]

    strlmap = np.zeros(x.shape)
    for i in range(len(x)):
        for j in range(len(y)):
            psf.shift_psf_off_axis(x[i, j], y[i, j])
            strlmap[i, j] = psf.strehl_ratio

    return strlmap


def on_axis_strehl_for_kernel_size(Narr=(128, 512, 2048), **kwargs):
    """Only for the on-axis kernel"""
    return [AnalyticalScaoPsf(N=N, **kwargs).strehl_ratio for N in Narr]

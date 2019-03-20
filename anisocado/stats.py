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
            psf.shift_off_axis(x[i, j], y[i, j])
            strlmap[i, j] = psf.strehl_ratio

    return strlmap


def on_axis_strehl_for_kernel_size(Narr=(128, 512, 2048), **kwargs):
    """Only for the on-axis kernel"""
    return [AnalyticalScaoPsf(N=N, **kwargs).strehl_ratio for N in Narr]


def make_psf_grid(r=14, dr=7, **kwargs):
    psf = AnalyticalScaoPsf(**kwargs)
    x, y = np.mgrid[-r:r+1:dr, -r:r+1:dr]

    psf_grid = []
    for i in range(len(x)):
        for j in range(len(y)):
            psf.shift_off_axis(x[i, j], y[i, j])
            psf_grid += [psf.kernel]

    return psf_grid


def make_image_of_psf_grid():
    psf_grid = make_psf_grid(wavelength=2.15, N=128)

    plt.figure(figsize=(10, 10))
    i = 0
    for y in range(5):
        for x in range(5):
            plt.subplot(5, 5, 1+x+5*(4-y))
            plt.imshow(psf_grid[i], origin="l", norm=LogNorm())
            plt.axis("off")
            plt.title("({}, {})".format((7*x-14), (7*x-14)))
            i += 1
    plt.suptitle("Ks-band (2.15um) SCAO PSFs")
    plt.show()

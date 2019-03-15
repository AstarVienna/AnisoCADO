import pytest

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from anisocado import AnalyticalScaoPsf

PLOTS = False


class TestMakeStrehlMap:
    def test_make_a_strehl_map_over_a_grid_on_the_focal_plane(self):
        psf = AnalyticalScaoPsf(N=128, wavelengthIR=2.e-6)
        x, y = np.mgrid[-25:26:1, -25:26:1]

        sr = np.zeros(x.shape)
        for i in range(len(x)):
            for j in range(len(y)):
                psf.shift_psf_off_axis(x[i, j], y[i, j])
                sr[i, j] = psf.strehl_ratio

        if PLOTS:
            plt.imshow(sr, norm=LogNorm())
            plt.colorbar()
            plt.show()


class TestStrehlRatiosForDifferentSizedKernels:
    def test_how_much_does_the_SR_change_for_kernel_size_increase(self):
        from anisocado.stats import on_axis_strehl_for_kernel_size

        side_length = [64, 128, 256, 512, 1024]
        sr = on_axis_strehl_for_kernel_size(side_length, wavelengthIR=2.15e-6)
        if PLOTS:
            plt.plot(side_length, sr)
            plt.show()

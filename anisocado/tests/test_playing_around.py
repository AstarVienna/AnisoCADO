import pytest

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from anisocado import AnalyticalScaoPsf
from anisocado.misc import on_axis_strehl_for_kernel_size

PLOTS = False


class TestMakeStrehlMap:
    def test_make_a_strehl_map_over_a_grid_on_the_focal_plane(self):
        psf = AnalyticalScaoPsf(N=128, wavelength=2.15)
        x, y = np.mgrid[-50:51:20, -50:51:20]

        sr = np.zeros(x.shape)
        for i in range(len(x)):
            for j in range(len(y)):
                psf.shift_off_axis(x[i, j], y[i, j])
                sr[i, j] = psf.strehl_ratio
        print(np.max(sr))

        if PLOTS:
            plt.imshow(sr)  # norm=LogNorm())
            plt.contourf(x, y, sr, np.arange(0, 1, 0.05))
            plt.colorbar()
            plt.show()


class TestStrehlRatiosForDifferentSizedKernels:
    def test_how_much_does_the_SR_change_for_kernel_size_increase(self):
        side_length = [64, 128, 256, 512, 1024]
        sr = on_axis_strehl_for_kernel_size(side_length, wavelength=2.15)
        if PLOTS:
            plt.plot(side_length, sr)
            plt.show()


class TestPsfsAlongXaxis:
    def test_strehl_ratio_along_x_axis(self):
        for wave in [0.8, 0.9, 1.1, 1.6, 2.15]:
            psf = AnalyticalScaoPsf(N=128, wavelength=wave)
            x, y = np.mgrid[0:41:1, 0:1:1]

            sr = np.zeros(x.shape)
            for i in range(x.shape[0]):
                for j in range(x.shape[1]):
                    psf.shift_off_axis(x[i, j], y[i, j])
                    sr[i, j] = psf.strehl_ratio

            if PLOTS:
                print(wave, np.sum(psf._kernel_sum))
                plt.plot(x.flatten(), sr, label="{} um".format(wave))

        if PLOTS:
            # plt.semilogy()
            plt.legend()
            plt.xlabel("Distance from NGS [arcsec]")
            plt.ylabel("Strehl Ratio")
            plt.show()


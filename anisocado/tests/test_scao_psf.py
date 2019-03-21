import pytest

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits

from anisocado.psf import AnalyticalScaoPsf

PLOTS = False


class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(AnalyticalScaoPsf(), AnalyticalScaoPsf)

    def test_initialises_with_bunch_of_keywords(self):
        psf = AnalyticalScaoPsf(N=256, wavelengthIR=1.2e-6)
        assert isinstance(psf, AnalyticalScaoPsf)

    def test_throws_warning_if_keyword_not_spelt_correctly(self):
        psf = AnalyticalScaoPsf(N=256, wavelengthIR=1.2e-6, hello="world")
        assert isinstance(psf, AnalyticalScaoPsf)

    def test_on_axis_psf_is_made(self):
        psf = AnalyticalScaoPsf(N=512)
        assert isinstance(psf.psf_on_axis, np.ndarray)

        if PLOTS:
            plt.imshow(psf.psf_on_axis.T, origin="lower", norm=LogNorm(), vmin=1e-7)
            plt.show()

    def test_kernel_sums_to_one(self):
        psf = AnalyticalScaoPsf(N=256)
        assert np.sum(psf.psf_latest) == pytest.approx(1)


class TestShiftPSF:
    def test_psf_blurs_when_shifted(self):
        psf = AnalyticalScaoPsf(N=1024)
        sr_orig = psf.strehl_ratio

        psf.shift_off_axis(0, 10)
        sr_last = psf.strehl_ratio

        assert sr_last < sr_orig

        if PLOTS:
            plt.subplot(121)
            plt.imshow(psf.psf_on_axis.T, origin="lower", norm=LogNorm(), vmin=3e-6)
            plt.colorbar()

            plt.subplot(122)
            plt.imshow(psf.psf_latest.T, origin="lower", norm=LogNorm(), vmin=3e-6)
            plt.colorbar()
            plt.show()

    def test_kernel_sums_to_one(self):
        psf = AnalyticalScaoPsf(N=512)
        psf.shift_off_axis(0, 0)
        assert np.sum(psf.psf_latest) == pytest.approx(1)


class TestHDUProperty:
    def test_returns_fits_imagehdu(self):
        psf = AnalyticalScaoPsf(N=512)
        psf.shift_off_axis(10, 10)
        hdu = psf.hdu
        assert isinstance(hdu, fits.ImageHDU)
        assert hdu.header["CDELT1"] == 4 / (3600 * 1000.)
        assert hdu.header["CRVAL1"] == 10 / 3600.

        if PLOTS:
            plt.imshow(hdu.data.T, origin="l", norm=LogNorm())
            plt.show()




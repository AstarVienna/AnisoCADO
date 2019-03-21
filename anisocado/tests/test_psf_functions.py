import pytest

import numpy as np
from astropy.io import fits

from anisocado import misc

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


@pytest.fixture(scope="class")
def basic_fv_psf():
    n = 30
    coords = [(-n, n),  (0, n),  (n, n),
              (-n, 0),  (0, 0),  (n, 0),
              (-n, -n), (0, -n), (n, -n)]
    waves = [0.9, 1.1, 1.6, 2.15]
    hdu = misc.make_simcado_psf_file(coords=coords, wavelengths=waves, N=128)

    return hdu


class TestSimcadoPsfFile:
    def test_throws_error_with_no_input(self):
        with pytest.raises(TypeError):
            misc.make_simcado_psf_file()

    def test_returns_hdulist_for_basic_input(self):
        hdu = misc.make_simcado_psf_file(coords=[(0, 0)], wavelengths=[2.15],
                                                   N=128)
        assert isinstance(hdu, fits.HDUList)
        assert isinstance(hdu[1], fits.BinTableHDU)
        assert len(hdu) == 3
        assert hdu[2].data.shape == (1, 128, 128)

    def test_returns_full_hdulist_for_full_input(self):
        n = 30
        coords = [(0, 0), (-n, -n), (n, -n), (n, n), (-n, n)]
        waves = [0.9, 1.1, 1.6, 2.15]
        hdu = misc.make_simcado_psf_file(coords=coords, wavelengths=waves,
                                         N=128)
        assert len(hdu) == len(waves) + 2
        assert hdu[2].data.shape == (len(coords), 128, 128)
        assert hdu[2].header["CDELT1"] == 0.004 / 3600.


@pytest.mark.usefixtures("basic_fv_psf")
class TestFvpsfFileConsistencyChecks:
    def test_has_correct_keywords_in_ext0(self, basic_fv_psf):
        for key in ["AUTHOR", "DATE_CRE", "DATE_MOD", "SOURCE", "STATUS",
                    "ETYPE", "ECAT", "EDATA"]:
            assert key in basic_fv_psf[0].header

    def test_has_catalogue_psf(self, basic_fv_psf):
        ecat = basic_fv_psf[0].header["ECAT"]
        cat_hdu = basic_fv_psf[ecat]
        for key in ["NUMPSFS", "CATTYPE", "CUNIT1"]:
            assert key in cat_hdu.header

        if cat_hdu.header["CATTYPE"] == "table":
            assert isinstance(cat_hdu.data, (fits.BinTableHDU, fits.FITS_rec))

        elif basic_fv_psf[ecat].header["CATTYPE"] == "image":
            assert isinstance(cat_hdu.data, fits.ImageHDU)
            for key in ["CRVAL1", "CRPIX1", "CDELT1"]:
                assert key in cat_hdu.header

    def test_has_the_right_number_of_layers_per_psf_hdu(self, basic_fv_psf):
        ecat = basic_fv_psf[0].header["ECAT"]
        cat_hdu = basic_fv_psf[ecat]
        n_layers = len(set(cat_hdu.data["layer"]))

        edata = basic_fv_psf[0].header["EDATA"]
        data_indexes = range(edata, len(basic_fv_psf))
        for ii in data_indexes:
            assert basic_fv_psf[ii].data.shape[0] == n_layers
            for key in ["CRVAL1", "CRPIX1", "CDELT1", "WAVE0"]:
                assert key in basic_fv_psf[ii].header

    def plot_psfs(self, basic_fv_psf):
        for psf_hdu in basic_fv_psf[2:]:
            for i in range(len(psf_hdu.data)):
                plt.subplot(3, 3, i + 1)
                plt.imshow(psf_hdu.data[i, :, :].T, origin="l", norm=LogNorm())

            plt.show()

    def print_coords(self):
        #coords = psf.field_positions_for_simcado_psf()
        #waves = [0.8, 0.9, 1.05, 1.2, 1.5, 1.7, 2.0, 2.2]
        #hdu = psf.make_simcado_psf_file()

        map = misc.make_strehl_map_from_table([(0, 0)])
        plt.imshow(map)
        plt.show()
        print(np.unique(map))


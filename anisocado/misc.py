import numpy as np
from astropy.io import fits
from astropy.table import Table

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


def make_simcado_psf_file(coords, wavelengths, header_cards=None, **kwargs):
    """
    Generate a set of Field-Varying PSF cubes for use with SimCADO

    Parameters
    ----------
    coords : list of tuples
        [arcsec] Sample positions of the PSF in field of view. (0, 0) is the
        centre of the field of view.

    wavelengths : list
        [um] The wavelengths for which the PSF should be sampled

    header_cards : dict, optional
        Any extra keyword-value pair to be added to the extension 0 header.

    kwargs
    ------
    Keyword-value pairs accepted by an ``AnalyticalScaoPsf`` object

    Returns
    -------
    hdulist : fits.HDUList
        A HDUList object which is formatted for use as a Field-Varying PSF in
        SimCADO

    """
    # accept list of coordinates
    # accept list of wavelengths
    # make ext0
    # make table for ext1
    # initialise PSF
    # shift to each positions
    # make hdus for each position
    # build hdulist
    # return hdulist

    x, y = np.array(coords).T
    layers = np.arange(len(x))

    keys = ["AUTHOR", "DATE_CRE", "DATE_MOD", "SOURCE", "STATUS"]
    ext0_dict = {key: "" for key in keys}
    ext0_dict["ETYPE"] = "FVPSF"
    ext0_dict["ECAT"] = (1, "The extension containing the catalogue data")
    ext0_dict["EDATA"] = (2, "The first extension with real data")
    ext0_dict.update({"WAVEEXT{}".format(i + 2): w
                      for i, w in enumerate(wavelengths)})

    pri_hdr = fits.PrimaryHDU()
    pri_hdr.header.update(ext0_dict)
    pri_hdr.header.update(header_cards)

    ext1_dict = {"NUMPSFS": len(x), "CATTYPE": "table", "CUNIT1": "arcsec"}

    tbl = Table(data=[x, y, layers], names=["x", "y", "layer"])
    cat_hdu = fits.table_to_hdu(tbl)
    cat_hdu.header.update(ext1_dict)

    psf_hdus = []
    for wave in wavelengths:
        print("Making psf cube for {} um".format(wave))
        psf = AnalyticalScaoPsf(wavelength=wave, **kwargs)
        kernel_cube = [psf.shift_off_axis(dx, dy) for dx, dy in coords]

        psf_hdu = psf.hdu
        psf_hdu.data = np.array(kernel_cube)
        psf_hdu.header["WAVE0"] = wave
        psf_hdu.header["WAVEUNIT"] = "um"

        psf_hdus += [psf_hdu]

    hdulist = fits.HDUList([pri_hdr, cat_hdu] + psf_hdus)

    return hdulist


def field_positions_for_simcado_psf():
    coords = [(0, 0)]
    # for r in [7.5, 15, 22.5]:
    for r in [1, 2, 3, 5, 7, 11, 15, 23, 31]:
        for ang in np.arange(0, 360, 45):
            coords += [(r * np.cos(np.deg2rad(ang)),
                        r * np.sin(np.deg2rad(ang)))]

    return coords


def make_strehl_map_from_coords(coords):
    x, y = np.array(coords).T

    from scipy.interpolate import griddata
    map = griddata((x, y), np.arange(len(x)),
                   np.array(np.meshgrid(np.arange(-25, 26),
                                        np.arange(-25, 26))).T,
                   method="nearest")

    return map

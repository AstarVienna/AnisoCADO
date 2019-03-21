"""Originally from file _anisocado.py"""

import numpy
import numpy as np
from . import pupil_utils

#  ____  _____    _    ____  __  __ _____
# |  _ \| ____|  / \  |  _ \|  \/  | ____|
# | |_) |  _|   / _ \ | | | | |\/| |  _|
# |  _ <| |___ / ___ \| |_| | |  | | |___
# |_| \_\_____/_/   \_\____/|_|  |_|_____|

"""
Hello.
This files contains some useful functions related to psf generation.
They are just raw, so that you can insert them into your own classes
as desired.
The functions are commented to help you understand what's going on, when
you need to modify them.

Don't read that code anyway.
Go right now to [the file _anisocado.py].
It contains examples that will show you how to use the functions, what they do,
etc.
"""


def defineDmFrequencyArea(kx, ky, rotdegree, dactu=0.5403):
    """
    <kx>        : spatial frequency produced by the
                  function computeSpatialFreqArrays() (metres^-1)
    <ky>        : idem kx
    <rotdegree> : rotation of M4 (degrees)
    <dactu>     : value of actuator pitch of M4 (metres)

    The function returns the 2D domain of spatial frequencies which can be
    impacted by M4, i.e. the M4 compensation domain.
    Underlying assumtions are than M4 is of an hexagonal pattern, with an
    inter-actuator distance <dactu> set by default to 54.03 cm.

    Result is returned with a frequency corner-centred representation.

    Standalone example::

        N = 512               # output will be 512x512
        pixelSize = 4.2       # 4.2 mas
        wavelength = 1.65e-6  # metres
        kx, ky = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, 0)
        plt.imshow(np.fft.fftshift(M4).T, origin='l')

    """

    # cut-off frequency
    fc = (1./np.sqrt(3)) / dactu

    # mask frequency definition
    A = np.pi/3  # 60 degrees
    A0 = rotdegree * np.pi / 180
    msk = np.abs(np.cos(A0)*ky+np.sin(A0)*kx)<fc
    msk = np.logical_and(msk, np.abs(np.cos(A+A0)*ky+np.sin(A+A0)*kx)<fc)
    msk = np.logical_and(msk, np.abs(np.cos(2*A+A0)*ky+np.sin(2*A+A0)*kx)<fc)
    k = np.sqrt(kx**2 + ky**2)
    msk = np.logical_and(msk, k < (fc*1.115))

    return msk


def computeSpatialFreqArrays(N, pixelSize, wavelength):
    """
    <N>             : size of the output image (NxN)
    <pixelSize>     : size of the pixels of the psf image (mas)
    <wavelength>    : wavelength (metres)

    The function returns a tuple of spatial frequencies (m^-1) in X and Y
    together with the pixel scale of these (the value will be useful later
    on in other procedures).
    Results of 2D arrays are returned with a Fourier corner-centred
    representation.

    Standalone example::

        N = 512               # output will be 512x512
        pixelSize = 4.2       # 4.2 mas
        wavelength = 1.65e-6  # metres
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)

    """

    # array of indices centred 'as expected' by Fourier frequencies, in 1D
    k1d = np.fft.fftshift(np.arange(N) - (N//2))

    # Proper scaling to transform indices in spatial frequencies.

    # dX = pixelSize * 4.84813681109536e-09  # from mas to radians
    dX = pixelSize * 4.84813681109536e-06  # from arcsec to radians
    uk = dX / wavelength
    k1d = k1d * uk  # now this is a spatial freq in metres^-1

    # now creating 2D arrays of spatial frequency
    kx, ky = np.meshgrid(k1d, k1d, indexing='ij') # for convention [x,y]

    return kx, ky, uk


def computeWiener(kx, ky, L0, r0):
    """
    <kx>  : spatial frequency produced by
            the function computeSpatialFreqArrays() (metres^-1)
    <ky>  : idem kx
    <L0>  : value of the outer scale (metres)
    <r0>  : value of Fried parameter r0 (metres)

    The function returns the 2D spectrum of the turbulence with an outer
    scale L0. It is expressed in rad^2.m^2.

    Result is returned with a frequency corner-centred representation.

    Standalone example::

        N = 512               # output will be 512x512
        pixelSize = 4.2       # 4.2 mas
        wavelength = 1.65e-6  # metres.  H band.
        L0 = 25.              # 25 metres outer scale
        r0 = 0.6              # r0 = 60cm in H band
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        W = computeWiener(kx, ky, L0, r0)

    """

    # computation of Wiener spectrum expressed in radians^2 (at the wavelength
    # where r0(lambda) is expressed !)
    Wiener = (kx**2 + ky**2 + 1./L0**2.)**(-11./6)
    Wiener *= 0.0228956 * r0**(-5./3)

    # frequency 0 set to 0. It's wrong anyway.
    Wiener[0, 0] = 0.

    return Wiener


def anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy, wavelength,
                          kx, ky, Wiener, M4):
    """
    <Cn2h>          : list of normalised (np.sum(Cn2h)==1.0) strengh of
                      turbulence of each layer (no unit)
    <layerAltitude> : list of altitudes of the turbulence layers (metres)
    <L0>            : value of the outer scale (metres)
    <offx>          : off-axis distance of the star along X axis (arcsec)
    <offy>          : idem along Y axis (arcsec)
    <wavelength>    : wavelength (metres)
    <kx>            : spatial frequency produced by
                      the function computeSpatialFreqArrays() (metres^-1)
    <ky>            : idem kx
    <Wiener>        : turbulent spectrum from function computeWiener()
    <M4>            : frequency domain of M4 (boolean) computed by the
                      function M4 = defineDmFrequencyArea(kx, ky, 0)

    Computes the phase spectrum (spatial frequencies) contribution due to the
    anisoplanatism in the AO compensation for a source located off-axis from
    the SCAO guide star by some amount (offx, offy).
    In input, the distribution of the turbulence in altitude is given.

    Result is returned with a frequency corner-centred representation.

    Example::

        N = 512
        pixelSize = 4.2
        Cn2h = [0.3, 0.2, 0.2, 0.1]
        layerAltitude = [0,1000,8000,15000]
        L0 = 25.0
        offx, offy = (34., 12.)
        wavelength = 1.65e-6
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, 0)
        W = computeWiener(kx, ky, L0, 1.0)
        f = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy,
                                  wavelength, kx, ky, W, M4)

    """

    # number of turbulent layers involved in that computation
    nlayers = len(Cn2h)

    # conversion arcsec to radians
    RASC = 206264.8062471

    # prepare memory alloc for transfer function of anisoplanatism
    Haniso = np.zeros(kx.shape)

    # loop over turbulent layers, summing transfer function of each layer
    for i in range(nlayers):
        dx = layerAltitude[i] * offx / RASC # shift in metres on the layer in X
        dy = layerAltitude[i] * offy / RASC # idem, in Y
        tmp = (2j*np.pi*dx)*kx + (2j*np.pi*dy)*ky
        Haniso += Cn2h[i] * np.abs(1 - np.exp( tmp ))**2

    # now applying the transfer function on the Wiener spectrum, only in the
    # spatial frequency range of M4
    Waniso = np.zeros(Haniso.shape)
    Waniso[M4] = Haniso[M4] * Wiener[M4]

    return Waniso


def fittingSpectrum(Wiener, M4):
    """
    <Wiener>        : turbulent spectrum from function computeWiener()
    <M4>            : frequency domain of M4 (boolean) computed by the
                      function M4 = defineDmFrequencyArea(kx, ky, 0)

    Returns the spatial spectrum of the (so-called) fitting error, i.e. the
    residual phase after a full, perfect, ideal, instantaneous, super-duper,
    hyper-clean, theoretical compensation of M4. It is expressed in rad^2.m^2.
    Result is returned with a frequency corner-centred representation.

    Example::

        N = 512
        pixelSize = 4.2
        L0 = 25.0
        offx, offy = (34., 12.)
        wavelength = 1.65e-6
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, 0)
        r0 = 0.6
        W = computeWiener(kx, ky, L0, r0)
        f = fittingSpectrum(W, M4)

    """

    Wfit = Wiener.copy()
    Wfit[M4] = 0.0 # M4 cancels whatever is in its compensation domain
    return Wfit


def otherSpectrum(nmRms, M4, uk, wavelength):
    """
    <nmRms>  : number of nm rms
    <M4>     : frequency domain of M4 (boolean) computed by the
                      function M4 = defineDmFrequencyArea(kx, ky, 0)
    <uk> : # size of the 'spatial frequency pixel' in m^-1
    <wavelength>    : wavelength (metres)

    Example::

        N = 512
        pixelSize = 4.2
        L0 = 25.0
        offx, offy = (34., 12.)
        wavelength = 1.65e-6
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, 0)
        nmRms = 250.
        f = otherSpectrum(nmRms, M4, uk, wavelength)

    """

    fact = 2 * np.pi * nmRms * 1e-9 / uk / wavelength
    fact = fact**2
    tot = np.sum(M4)
    Wothers = np.zeros(M4.shape)
    Wothers[M4] = fact / tot

    return Wothers


def aliasingSpectrum(kx, ky, r0, L0, M4, dssp=0.4015):
    """
    Example::

        N = 512
        pixelSize = 4.2
        L0 = 25.0
        r0 = 0.6
        wavelength = 1.65e-6
        rotdegree = 10.
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, rotdegree)
        W = aliasingSpectrum(kx, ky, r0, L0, M4)

    """

    ke = 1.0 / dssp  # computes the sampling spatial-frequency of the WFS

    kxt = kx[M4]
    kyt = ky[M4]
    Wt  = ((kxt-ke)**2 + kyt**2 + 1./L0**2.)**(-11./6)
    Wt += ((kxt+ke)**2 + kyt**2 + 1./L0**2.)**(-11./6)
    Wt += (kxt**2 + (kyt-ke)**2 + 1./L0**2.)**(-11./6)
    Wt += (kxt**2 + (kyt+ke)**2 + 1./L0**2.)**(-11./6)
    Wt *= 0.0228956 * r0**(-5./3)
    Walias = np.zeros(M4.shape)
    Walias[M4] = Wt

    return Walias


def computeBpSpectrum(kx, ky, V, Fe, tret, gain, Wiener, M4):
    """
    Example::

        N = 512
        pixelSize = 4.2
        L0 = 25.0
        r0 = 0.6
        wavelength = 1.65e-6
        rotdegree = 10.
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, rotdegree)
        W = computeWiener(kx, ky, L0, r0)

        V = 10.
        Fe = 500.
        tret = 0.004
        gain = 0.3
        f = computeBpSpectrum(kx, ky, V, Fe, tret, gain, W, M4)

    """

    k = np.sqrt(kx*kx + ky*ky)
    nu = k * V / np.sqrt(2)    # pourquoi un sqrt(2) ? je ne saurais dire ...!!!
    Wbp = hcor(nu, Fe, tret, gain, 500) * Wiener
    Wbp[np.logical_not(M4)] = 0.
    return Wbp


def hcor(freq, Fe, tret, G, BP, an=True):
    """
    ***** Fonction extraite de STYC le 9 Jan 2013 *****

    The option an=1 sets the integrator to "analog". Doing this, an
    extra 1/2 frame delay is added compared to case of the numeric
    integrator an=0.

    <tret> is the delay expressed as a *time in seconds*, between the
    end of the integration and the start of the command.

    """

    Te = 1. / Fe
    p = 1j * 2 * np.pi * freq + 1e-12

    Hint = 1./(1-np.exp(-p*Te)) # numeric integrator
    Hccd = (1.-np.exp(-p*Te))/(p*Te)  # echant bloqueur avec retard 1/2 trame
    Hdac = Hccd                    # echant bloqueur avec retard 1/2 trame
    Hret = np.exp(-p*tret)
    Hmir = 1./(1. + 1j*freq/BP)
    Hbo = Hint * Hccd * Hdac * Hret * Hmir
    Hcor   = 1./abs(1 + Hbo*G)**2

    return Hcor


def convertSpectrum2Dphi(W, uk):
    """
    <W>  : spatial spectrum to be converted into phase structure function
           in rd^2.m^2
    <uk> : # size of the 'spatial frequency pixel' in m^-1

    Converts the 2D spectrum into a phase structure function Dphi.
    Uses Dphi(r) = $ $ (1-cos(2.pi.k.r)) W(k) d2k
    Computation of Dphi is in radians^2 at the wavelength of r0.
    """

    W[0, 0] = 0.0
    W[0, 0] = -np.sum(W)
    Dphi = 2*np.abs(np.fft.fft2(W)) *  (uk**2)
    return Dphi


def fake_generatePupil(N, deadSegments, rotdegree, pixelSize, wavelength):
    """
    <N>            : size of the output image, that is made to match the size
                     of the (square) psf image to be processed. In other
                     words, N = psf.shape[0]
    <deadSegments> : number of hexa segments of M1 that are missing
    <rotdegree>    : pupil rotation in degrees
    <pixelSize>    : size of the pixels of the psf image (mas)
    <wavelength>   : wavelength (metres)

    Examples::

        N = 512
        deadSegments = 3
        rotdegree = 14.
        pixelSize = 4.2
        wavelength = 1.65e-6
        pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize,
                                 wavelength)

    """
    nseg = pupil_utils.getEeltSegmentNumber()
    refl = np.ones(nseg)+np.random.randn(nseg)/20.
    if deadSegments:
        refl[(np.random.rand(deadSegments)*nseg).astype(int)] = 0.
    i0 = N/2+0.5
    j0 = N/2+0.5

    # field of view of the psf image in rd
    FoV = N * pixelSize * 4.84813681109536e-06  # from arcsec to radians
    # original line used mas
    # FoV = N * pixelSize * 4.84813681109536e-09  # from mas to radians

    # pixel scale of pupil image
    pixscale = wavelength / FoV   # expressed in metres
    dspider = 0.53
    gap = 0.02
    pup = pupil_utils.generateEeltPupilReflectivity(refl, N, dspider, i0, j0,
                                                    pixscale, gap, rotdegree,
                                                    softGap=True)
    return pup


def computeEeltOTF(pup):
    """
    """
    # Computation of telescope OTF
    Nx, Ny = pup.shape
    FTOtel = np.fft.fft2( np.abs(np.fft.fft2(pup))**2 ).real
    FTOtel /= np.sum(pup)**2 * Nx * Ny
    return FTOtel


def core_generatePsf(Dphi, FTOtel):
    """
    Examples
    --------
    ::
        N = 1024
        pixelSize = 4.2

        # atmospheric profile (old ESO profile before 2010)
        layerAltitude = [47., 140, 281, 562, 1125, 2250, 4500, 9000, 18000.]

        # from ref. E-SPE-ESO-276-0206_atmosphericparameters
        Cn2h = [52.24, 2.6, 4.44, 11.60, 9.89, 2.95, 5.98, 4.30, 6]
        Cn2h = np.array(Cn2h)
        Cn2h /= np.sum(Cn2h)

        # outer scale
        L0 = 25.0
        offx, offy = (15., 20.)
        wavelength = 1.65e-6
        rotdegree = 10.0
        kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
        M4 = defineDmFrequencyArea(kx, ky, rotdegree)
        r0 = 0.6

        # This is the turbulent spectrum ....
        W = computeWiener(kx, ky, L0, r0)

        # And here are some of the PSF-destroyers
        Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy,
                                       wavelength, kx, ky, W, M4)
        Wfit = fittingSpectrum(W, M4)
        nmRms = 100.
        Wother = otherSpectrum(nmRms, M4, uk, wavelength)
        Dphi = convertSpectrum2Dphi(Waniso + Wfit + Wother, uk)

        # Here, i need to generate a kind of PSF or telescope OTF
        deadSegments = 3
        pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize,
                                 wavelength)
        FTOtel = computeEeltOTF(pup)

        psf = core_generatePsf(Dphi, FTOtel)

        print(psf.max())
        window = 100
        plt.imshow( psf[N//2-window:N//2+window, N//2-window:N//2+window]**0.3 )

    """

    # total FTO
    FTO = np.exp(-0.5*Dphi) * FTOtel

    # PSF
    psf = np.fft.fftshift( np.fft.fft2(FTO).real )
    return psf


def createAdHocScaoPsf(N, pixelSize, wavelengthIR, rotdegree, r0Vis, nmRms):
    """
    <N>            : size of the output image, that is made to match the size
                     of the (square) psf image to be processed. In other
                     words, N = psf.shape[0]
    <pixelSize>    : size of the pixels of the psf image (mas)
    <wavelengthIR> : IR wavelength (metres)
    <rotdegree>    : pupil rotation (degrees)
    <r0Vis>        : value of r0 in the visible (metres)
    <nmRms>        : number of nm rms affecting the psf (nm)

    Do not use that function.
    It's just there to create a quicky random shitty psf.
    It's there because i needed a random shitty psf.
    So i wrote it.
    And it's still there.

    Example::

        N = 512
        pixelSise = 4.2         # mas
        wavelengthIR = 1.65e-6  # metres
        rotdegree = 10.
        r0Vis = 0.12
        nmRms = 150.
        psf, pup = createAdHocScaoPsf(N, pixelSize, wavelengthIR, rotdegree,
                                      r0Vis, nmRms)

    """

    # let's compute r0 in the IR using the
    # r0 chromatic translation formula
    wavelengthVis = 500e-9
    r0IR = r0Vis * (wavelengthIR / wavelengthVis)**(6/5.)

    kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelengthIR)
    M4 = defineDmFrequencyArea(kx, ky, rotdegree)

    # I hardcode an outer scale of 25 metres
    L0 = 25.

    # This is the turbulent spectrum ....
    W = computeWiener(kx, ky, L0, r0IR)

    # And here are some of the PSF-destroyers
    Wfit = fittingSpectrum(W, M4)
    nmRms = 200.
    Wother = otherSpectrum(nmRms, M4, uk, wavelengthIR)
    Dphi = convertSpectrum2Dphi(Wfit + Wother, uk)

    # Here, i need to generate a kind of PSF or telescope OTF
    deadSegments = 3
    pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize,
                             wavelengthIR)
    FTOtel = computeEeltOTF(pup)

    psf = core_generatePsf(Dphi, FTOtel)

    return psf, pup


def r0Converter(r0, lambda1, lambda2):
    """
    Converts a r0 defined at some wavelength lambda1,
    into a r0' at lambda2.
    Lambda1 and 2 shall be in the SAME unit.
    Returned r0 will be in the SAME unit than the input one.

    Example::

        r0_K = r0Converter(0.12, 500, 2200)

    """
    return r0 * (lambda2/lambda1)**(6/5.)


def airmassImpact(r0_at_zenith, zenith_distance):
    """
    <r0_at_zenith>    : r0 at zenith (any unit is ok)
    <zenith_distance> : zenith distance (degrees)

    The seeing/r0 actually perceived by a telescope depends on the zenith
    distance. This function converts a r0 given at zenith into the real r0
    actually observed by the telescope.
    """

    z = zenith_distance * np.pi / 180  # the same, in radians
    r0 = r0_at_zenith * np.cos(z)**(3./5)

    return r0


def get_atmospheric_turbulence(myProfile='officialEsoMedian'):
    """
    Returns the relative level of turbulence at a given height

    Note: The np.sum(Cn2h) = 1.0

    The following 3 profile are currently available:

    * ``oldEso``

        The old ESO profile from before 2010. Cn2h paramateres are taken from
        ref. E-SPE-ESO-276-0206_atmosphericparameters

    * ``officialEsoMedian``

        Median Armazones atmospheric profile directly copied from doc.
        ESO-258292 "Relevant Atmospheric Parameters for E-ELT AO Analysis and
        Simulations".

    * ``gendron``

        The Gendron profile. Short, fast. Saves CPU. Carbon efficient.
        The np.sum(Cn2h) = 1.0 is hyper-guaranteed here.


    Parameters
    ----------
    myProfile : str, optional
        Profile name: ['oldEso', 'officialEsoMedian', 'gendron']

    Returns
    -------
    layerAltitude : list of floats
        [m] height of layer above ground
    Cn2h : list of floats
        Relative strength of turbulence

    """

    layerAltitude, Cn2h = [], []

    if myProfile == 'oldEso':
        layerAltitude = [47., 140, 281, 562, 1125, 2250, 4500, 9000, 18000.]
        Cn2h = [0.5224, 0.026, 0.0444, 0.116, 0.0989,
                0.0295, 0.0598, 0.043, 0.06]

    elif myProfile == 'officialEsoMedian':
        layerAltitude = [30, 90, 150, 200, 245, 300, 390, 600, 1130, 1880, 2630,
                         3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500, 11500,
                         12500, 13500, 14500, 15500, 16500, 17500, 18500, 19500,
                         20500, 21500, 22500, 23500, 24500, 25500, 26500]
        Cn2h = [24.2, 12, 9.68, 5.9, 4.73, 4.73, 4.73, 4.73, 3.99, 3.24, 1.62,
                2.6, 1.56, 1.04, 1, 1.2, 0.4, 1.4, 1.3, 0.7, 1.6, 2.59, 1.9,
                0.99, 0.62, 0.4, 0.25, 0.22, 0.19, 0.14, 0.11, 0.06, 0.09, 0.05,
                0.04]
        Cn2h = np.array(Cn2h)
        Cn2h /= np.sum(Cn2h)

    elif myProfile == 'gendron':
        Cn2h = [1.0]
        layerAltitude = [4414.]

    return layerAltitude, Cn2h


def clean_psf(psf, threshold):
    psf[psf < threshold] = 0
    edge_threshold = np.median([psf[:, 0], psf[0, :]])
    psf[psf < edge_threshold] = 0
    psf /= np.sum(psf)

    return psf


def round_edges(kernel, edge_width=10):
    n = edge_width
    falloff = np.cos(1.5708 * np.arange(n) / (n-1)).reshape([1, n])

    kernel[:n, :] *= falloff.T[::-1, :]
    kernel[-n:, :] *= falloff.T
    kernel[:, :n] *= falloff[:, ::-1]
    kernel[:, -n:] *= falloff

    return kernel




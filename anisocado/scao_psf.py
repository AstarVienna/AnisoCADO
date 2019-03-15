import warnings
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits

from .utils import *


class AnalyticalScaoPsf:
    """
    A class to generate SCAO PSFs for MICADO at the ELT

    Parameters
    ----------
    N : int
        [pixel] Side-length of the kernel array
    pixelSize : float
        [mas] On-sky pixel scale
    wavelengthIR : float
        [m] Wavelength for which the PSF should be generated
    rotdegree : float
        [deg] Rotation angle of the pupil w.r.t the plane of the optical axis
    seeing : float
        [arcsec]
    nmRms : float
        [nm] Residual wavefront error of the system
    L0 : float
        [m] Outer scale
    profile_name : str
        ["oldEso", "officialEsoMedian", "gendron"]. Names of specific turbulence
        profiles for which presets exist
    zenDist : float
        [degree] Zenith distance
    deadSegments : int
        Number of segments missing from the ELT primary mirror
    V : float
        [m/s] Wind speed
    Fe : float
        [Hz] AO sampling frequency of the system. Default is 500 Hz
    tret : float
        [s] Delay in the AO loop. Default is 4 ms
    gain : float
        Closed-loop gain default is 0.3
    dactu : float
        [m] Interactuator distance on M4. Default in 0.5403
    _last_x, _last_y : float
        [arcsec] shifts used to generate the psf_latest kernel


    Derived Attributes
    ------------------
    r0Vis : float
        [m] Fried parameter at 500 nm.
    r0IR : float
        [m] Fried parameter at 500 nm.
    _on_axis_psf : array
        The on-axis PSF kernel
    _last_psf : array
        The last PSF generated or shifted is kept in memory
    pup : array
        An image of the telescope pupil, i.e. ELT M1
    layerAltitude: list, array
        [m] The heights of the different turbulent layers
    Cn2h : array, list
        [0..1] The relative strength of each turbulent layer
    kx, ky : 2D meshgrid arrays
        [m-1] spatial frequencies - see ``utils.computeSpatialFreqArrays``
    uk : array
        [mas / m]
    M4 : array
        2D domain of spatial frequencies which can be impacted by M4, i.e. the
        M4 compensation domain. See ``utils.defineDmFrequencyArea``
    W :
        [rad2 m2] 2D spectrum of the turbulence with an outer scale L0.
        See ``utils.computeWiener``

    """

    def __init__(self, **psf_on_axis):
        self.kwarg_names = ["N", "pixelSize", "wavelengthIR", "rotdegree",
                            "seeing", "r0Vis", "r0IR", "nmRms", "L0",
                            "profile_name", "zenDist", "deadSegments",
                            "V", "Fe", "tret", "gain", "dactu", "x_last",
                            "y_last"]

        # Variable attributes
        self.N = 512
        self.pixelSize = 4          # mas
        self.wavelengthIR = 2.15e-6 # metres
        self.rotdegree = 0.         # deg
        self.seeing = 0.8           # arcsec
        self.nmRms = 100.           # nm
        self.L0 = 25.               # meters - from doc. ESO-258292.
        self.profile_name = "gendron"
        self.zenDist = 30.          # degree
        self.deadSegments = 5       # there are some missing segments tonight !
        self.V = 10.                # wind is 10 m/s
        self.Fe = 500.              # sampling frequency of the system is 500 Hz
        self.tret = 0.004           # delay in the loop is 4 ms
        self.gain = 0.3             # closed-loop gain is 0.3
        self.dactu = 0.5403         # [m] distance betweem actuators on M4
        self.x_last = 0            # [arcsec] x shift used to make psf_latest
        self.y_last = 0            # [arcsec] x shift used to make psf_latest

        # Derived attributes
        self.r0Vis = None           # meters
        self.r0IR = None            # meters
        self.psf_on_axis = None
        self.psf_latest = None
        self.pup = None
        self.layerAltitude = None
        self.Cn2h = None
        self.kx = None
        self.ky = None
        self.uk = None
        self.M4 = None
        self.W = None

        self.__dict__.update(psf_on_axis)
        self.update()

    def update(self, **kwargs):
        """
        Updates the parameter needed to generate a PSF and/or shift if off-axis

        Valid parameter names can be found in ``self.kwarg_names``

        """
        for key in kwargs:
            if key not in self.kwarg_names:
                warnings.warn("{} not found in self.kwarg_names".format(key))

        self.__dict__.update(kwargs)

        # and layers appear further away with zenith distance
        if self.profile_name is not None:
            layerAltitude, Cn2h = get_atmospheric_turbulence(self.profile_name)
            self.layerAltitude = np.array(layerAltitude)
            self.Cn2h = Cn2h
        self.layerAltitude *= 1 / np.cos(self.zenDist * np.pi / 180.)

        # r0Vis is in metres here (0.103 is in metres.arcsec)
        # convert r0 from 500nm to IR
        # apparent seeing degrades with airmass
        self.r0Vis = 0.103 / self.seeing
        self.r0IR = r0Converter(self.r0Vis, 500e-9, self.wavelengthIR)
        self.r0IR = airmassImpact(self.r0IR, self.zenDist)

        # generate the ELT pupil
        self.pup = fake_generatePupil(self.N, self.deadSegments, self.rotdegree,
                                      self.pixelSize, self.wavelengthIR)

        # Let's go. Let's define some basic parameters
        # (arrays of spatial frequencies)
        self.kx, self.ky, self.uk = computeSpatialFreqArrays(self.N,
                                                             self.pixelSize,
                                                             self.wavelengthIR)
        self.M4 = defineDmFrequencyArea(self.kx, self.ky, self.rotdegree,
                                        self.dactu)

        # This is the turbulent spectrum ....
        self.W = computeWiener(self.kx, self.ky, self.L0, self.r0IR)

        if self.psf_on_axis is None:
            self.psf_on_axis = self.make_psf()
            self.psf_latest = self.psf_on_axis

    def make_psf(self):
        """
        Generates a analytical SCAO PSF for a long (>10 sec) exposure

        Parameters need to be set be setting the attibute directly, or by
        calling ``self.update()`` with the desired keyword-value pair passed as
         a kwarg. Valid keywords can be found in ``self.kwarg_names``.

        Returns
        -------
        psf : array
            The PSF kernel

        """
        # Initial setup happens in update()

        dx, dy = 0, 0
        # And here are some of the PSF-destroyers - (English: Wavefront errors)
        Waniso = anisoplanaticSpectrum(self.Cn2h, self.layerAltitude, self.L0,
                                       dx, dy, self.wavelengthIR,
                                       self.kx, self.ky, self.W, self.M4)
        Wfit = fittingSpectrum(self.W, self.M4)
        Walias = aliasingSpectrum(self.kx, self.ky, self.r0IR, self.L0, self.M4)
        Wbp = computeBpSpectrum(self.kx, self.ky, self.V, self.Fe, self.tret,
                                self.gain, self.W, self.M4)
        Wother = otherSpectrum(self.nmRms, self.M4, self.uk, self.wavelengthIR)
        # Wnoise = noiseSpectrum(Rmag, .. kx, ky)  available one day ...

        # Now, you sum up every contributor, and produce a phase structure funcn
        Wtotal = Waniso + Wfit + Wother + Walias + Wbp  # + Wnoise
        Dphi = convertSpectrum2Dphi(Wtotal, self.uk)

        # And you "blur" the nice Airy pattern using that phase structure funcn
        FTOtel = computeEeltOTF(self.pup)
        psf = core_generatePsf(Dphi, FTOtel)
        psf /= np.sum(psf)
        # psf = clean_psf(psf, 1E-7)
        self.psf_latest = psf

        return psf

    def shift_psf_off_axis(self, dx, dy):
        """
        Shifts the on-axis PSF off axis by an amount ``(dx, dy)``

        Parameters
        ----------
        dx, dy : float
            [arcsec] Offset in each of the dimensions relative to the plane of
            the optical axis

        Returns
        -------
        psf : array
            The PSF kernel

        """
        self.x_last = dx
        self.y_last = dy

        # Original setup starts in update()

        # and all this will be used to run the function below, that will
        # compute the spatial spectrum of the phase due to anisoplanatism
        Waniso = anisoplanaticSpectrum(self.Cn2h, self.layerAltitude, self.L0,
                                       dx, dy, self.wavelengthIR,
                                       self.kx, self.ky, self.W, self.M4)

        # Transforming this spectrum into a phase structure function
        Dphi = convertSpectrum2Dphi(Waniso, self.uk)

        # Here, the on-axis psf comes into play ... I take its Fourier transform
        # it's complex.
        fto = np.fft.fft2(np.fft.fftshift(self.psf_on_axis)) / self.N ** 2

        psf = core_generatePsf(Dphi, fto)
        psf /= np.sum(psf)
        self.psf_latest = psf

        return psf

    def make_short_exposure_psf(self, DIT=1.0, step=0.5):
        """
        Returns a PSF for an 'short' exposure time of ``DIT``

        The PSF kernel will be a single 2D array made from N stacked
        instantaneous PSFs, where the instantaneous PSFs are generated at time
        intervals during the DIT length determined by the wind speed,
        ``self.V``, and the phase-screen ``step`` length.

        Parameters
        ----------
        DIT : float
            [s] Default is 1.0 sec. Exposure time for the PSF

        step : float
            [m] Sample step length for atmospheric phase screen
            Default is 0.5m - the length of the M4 actuator pitch

        Returns
        -------
        psfLE : array

        """

        # The dirty one.
        # Let's try to simulate the fluctuations due to short exposures.
        Waniso = anisoplanaticSpectrum(self.Cn2h, self.layerAltitude, self.L0,
                                       self.x_last, self.x_last,
                                       self.wavelengthIR, self.kx, self.ky,
                                       self.W, self.M4)
        Wfit = fittingSpectrum(self.W, self.M4)
        Walias = aliasingSpectrum(self.kx, self.ky, self.r0IR, self.L0, self.M4)
        Wbp = computeBpSpectrum(self.kx, self.ky, self.V, self.Fe, self.tret,
                                self.gain, self.W, self.M4)
        Wother = otherSpectrum(self.nmRms, self.M4, self.uk, self.wavelengthIR)
        # Wnoise = noiseSpectrum(Rmag, .. kx, ky)  available one day ...

        # I start from "use case 2", and I sum the contributors to the phase err
        # What I get is the total power spectrum of the perturbed phase.
        Wtotal = Waniso + Wfit + Wother + Walias + Wbp

        # So, i'm gonna do some random draw of a phase that follows the
        # statistics of the spectrum WW. For that, i'm gonna use sqrt(WW) as
        # the modulus of the FFT of the phase, and generate a random phase
        # chosen uniformly between 0 and 2.pi. I will then do a FFT of that, in
        # order to get the phase.
        Wtotal[0, 0] = 0  # because i don't care about piston mode
        Wtotal = np.sqrt(Wtotal)
        rand_phase = np.random.rand(self.N, self.N)
        tmp = np.fft.fft2(Wtotal * np.exp(2j * np.pi * rand_phase)) * self.uk
        ph1 = tmp.real * np.sqrt(2)
        ph2 = tmp.imag * np.sqrt(2)

        # now i compute widthScreen, the size of the pixels of the phase screens
        # I have generated.
        widthScreen = 1. / self.uk  # in metres
        ud = widthScreen / self.N  # size of the pixels of the phase screen

        # With such a wide screen, and using a wind speed of V m/s, then I can
        # simulate an exposure time of  (widthScreen/V) seconds.
        stepPix = int(np.round(step / ud))
        stepTime = (stepPix * ud) / self.V
        niter = int(np.round(DIT / stepTime))
        psfLE = 0  # psf Long Exposure
        normFactor = np.sum(self.pup) ** 2
        for i in range(niter):
            inv_pupil = np.fft.fft2(self.pup * np.exp(1j * ph1))
            psfSE = np.fft.fftshift(np.abs(inv_pupil) ** 2)
            psfSE /= normFactor
            # print(psfSE.max())
            psfLE += psfSE
            ph1 = np.roll(ph1, stepPix, axis=0)
        psfLE /= niter
        psfLE /= np.sum(psfLE)
        self.psf_latest = psfLE

        return psfLE

    @property
    def strehl_ratio(self):
        return np.max(self.psf_latest)

    @property
    def kernel(self):
        return self.psf_latest

    @property
    def hdu(self):
        return self.get_hdu()

    def get_hdu(self):
        w, h = self.psf_latest.shape

        hdr = fits.Header()
        hdr["CDELT1"] = self.pixelSize / (3600. * 1000.)
        hdr["CDELT2"] = self.pixelSize / (3600. * 1000.)   #because pixelSize is in mas
        hdr["CRVAL1"] = self.x_last / 3600.
        hdr["CRVAL2"] = self.y_last / 3600.
        hdr["CRPIX1"] = w / 2.
        hdr["CRPIX2"] = h / 2.
        hdr["CTYPE1"] = "RA---TAN"
        hdr["CTYPE2"] = "DEC--TAN"
        hdr["CUNIT1"] = "degree"
        hdr["CUNIT2"] = "degree"
        hdr["WAVE0"] = (self.wavelengthIR * 1e6, "[um] Central wavelength")

        dic = {key: self.__dict__[key] for key in self.kwarg_names}
        hdr.update(dic)
        hdu_psf = fits.ImageHDU(data=self.psf_latest, header=hdr)

        return hdu_psf

    def plot_psf(self, which="psf_latest"):
        plt.imshow(getattr(self, which).T, origin='l', norm=LogNorm())
        print('Strehl ratio of {} is {}'.format(which, self.psf_latest.max()))


AnalyticalScaoPsf().hdu.writeto()
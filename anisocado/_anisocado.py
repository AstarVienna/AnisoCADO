#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:48:47 2019

@author: Eric Gendron
Edited by Kieran Leschinski
"""

import matplotlib.pyplot as plt
from anisocado.psf_utils import *
from anisocado.psf_utils import get_atmospheric_turbulence

#  _   _             ____                  _
# | | | |___  ___   / ___|__ _ ___  ___   / |
# | | | / __|/ _ \ | |   / _` / __|/ _ \  | |
# | |_| \__ \  __/ | |__| (_| \__ \  __/  | |
#  \___/|___/\___|  \____\__,_|___/\___|  |_|


def shift_scao_psf(plots=False):
    """
    Problem:
    You have an on-axis PSF.
    You want to 'move' it off-axis, let's say (+15, +20) arcsec.


    For that, you will need to know:
        - the Cn2h profile
        - the global r0
        - the global L0

    """
    ###########
    # Setup

    # Let's take an example. I create a PSF, that could be coming from
    # a SCAO simulation. It has 512x512 pixels, sampled with 4.2 mas in H band.
    # Pupil is rotated by 10 deg. The r0 was 12cm.
    #
    # That psf will be the starting point.
    N = 512
    pixelSize = 4.2         # mas
    wavelengthIR = 1.65e-6  # metres
    rotdegree = 10.         # deg
    r0Vis = 0.12
    nmRms = 150.
    psf, pup = createAdHocScaoPsf(N, pixelSize, wavelengthIR, rotdegree, r0Vis,
                                  nmRms)

    if plots:
        # I can even look at it:
        plt.imshow(psf.T, origin='l')
        print('Strehl ratio of initial psf is ', psf.max())

    # OK. That's the starting point.......................

    # Now I need to know the atmospheric properties, in particular the Cn2h
    # profile. Let me offer you a little selection of atmospheric profiles.
    profile_name = "gendron"
    layerAltitude, Cn2h = get_atmospheric_turbulence(profile_name)

    # Let us define the outer scale value. 25 metres is the Armazones median
    # value from doc. ESO-258292.
    # Choose 15 m for a lucky observer. 50 m for the looser.
    L0 = 25.   # 25 metres is the Armazones median value from doc. ESO-258292.

    # Now, the seeing.
    # Here, we're cheating, we already know the r0 is 12cm because we've
    # generated the PSF with this.
    r0Vis = 0.12  # I know, we know it already ...

    # I also need to know where's the off-axis star I want to simulate
    # in arcsecs (this one is for those who'd like to check nothing will change
    # at the end)
    # offx, offy = (0, 0)
    offx, offy = (0., 16.)   # in arcsecs

    ####################
    # Generate PSF

    # Then let's start the work. I will create spatial frequency arrays.
    kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelengthIR)

    # convert r0 in the infra-red
    r0IR = r0Converter(r0Vis, 500e-9, wavelengthIR)

    # and create M4 working domain
    M4 = defineDmFrequencyArea(kx, ky, rotdegree)

    # and finally the turbulent spectrum ....
    W = computeWiener(kx, ky, L0, r0IR)

    # and all this will be used to run that function below, that will compute
    #  the spatial spectrum of the phase due to anisoplanatism
    Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy,
                                   wavelengthIR, kx, ky, W, M4)

    # Transforming this spectrum into a phase structure function
    Dphi = convertSpectrum2Dphi(Waniso, uk)

    # Here, the on-axis psf comes into play ... I take its Fourier transform
    fto = np.fft.fft2(np.fft.fftshift(psf)) / N**2  # it's complex.

    psf_aniso = core_generatePsf(Dphi, fto)

    print('Strehl off-axis is', psf_aniso.max())
    plt.imshow(psf_aniso.T, origin='l' )

    return psf_aniso

#  _   _             ____                 ____
# | | | |___  ___   / ___|__ _ ___  ___  |___ \
# | | | / __|/ _ \ | |   / _` / __|/ _ \   __) |
# | |_| \__ \  __/ | |__| (_| \__ \  __/  / __/
#  \___/|___/\___|  \____\__,_|___/\___| |_____|


def exnihilo_scao_psf():
    """

    You want to generate SCAO PSFs ex-nihilo, using a simple, approximative
    simulation software.
    You are aware that the simulated PSFs will be perfectly smooth (infinitely
    converged), they do not reflect the fluctuations associated to short
    exposures.

    For this you need to know all the parameters of the simulation.

    """

    ###########
    # Setup

    # nb of pixels of the image to be simulated.
    N = 1024
    pixelSize = 3.1  # mas

    # this is where you observe
    wavelengthIR = 2.2e-6   # oh! a K band ...

    # Now you need to tell what turbulence looks like.
    # Here you can imagine to take those numbers into ESO doc, statistics,
    # you can also include some dependence of r0 with respect to airmass
    #
    # If everything is as usual, you should end up with too many parameters that
    # nobody knows about, and someone will tell you "ok it's very nice, but can
    # you please simplify this ?"
    layerAltitude = [47., 140, 281, 562, 1125, 2250, 4500, 9000, 18000.]
    # from ref. E-SPE-ESO-276-0206_atmosphericparameters
    Cn2h = [0.5224, 0.026, 0.0444, 0.116, 0.0989, 0.0295, 0.0598, 0.043, 0.06]
    L0 = 25.   # 25 metres is the Armazones median value
    seeing = 0.8          # in arcseconds
    r0Vis = 0.103 / seeing  # r0Vis is in metres here, 0.103 is in metres.arcsec
    r0IR = r0Converter(r0Vis, 500e-9, wavelengthIR) # convert r0 at 500nm to IR

    # Just to use that wonderful function, i decide that the seeing given here
    # was expressed at zenith, while our telescope observes at 30° from zenith.
    # This will transform our r0 into the effective/actual one seen at 30°.
    # In addition, the turbulent layers will all appear further away from the
    # telescope, so that their apparent distance grows with airmass --> this is
    # very bad for anisoplanatism !..
    zenDist = 30.  # I observe at 30° from zenith
    r0IR = airmassImpact(r0IR, zenDist)  # apparent seeing degrades with airmass

    layerAltitude = np.array(layerAltitude)
    layerAltitude *= 1/np.cos(zenDist*np.pi/180)  # layers appear further away

    # Also you may want to say something about how the pupil is rotated wrt
    # your image
    rotdegree = 10.0

    # And you also need to generate the EELT pupil properly
    deadSegments = 5  # there are some missing segments tonight !
    pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize,
                             wavelengthIR)

    # For temporal aspects you need to know the characteristics of your system
    V = 10.        # wind is 10 m/s
    Fe = 500.      # sampling frequency of the system is 500 Hz
    tret = 0.004   # delay in the loop is 4 ms
    gain = 0.3     # closed-loop gain is 0.3

    # Here is the position of the object in the field
    offx, offy = (10., 0.)


    #################
    # PSF generation

    # Let's go. Let's define some basic parameters (arrays of spatial
    # frequencies)
    kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelengthIR)
    M4 = defineDmFrequencyArea(kx, ky, rotdegree)

    # This is the turbulent spectrum ....
    W = computeWiener(kx, ky, L0, r0IR)

    # And here are some of the PSF-destroyers - (English: Wavefront errors)
    Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy,
                                   wavelengthIR, kx, ky, W, M4)
    Wfit = fittingSpectrum(W, M4)
    Walias = aliasingSpectrum(kx, ky, r0IR, L0, M4)
    Wbp = computeBpSpectrum(kx, ky, V, Fe, tret, gain, W, M4)
    nmRms = 100.
    Wother = otherSpectrum(nmRms, M4, uk, wavelengthIR)

    # THE missing term = noise
    # Wnoise = noiseSpectrum(Rmag, .. kx, ky)  available one day ...

    # Now, you sum up every contributor, and produce a phase structure function
    Dphi = convertSpectrum2Dphi(Waniso + Wfit + Wother + Walias + Wbp, uk)

    # And you "blur" the nice Airy pattern using that phase structure function
    FTOtel = computeEeltOTF(pup)
    psf = core_generatePsf(Dphi, FTOtel)

    print('Strehl is ', psf.max())
    plt.imshow( np.log(psf) )

    return psf


#  _   _             ____                 _____
# | | | |___  ___   / ___|__ _ ___  ___  |___ /
# | | | / __|/ _ \ | |   / _` / __|/ _ \   |_ \
# | |_| \__ \  __/ | |__| (_| \__ \  __/  ___) |
#  \___/|___/\___|  \____\__,_|___/\___| |____/


def instantaneous_scao_psf():
    """

    The dirty one.
    Let's try to simulate the fluctuations due to short exposures.

    """
    # I start from "use case 2", and I sum the contributors to the phase error.
    # What I get is the total power spectrum of the perturbed phase.
    WW = Waniso + Wfit + Wother + Walias + Wbp

    # So, i'm gonna do some random draw of a phase that follows the statistics
    # of the spectrum WW. For that, i'm gonna use sqrt(WW) as the modulus of the
    # FFT of the phase, and generate a random phase chosen uniformly between 0
    # and 2.pi. I will then do a FFT of that, in order to get the phase.
    WW[0, 0] = 0        # because i don't care about piston mode
    WW = np.sqrt(WW)
    tmp = np.fft.fft2(WW * np.exp(2j * np.pi*np.random.rand(N, N))) * (uk)
    ph1 = tmp.real * np.sqrt(2)
    ph2 = tmp.imag * np.sqrt(2)

    # now i compute widthScreen, the size of the pixels of the phase screens I
    # have generated.
    widthScreen = 1. / uk  # in metres
    ud = widthScreen / N   # size of the pixels of the phase screen

    # With such a wide screen, and using a wind speed of V m/s, then I can
    # simulate an exposure time of  (widthScreen/V) seconds.
    # I recommend to sum psf snapshots every 50 cm (actuator pitch of M4).
    step = 0.50
    stepPix = int(np.round(step / ud))
    stepTime = (stepPix * ud) / V
    DIT = 1.0  # 1 second integration time
    niter = int(np.round(DIT / stepTime))
    psfLE = 0  # psf Long Exposure
    normFactor = np.sum(pup)**2
    for i in range(niter):
        psfSE = np.fft.fftshift(np.abs(np.fft.fft2(pup * np.exp(1j*ph1)))**2)
        psfSE /= normFactor
        print(psfSE.max())
        psfLE += psfSE
        ph1 = np.roll(ph1, stepPix, axis=0)
    psfLE /= niter

    # Here, possibilities are infinite ..
    # You can add some static aberrations, etc etc,
    # and generate all the PSFs you want.








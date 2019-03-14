#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:48:47 2019

@author: gendron
"""


import numpy as np
import matplotlib.pyplot as plt
plt.ion()



"""


 ____  _____    _    ____  __  __ _____
|  _ \| ____|  / \  |  _ \|  \/  | ____|
| |_) |  _|   / _ \ | | | | |\/| |  _|
|  _ <| |___ / ___ \| |_| | |  | | |___
|_| \_\_____/_/   \_\____/|_|  |_|_____|



Hello.
This files contains some useful functions related to psf generation.
They are just raw, so that you can insert them into your own classes
as desired.
The functions are commented to help you understand what's going on, when
you need to modify them.


Don't read that code anyway.
Go right now to line 600 of the file.
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
    
    Standalone example:
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
    A = np.pi/3 # 60 degrees
    A0 = rotdegree * np.pi / 180
    msk = np.abs(np.cos(A0)*ky+np.sin(A0)*kx)<fc
    msk = np.logical_and(msk, np.abs(np.cos(A+A0)*ky+np.sin(A+A0)*kx)<fc)
    msk = np.logical_and(msk, np.abs(np.cos(2*A+A0)*ky+np.sin(2*A+A0)*kx)<fc)
    k = np.sqrt(kx**2 + ky**2)
    msk = np.logical_and(msk, k<(fc*1.115))
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
    
    Standalone example:
    N = 512               # output will be 512x512
    pixelSize = 4.2       # 4.2 mas
    wavelength = 1.65e-6  # metres
    kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelength)
    """
    # array of indices centred 'as expected' by Fourier frequencies, in 1D
    k1d = np.fft.fftshift(np.arange(N) - (N//2))
    # Proper scaling to transform indices in spatial frequencies.
    dX = pixelSize * 4.84813681109536e-09  # from mas to radians
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
    
    Standalone example:
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
    Wiener[0, 0] = 0.00
    # done
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
    f = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy, wavelength, kx, ky, W, M4)

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
    Usage:
    
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
    nu = k * V / np.sqrt(2)    # pourquoi un sqrt(2) ? je ne saurais dire ...!!!!
    Wbp = hcor(nu, Fe, tret, gain, 500) * Wiener
    Wbp[np.logical_not(M4)] = 0.00
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
    
    Sandalone usage:
    N = 512
    deadSegments = 3
    rotdegree = 14.
    pixelSize = 4.2
    wavelength = 1.65e-6
    pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize, wavelength)
    """
    nseg = getEeltSegmentNumber()
    refl = np.ones(nseg)+np.random.randn(nseg)/20.
    if deadSegments:
        refl[(np.random.rand(deadSegments)*nseg).astype(int)] = 0.
    i0 = N/2+0.5
    j0 = N/2+0.5
    # field of view of the psf image in rd
    FoV = N * pixelSize * 4.84813681109536e-09  # from mas to radians
    # pixel scale of pupil image
    pixscale = wavelength / FoV   # expressed in metres
    dspider = 0.53
    gap = 0.02
    pup = generateEeltPupilReflectivity(refl, N, dspider, i0, j0, pixscale,
                                        gap, rotdegree, softGap=True)
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
    
    Standalone example:
    
    N = 1024
    pixelSize = 4.2
    # atmospheric profile (old ESO profile before 2010)
    layerAltitude = [47., 140, 281, 562, 1125, 2250, 4500, 9000, 18000.]
    Cn2h = [52.24, 2.6, 4.44, 11.60, 9.89, 2.95, 5.98, 4.30, 6] # from ref. E-SPE-ESO-276-0206_atmosphericparameters
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
    Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy, wavelength, kx, ky, W, M4)
    Wfit = fittingSpectrum(W, M4)
    nmRms = 100.
    Wother = otherSpectrum(nmRms, M4, uk, wavelength)
    Dphi = convertSpectrum2Dphi(Waniso + Wfit + Wother, uk)

    # Here, i need to generate a kind of PSF or telescope OTF
    deadSegments = 3
    pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize, wavelength)
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
    
    N = 512
    pixelSise = 4.2         # mas
    wavelengthIR = 1.65e-6  # metres
    rotdegree = 10.
    r0Vis = 0.12
    nmRms = 150.
    psf, pup = createAdHocScaoPsf(N, pixelSize, wavelengthIR, rotdegree, r0Vis, nmRms)
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
    pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize, wavelengthIR)
    FTOtel = computeEeltOTF(pup)
    
    psf = core_generatePsf(Dphi, FTOtel)
    return psf, pup







def r0Converter(r0, lambda1, lambda2):
    """
    Converts a r0 defined at some wavelength lambda1,
    into a r0' at lambda2.
    Lambda1 and 2 shall be in the SAME unit.
    Returned r0 will be in the SAME unit than the input one.
    
    Example:
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
















"""

 _   _             ____                  _ 
| | | |___  ___   / ___|__ _ ___  ___   / |
| | | / __|/ _ \ | |   / _` / __|/ _ \  | |
| |_| \__ \  __/ | |__| (_| \__ \  __/  | |
 \___/|___/\___|  \____\__,_|___/\___|  |_|
                                           



Problem:
You have an on-axis PSF.
You want to 'move' it off-axis, let's say (+15, +20) arcsec.


For that, you will need to know:
    - the Cn2h profile
    - the global r0
    - the global L0



"""

if False:
    # Let's take an example. I create a PSF, that could be coming from
    # a SCAO simulation. It has 512x512 pixels, sampled with 4.2 mas in H band.
    # Pupil is rotated by 10°. The r0 was 12cm.
    #
    # That psf will be the starting point.
    N = 512
    pixelSize = 4.2         # mas
    wavelengthIR = 1.65e-6  # metres
    rotdegree = 10.         # deg
    r0Vis = 0.12
    nmRms = 150.
    psf, pup = createAdHocScaoPsf(N, pixelSize, wavelengthIR, rotdegree, r0Vis, nmRms)
    
    #I can even look at it:
    plt.imshow(psf.T, origin='l')
    print('Strehl ratio of initial psf is ', psf.max())
    
    # OK. That's the starting point.......................
    
    
    # Now I need to know the atmospheric properties, in particular the Cn2h profile.
    # Let me offer you a little selection of atmospheric profiles.
    myProfile = 'officialEsoMedian'
    if myProfile=='oldEso':
        # atmospheric profile (old ESO profile before 2010)
        # The np.sum(Cn2h) = 1.0
        Cn2h = [0.5224, 0.026, 0.0444, 0.116, 0.0989, 0.0295, 0.0598, 0.043, 0.06] # from ref. E-SPE-ESO-276-0206_atmosphericparameters
        layerAltitude = [47., 140, 281, 562, 1125, 2250, 4500, 9000, 18000.]
    elif myProfile=='officialEsoMedian':
        # median Armazones atmospheric profile
        # directly copied from doc. ESO-258292 "Relevant Atmospheric Parameters
        # for E-ELT AO Analysis and Simulations".
        # The np.sum(Cn2h) = 1.0
        Cn2h = [24.2,12,9.68,5.9,4.73,4.73,4.73,4.73,3.99,3.24,1.62,2.6,1.56,1.04,1,1.2,0.4,1.4,1.3,0.7,1.6,2.59,1.9,0.99,0.62,0.4,0.25,0.22,0.19,0.14,0.11,0.06,0.09,0.05,0.04]
        layerAltitude = [30,90,150,200,245,300,390,600,1130,1880,2630,3500,4500,5500,6500,7500,8500,9500,10500,11500,12500,13500,14500,15500,16500,17500,18500,19500,20500,21500,22500,23500,24500,25500,26500]
        Cn2h = np.array(Cn2h)
        Cn2h /= np.sum(Cn2h)
    elif myProfile=='gendron':
        # The Gendron profile. Short, fast. Saves CPU. Carbon efficient.
        # The np.sum(Cn2h) = 1.0 is hyper-garanteed here.
        Cn2h = [1.0]
        layerAltitude = [4414.]
    
    # Let us define the outer scale value. 25 metres is the Armazones median
    # value from doc. ESO-258292.
    # Choose 15 m for a lucky observer. 50 m for the looser.
    L0 = 25.   # 25 metres is the Armazones median value from doc. ESO-258292.
    
    # Now, the seeing.
    # Here, we're cheating, we already know the r0 is 12cm because we've 
    # generated the PSF with this.
    r0Vis = 0.12  # I know, we know it already ...
    
    # I also need to know where's the off-axis star I want to simulate
    offx, offy = (0, 0)   # in arcsecs (this one is for those who'd like to check nothing will change at the end)
    offx, offy = (0., 16.)   # in arcsecs
    
    
    # Then let's start the work. I will create spatial frequency arrays.
    kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelengthIR)
    # convert r0 in the infra-red
    r0IR = r0Converter(r0Vis, 500e-9, wavelengthIR) # convert r0 from 500nm to IR
    # and create M4 working domain
    M4 = defineDmFrequencyArea(kx, ky, rotdegree)
    # and finally the turbulent spectrum ....
    W = computeWiener(kx, ky, L0, r0IR)
    # and all this will be used to run that function below, that will compute the
    # spatial spectrum of the phase due to anisoplanatism
    Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy, wavelengthIR, kx, ky, W, M4)
    # Transforming this spectrum into a phase structure function
    Dphi = convertSpectrum2Dphi(Waniso, uk)
    # Here, the on-axis psf comes into play ... I take its Fourier transform
    fto = np.fft.fft2(np.fft.fftshift(psf)) / N**2  # it's complex.
    
    psfaniso = core_generatePsf(Dphi, fto)
    
    print('Strehl off-axis is', psfaniso.max())
    plt.imshow( psfaniso.T, origin='l' )
    



"""

 _   _             ____                 ____
| | | |___  ___   / ___|__ _ ___  ___  |___ \
| | | / __|/ _ \ | |   / _` / __|/ _ \   __) |
| |_| \__ \  __/ | |__| (_| \__ \  __/  / __/
 \___/|___/\___|  \____\__,_|___/\___| |_____|

 
You want to generate SCAO PSFs ex-nihilo, using a simple, approximative
simulation software.
You are aware that the simulated PSFs will be perfectly smooth (infinitely
converged), they do not reflect the fluctuations associated to short exposures.

For this you need to know all the parameters of the simulation.

"""

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
Cn2h = [0.5224, 0.026, 0.0444, 0.116, 0.0989, 0.0295, 0.0598, 0.043, 0.06] # from ref. E-SPE-ESO-276-0206_atmosphericparameters
L0 = 25.   # 25 metres is the Armazones median value
seeing = 0.8          # in arcseconds
r0Vis = 0.103 / seeing  # r0Vis is in metres here (0.103 is in metres.arcsec)
r0IR = r0Converter(r0Vis, 500e-9, wavelengthIR) # convert r0 at 500nm to IR

# Just to use that wonderful function, i decide that the seeing given here
# was expressed at zenith, while our telescope observes at 30° from zenith.
# This will transform our r0 into the effective/actual one seen at 30°.
# In addition, the turbulent layers will all appear further away from the
# telescope, so that their apparent distance grows with airmass --> this is
# very bad for anisoplanatism !..
zenDist = 30.  # I observe at 30° from zenith
r0IR = airmassImpact(r0IR, zenDist)   # apparent seeing degrades with airmass

layerAltitude = np.array(layerAltitude)
layerAltitude *= 1/np.cos(zenDist*np.pi/180) # and layers appear further away


# Also you may want to say something about how the pupil is rotated wrt
# your image
rotdegree = 10.0

# And you also need to generate the EELT pupil properly
deadSegments = 45  # there are some missing segments tonight !
pup = fake_generatePupil(N, deadSegments, rotdegree, pixelSize, wavelengthIR)

# For temporal aspects you need to know the characteristics of your system
V = 10.        # wind is 10 m/s
Fe = 500.      # sampling frequency of the system is 500 Hz
tret = 0.004   # delay in the loop is 4 ms
gain = 0.3     # closed-loop gain is 0.3

# Here is the position of the object in the field
offx, offy = (10., 0.)

# Let's go. Let's define some basic parameters (arrays of spatial frequencies)
kx, ky, uk = computeSpatialFreqArrays(N, pixelSize, wavelengthIR)
M4 = defineDmFrequencyArea(kx, ky, rotdegree)
# This is the turbulent spectrum ....
W = computeWiener(kx, ky, L0, r0IR)

# And here are some of the PSF-destroyers
Waniso = anisoplanaticSpectrum(Cn2h, layerAltitude, L0, offx, offy, wavelengthIR, kx, ky, W, M4)
Wfit = fittingSpectrum(W, M4)
Walias = aliasingSpectrum(kx, ky, r0IR, L0, M4)
Wbp = computeBpSpectrum(kx, ky, V, Fe, tret, gain, W, M4)    
nmRms = 100.
Wother = otherSpectrum(nmRms, M4, uk, wavelengthIR)

# THE missing term = noise
#Wnoise = noiseSpectrum(Rmag, .. kx, ky)  available one day ...

# Now, you sum up every contributor, and produce a phase structure function
Dphi = convertSpectrum2Dphi(Waniso + Wfit + Wother + Walias + Wbp, uk)

# And you "blur" the nice Airy pattern using that phase structure function
FTOtel = computeEeltOTF(pup)
psf = core_generatePsf(Dphi, FTOtel)

print('Strehl is ', psf.max())
plt.imshow( np.log(psf) )






"""

 _   _             ____                 _____
| | | |___  ___   / ___|__ _ ___  ___  |___ /
| | | / __|/ _ \ | |   / _` / __|/ _ \   |_ \
| |_| \__ \  __/ | |__| (_| \__ \  __/  ___) |
 \___/|___/\___|  \____\__,_|___/\___| |____/


The dirty one.
Let's try to simulate the fluctuations due to short exposures.

 
"""
# I start from "use case 2", and I sum the contributors to the phase error.
# What I get is the total power spectrum of the perturbed phase.
WW = Waniso + Wfit + Wother + Walias + Wbp


# So, i'm gonna do some random draw of a phase that follows the statistics of
# the spectrum WW. For that, i'm gonna use sqrt(WW) as the modulus of the
# FFT of the phase, and generate a random phase chosen uniformly between 0
# and 2.pi. I will then do a FFT of that, in order to get the phase.
WW[0,0] = 0        # because i don't care about piston mode
WW = np.sqrt(WW)
tmp = np.fft.fft2(WW * np.exp(2j*np.pi*np.random.rand(N,N))) * (uk)
ph1 = tmp.real * np.sqrt(2)
ph2 = tmp.imag * np.sqrt(2)

# now i compute widthScreen, the size of the pixels of the phase screens I have 
# generated.
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
    psfSE = np.fft.fftshift(np.abs(np.fft.fft2( pup * np.exp(1j*ph1) ))**2) 
    psfSE /= normFactor
    print(psfSE.max())
    psfLE += psfSE
    ph1 = np.roll(ph1, stepPix, axis=0)
psfLE /= niter

# Here, possibilities are infinite ..
# You can add some static aberrations, etc etc,
# and generate all the PSFs you want.








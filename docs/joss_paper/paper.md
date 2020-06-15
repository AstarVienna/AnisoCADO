---
title: 'AnisoCADO - a python package for analytically generating adaptive optics point spread functions for the Extremely Large Telescope'

tags:
  - Python
  - astronomy
  - simulations
  - point spread functions
  - Extreme Large Telescope
  
authors:
  - name: Kieran Leschinski
    orcid: 0000-0003-0441-9784
    affiliation: 1
  - name: Eric Gendron
    affiliation: 2
    
affiliations:
 - name: Department of Astrophysics, University of Vienna
   index: 1
 - name: Observatoire de Paris
   index: 2
   
date: 19 May 2020

bibliography: paper.bib

---

# Summary

AnisoCADO is a Python package for generating images of the point spread function (PSF) for the european extremely large telescope. 
The code allows the user to set a large range of the most important atmospheric and observational parameters that influence the shape and strehl ratio of the resulting PSF, including but not limited to: the atmospheric turbulence profile, the guide star position for a single conjugate adaptive optics (SCAO) solution, differential telescope pupil transmission, etc.


# Statement of need

## Adaptive optics are mandatory for the next generation of ground-based telescopes
The larger the telescope aperture, the smaller the diffraction limit of the observations. 
For space-based telescope this statement is always true. 
However the resolution of ground based telescopes is limited by the blur caused by turbulence in the atmosphere - known as atmospheric Seeing. 
This blurring can be (mostly) removed by measuring the deformation of the wavefront of the incoming light, and applying an equal and opposite deformation to the surface of one or more of the mirrors along a telescope's optical path.
The current fleet of large (8-10m) telescopes were built to primarily operate at the edge of the natural seeing limit (FWHM~0.5 arcseconds @ 1um). 
Over the last two decades some have received upgrades in the form of active and adaptive mirrors in order to achieve up to 20x increase in resolution afforded by the physical diffraction limit of a ~10m primary mirror (FWHM~0.03 arcseconds @ 1um).
The next generation of "extremely large" telescopes will have primary mirrors on the order of 30-40m, with theoretical diffraciton limits on the order of 50x smaller than the natural Seeing limit.
In order for these telescopes to resolve structures at scales of the diffraction limit, they must, by design, include adaptive optics systems.

## Diffraction limited point-spread-functions are complex beasts
The point spread function (PSF) of an optical system is the description of the spatial distribution of light from an infinitely small point source after passing through an optical system (e.g. layers of the atmosphere, mirrors of a telescope).
Due to the random nature of atmospheric turbulence, the PSF of a star in a Seeing-limited observation is well approximated by a ("nice") smooth Gaussian-like function.
The PSF of a diffraction limited telescope system using an adaptive-optics correction is a complex ("ugly") function that depends on a veritible zoo of atmospheric, observational, and technical parameters. 
From an astronomers point of view, the consequences of a poor adaptive optics solution means the difference between a successful and a failed observation run.
Therfore it is imperative that the consequences of such large variations in the PSF are accounted for in advance by those proposing to observe with this next generation of billion-dollar telescopes.

## AnisoCADO - Anisoplanatism for MICADO





# Acknowledgments
AnisoCADO uses the following packages as dependencies: 
Numpy [@numpy],
Matplotlib [@numpy],
Astropy, a community-developed core Python package for astronomy [@astropy2018],

This development of this project was funded by the project IS538003 of the Hochschulraum-strukturmittel (HRSM) provided by the Austrian Government and administered by the University of Vienna.


# References




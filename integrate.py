#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit
from srt_functions import *

def calc_beam_area(beam1,beam2):
	return (np.pi * (beam1*np.pi/180.0)*(beam2*np.pi/180.0))/(4.0*np.log(2.0))

from photutils import CircularAperture, aperture_photometry, CircularAnnulus
positions = [(64., 64.)]
aperture = CircularAperture(positions, r=25.)
annulus_aperture = CircularAnnulus(positions, r_in=25, r_out=35)

m51high = fits.open('/Users/mpeel/Desktop/SPLAT_M51high_SWT_I_AVERAGE_DENOISE.FITS')
print(m51high.info())
print(m51high[0].header['CDELT1']*60)
print(m51high[0].data)
print(np.shape(m51high[0].data[0]))
m51sum = np.nansum(m51high[0].data[0][m51high[0].data[0]>-0.0001])
print(m51sum)
print(np.nanmax(m51high[0].data[0]))
print(np.nanmin(m51high[0].data[0]))

phot_table = aperture_photometry(m51high[0].data[0], aperture)
print(phot_table)
bg_table = aperture_photometry(m51high[0].data[0], annulus_aperture)
print(bg_table)
exit()

beamsize = 0.8
beamarea = calc_beam_area(beamsize/60,beamsize/60)
print(beamarea)
pixelsize = 0.25
pixelarea = ((pixelsize/60)*(np.pi/180.0))**2
print(pixelarea)
print(pixelarea/beamarea)
print(m51sum*pixelarea/beamarea)
highval = m51sum*pixelarea/beamarea

m51low = fits.open('/Users/mpeel/Desktop/SPLAT_M51low_SWT_I_AVERAGE_DENOISE.FITS')
print(m51low.info())
print(m51low[0].header['CDELT1']*60)
print(m51low[0].data)
print(np.shape(m51low[0].data[0]))
m51sum = np.nansum(m51low[0].data[0])
print(m51sum)
print(np.nanmax(m51low[0].data[0]))

beamsize = 0.8
beamarea = calc_beam_area(beamsize/60,beamsize/60)
print(beamarea)
pixelsize = 0.25
pixelarea = ((pixelsize/60)*(np.pi/180.0))**2
print(pixelarea)
print(pixelarea/beamarea)
print(m51sum*pixelarea/beamarea)
lowval = m51sum*pixelarea/beamarea

print(np.log10(highval/lowval)/np.log10(24.6/18.6))

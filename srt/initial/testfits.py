#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit

col1 = np.arange(0,100,0.1)
col2 = np.arange(0,100,0.1)
col3 = np.arange(0,100,0.1)

col1 = fits.Column(name='ch1',format='E',array=np.asarray(ra).astype(float))
col2 = fits.Column(name='ch2',format='E',array=np.asarray(dec).astype(float))
col3 = fits.Column(name='ch3',format='E',array=np.asarray(timestream).astype(float))
=cols = fits.ColDefs([col1, col2, col3])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('test.fits',overwrite=True)

inputfits = fits.open('test.fits')
cols = inputfits[1].columns
print(cols.info())

col1_new = np.degrees(inputfits[1].data.field(3))
dec = np.degrees(inputfits[1].data.field(4))
timestream = inputfits[1].data.field(2)
print(ra)
print(dec)
print(timestream)
print(len(ra))
exit()

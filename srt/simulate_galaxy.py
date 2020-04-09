# Do a quick simulation of a galaxy with SRT observations

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from scipy.interpolate import interp2d

basedir='/Users/mpeel/Desktop/nearbygal/'
filename = 'NGC_6946_I_6cm_pb1991.fits'
# filename = 'NGC_6946_I_103aE_dss1.fits'
# filename = 'NGC_6946_NA_MOM0_I_HI_wbb2008.fits'
galaxy='NGC6946_0.3'
signallevel=280.0
signallevel_c = 476.0 
rescalemap = 1.0#18.0/25.6
dataformat = 0
declination = 0.0#60.0*np.pi/180.0

filename='NGC_891_RPcIm_I_5cm_w2015.fits'
galaxy='NGC891_0.3'
signallevel=100.0
signallevel_c = 203.4 
rescalemap = 1.0
dataformat = 0

filename='MESSIER_051_I_6cm_kwb1984.fits'
galaxy='M51_0.3'
signallevel=130.0
signallevel_c = 340.8
rescalemap = 1.0
dataformat = 0

# filename = 'NGC_4631_RPcIm_I_5cm_w2015.fits'
# galaxy = 'NGC 4631'
# signallevel = 280.0*0.72
# rescalemap = 1.0
# dataformat = 1

# filename = 'MESSIER_082_I_SPIRE_500_b2012.fits'
# filename = '2MASS_MESSIER_082_K.fits'
# galaxy = 'M82'
# signallevel = 1000.0
# rescalemap = 1.0
# dataformat = 1

# filename = 'NGC_5457.xmo9_I_HI_b1995.fits'
# filename = 'NGC_5457_I_250um_hipe_k2011.fits'
# # filename = 'NGC_5457_I_MIPS70_d2009.fits'
# galaxy = 'M101'
# signallevel = 140.0
# signallevel_c = 242.0
# rescalemap = 1.0
# dataformat = 1

# filename = 'NGC_2403_I_SPIRE_250_b2012.fits'
# filename = 'NGC_2403_RO_MOM0_I_HI_wbb2008.fits'
# filename = 'NGC_2403.xmi65_I_HI_b1995.fits'
# galaxy = 'NGC2403__h1'
# signallevel = 100.0
# signallevel_c = 150.0
# rescalemap = 1.0
# dataformat = 2

# C-Band
# beamsize = 2.7
# noiselevel=4.0
# signallevel = signallevel_c
# galaxy = galaxy+'_C'
# K-band
beamsize = 0.8
noiselevel=0.3

# Perseus for test
# filename='perseus_test.fits'
# galaxy='Perseus'
# signallevel=400.0e3
# rescalemap = 1.0
# dataformat = 1
# beamsize = 0.8
# noiselevel=13

hdul = fits.open(basedir+filename)
wcs = WCS(hdul[0].header)
print(hdul.info())
print(hdul[0].header)
try:
	print(hdul[0].header['CDELT1']*60.0*1024*rescalemap)
except:
	print(hdul[0].header['CD1_1']*60.0*1024*rescalemap)

# exit()
print(np.shape(hdul[0].data))
if dataformat == 0:
	data = np.asarray(hdul[0].data[0][0])
	# data = np.flip(data,axis=1)
	# data = np.flip(data,axis=0)
	slices = ('x','y',0,0)
elif dataformat == 2:
	data = np.asarray(hdul[0].data[0])
	# data = np.flip(data,axis=1)
	slices = ('x','y',0)
	# data = np.flip(data,axis=0)
else:
	data = np.asarray(hdul[0].data)
	# data = np.flip(data,axis=1)
	data[data < 0] = 0.0
	data[np.isfinite(data)==False] = 0.0
	slices = ('x','y')
plt.subplot(projection=wcs,slices=slices)
plt.imshow(data)
plt.colorbar()
plt.savefig(basedir+galaxy+'_sim.png')
plt.clf()

# if dataformat == 0:
try:
	val = int(np.abs(beamsize/(hdul[0].header['CDELT1']*60.0*rescalemap)))
except:
	val = int(np.abs(beamsize/(hdul[0].header['CD1_1']*60.0*rescalemap)))
print(val)
gauss = Gaussian2DKernel(val)
result = convolve_fft(data, gauss)
# else:
# 	result = data
plt.subplot(projection=wcs,slices=slices)
plt.imshow(result)
plt.colorbar()
plt.savefig(basedir+galaxy+'_sim_smth.png')
plt.clf()

W, H = result.shape[:2]
# if dataformat == 0:
new_W, new_H = (int(W/val),int(H/val))
print(W)
print(H)
print(new_W)
print(new_H)
xrange = lambda x: np.linspace(0, 1, x)
f = interp2d(xrange(H), xrange(W), result, kind="linear")
new_arr = f(xrange(new_H), xrange(new_W))
# else:
# 	new_W = W
# 	new_H = H
# 	new_arr = result
# rescale = ((val*W/new_W)**2)*280.0/np.sum(result)
rescale = signallevel/np.sum(new_arr)
print(rescale)
new_arr *= rescale
print(new_W*new_H)
new_arr[new_arr<0]=0
newheader = hdul[0].header
try:
	newheader['CDELT1'] *= W/new_W
	newheader['CDELT2'] *= H/new_H
except:
	newheader['CD1_1'] *= W/new_W
	newheader['CD2_2'] *= H/new_H
# newheader['']
# wcs.reset_wcs()
wcs = WCS(newheader)


plt.subplot(projection=wcs,slices=slices)
plt.axis('off')
plt.tight_layout()
plt.imshow(new_arr)
# plt.title(galaxy + ' convolved')
cbar = plt.colorbar()
cbar.set_label('mJy')
plt.savefig(basedir+galaxy+'_sim_smth_resample.png')
plt.clf()

newarr_size = int(20/0.8)
new_arr2 = np.zeros((newarr_size,newarr_size))
extra_width=int((newarr_size-new_W)/2)
extra_height=int((newarr_size-new_H)/2)
for i in range(new_W):
	for j in range(new_H):
		new_arr2[i+extra_width][j+extra_height] = new_arr[i][j]
plt.subplot(projection=wcs,slices=slices)
plt.rcParams["font.family"] = "serif"
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.axis('off')
plt.tight_layout()
plt.imshow(new_arr2)
# plt.title(galaxy + ' convolved')
cbar = plt.colorbar()
cbar.set_label('mJy',fontsize=20)
plt.savefig(basedir+galaxy+'_sim_smth_resample_20.png')
plt.clf()
print('Number of beams with S/N over 3:')
print(np.sum(new_arr2 > 3*noiselevel))
print('Number of beams with S/N over 5:')
print(np.sum(new_arr2 > 5*noiselevel))
print('Number of beams with S/N over 10:')
print(np.sum(new_arr2 > 10*noiselevel))

noise = np.random.normal(scale=1.0, size=new_W*new_H)*noiselevel
for i in range(0,new_W):
	for j in range(0,new_H):
		new_arr[i][j] += noise[i*new_W+j]

plt.subplot(projection=wcs,slices=slices)
# plt.title(galaxy + ' convolved and with noise')
plt.axis('off')
plt.tight_layout()
plt.imshow(new_arr)
cbar = plt.colorbar()
cbar.set_label('mJy',fontsize=14)
plt.savefig(basedir+galaxy+'_sim_smth_noise.png')
plt.clf()

noise = np.random.normal(scale=1.0, size=newarr_size*newarr_size)*noiselevel
for i in range(0,newarr_size):
	for j in range(0,newarr_size):
		new_arr2[i][j] += noise[i*newarr_size+j]

plt.subplot(projection=wcs,slices=slices)
# plt.title(galaxy + ' convolved and with noise')
plt.rcParams["font.family"] = "serif"
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.axis('off')
plt.tight_layout()
plt.imshow(new_arr2)
cbar = plt.colorbar()
cbar.set_label('mJy',fontsize=20)
plt.savefig(basedir+galaxy+'_sim_smth_noise_20.png')
plt.clf()

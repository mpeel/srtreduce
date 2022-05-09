# 2020-07-13 MP Start code to calculate moon positions for SRT for a given list of times.
# 2020-07-13 MP Rewrite (based on read_srt_data.py) to use times from the fits files.
import astropy as ap
import astropy.coordinates
import astropy.units as u
import numpy as np
import os
import astropy.io.fits as fits

def get_time_from_fits(filename):
	inputfits = fits.open(filename)
	return inputfits[4].data.field('time')

srt = ap.coordinates.EarthLocation(lat=39.4930*u.deg, lon=9.2451*u.deg, height=600*u.m)

# Look through the folders and find the moon data folders
basedir = '/Volumes/Toshiba5TB2/SRT/TauA/sardaraData/'
folderlist = os.listdir(basedir)
todolist = []
for folder in folderlist:
	if '.' not in folder:
		subfolderlist = os.listdir(basedir+folder)
		search = 'MOON'
		subfolderlist = [f for f in subfolderlist if search in f]
		for subfolder in subfolderlist:
			todolist.append(basedir+folder+'/'+subfolder)
		print(subfolderlist)
print(folderlist)
print(todolist)

# Now loop through the folders for the moon data
for inputdir in todolist:
	prefix=inputdir.replace(basedir,'').replace('/','_')
	ext = 'fits'
	inputlist = os.listdir(inputdir)
	print(inputlist)
	filelist = [f for f in inputlist if ext in f]
	trip = 0
	for filename in filelist:#[0:10]:
		if 'summary' in filename:
			continue
		print(filename)
		# Read in the time
		timelist = get_time_from_fits(inputdir+'/'+filename)
		times = ap.time.Time(timelist, format='mjd', scale='utc')
		# print(timelist)
		# Convert it
		bodypos = ap.coordinates.get_body('moon',times,srt)
		bodypos_azel = bodypos.transform_to(ap.coordinates.AltAz(location=srt,obstime=times))

		# Write it out
		np.savetxt('moonpos/'+filename.replace('.fits','_azelradec.dat'),np.transpose([times.jd,bodypos_azel.az.deg,bodypos_azel.alt.deg,bodypos.ra.deg,bodypos.dec.deg]),fmt='%10.8f',header='# JD	Az[deg]	El[deg]	RA[deg]	Dec[deg]')

# This was for the filelist
# filename = '20190215-214136-42-18-MOON_CROSS1_jd.dat'
# timelist = np.loadtxt(filename)
# times = ap.time.Time(timelist, format='jd', scale='utc')

# bodypos = ap.coordinates.get_body('moon',times,srt)
# bodypos_azel = bodypos.transform_to(ap.coordinates.AltAz(location=srt,obstime=times))

# np.savetxt(filename.replace('_jd.dat','_azelradec.dat'),np.transpose([times.jd,bodypos_azel.az.deg,bodypos_azel.alt.deg,bodypos.ra.deg,bodypos.dec.deg]),fmt='%10.8f',header='# JD	Az[deg]	El[deg]	RA[deg]	Dec[deg]')

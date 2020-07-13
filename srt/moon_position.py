# 2020-07-13 MP Start code to calculate moon positions for SRT for a given list of times.
import astropy as ap
import astropy.coordinates
import astropy.units as u
import numpy as np

srt = ap.coordinates.EarthLocation(lat=39.4930*u.deg, lon=9.2451*u.deg, height=600*u.m)
filename = '20190215-214136-42-18-MOON_CROSS1_jd.dat'
timelist = np.loadtxt(filename)
times = ap.time.Time(timelist, format='jd', scale='utc')

bodypos = ap.coordinates.get_body('moon',times,srt)
bodypos_azel = bodypos.transform_to(ap.coordinates.AltAz(location=srt,obstime=times))

np.savetxt(filename.replace('_jd.dat','_azelradec.dat'),np.transpose([times.jd,bodypos_azel.az.deg,bodypos_azel.alt.deg,bodypos.ra.deg,bodypos.dec.deg]),fmt='%10.8f',header='# JD	Az[deg]	El[deg]	RA[deg]	Dec[deg]')

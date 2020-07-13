import astropy as ap
import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
import matplotlib.pyplot as plt
import datetime
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)


srt = EarthLocation(lat=39.4928*u.deg, lon=9.245*u.deg, height=600*u.m)
utc_plus_one_hour = TimezoneInfo(utc_offset=1*u.hour)

source_names = ['ngc6946','m51','ngc891','3C48','3C286','3C138']
source_pos = []
for i in range(0,len(source_names)):
	source_pos.append(SkyCoord.from_name(source_names[i]))
	print(source_pos[i].to_string(style='hmsdms'))
# exit()
starttime = [Time('2020-11-24 10:37:00'),Time('2020-11-25 10:34:00'),Time('2020-11-27 03:26:00'),Time('2020-11-28 03:22:00'),Time('2020-11-30 17:30:00'),Time('2020-12-03 04:00:00'),Time('2020-12-03 13:45:00'),Time('2020-12-03 16:00:00'),Time('2020-12-04 15:56:00')]
endtime = [Time('2020-11-24 20:07:00'),Time('2020-11-25 20:04:00'),Time('2020-11-27 12:56:00'),Time('2020-11-28 12:52:00'),Time('2020-11-30 19:30:00'),Time('2020-12-03 06:00:00'),Time('2020-12-03 14:45:00'),Time('2020-12-04 01:30:00'),Time('2020-12-05 01:26:00')]
sampling = 1000
formatter = DateFormatter('%H:%M')
for k in range(0,len(starttime)):
	print('\n')
	print(starttime[k])
	fig, ax = plt.subplots(figsize=(8, 6))
	time_delta = (endtime[k].mjd - starttime[k].mjd)/sampling
	times_mjd = np.arange(starttime[k].mjd,endtime[k].mjd,time_delta)
	times = Time(times_mjd,format='mjd')

	sun = get_sun(times)
	altaz = sun.transform_to(AltAz(obstime=times,location=srt))
	plt.plot_date(times.plot_date,altaz.alt.deg,label='Sun',fmt='-')
	moon = get_moon(times)
	altaz = moon.transform_to(AltAz(obstime=times,location=srt))
	plt.plot_date(times.plot_date,altaz.alt.deg,label='Moon',fmt='-')

	for i in range(0,len(source_names)):
		# print(times)
		print(source_names[i])
		el40 = 0
		el85 = 0
		elm85 = 0
		elm40 = 0
		altaz = source_pos[i].transform_to(AltAz(obstime=times,location=srt))
		plt.plot_date(times.plot_date,altaz.alt.deg,label=source_names[i],fmt='-')
		# print(np.max(altaz.alt))
		# exit()
		for j in range(0,sampling):
			# print(altaz.alt[j].deg)
			if altaz.alt[j].deg > 40 and el40 == 0:
				print('el40 up: ' + str(times[j].datetime))
				el40 = 1
			if altaz.alt[j].deg > 85 and el85 == 0:
				print('el85 up: ' + str(times[j].datetime))
				el85 = 1
			if altaz.alt[j].deg < 85 and el85 == 1 and elm85 == 0:
				print('el85 down: ' + str(times[j].datetime))
				elm85 = 1
			if altaz.alt[j].deg < 40 and el40 == 1 and elm40 == 0:
				print('el40 down: ' + str(times[j].datetime))
				elm40 = 1
	plt.legend()
	plt.ylim([0,90])
	plt.xlabel('Time (UTC)')
	plt.ylabel('Elevation (deg)')
	plt.title(starttime[k].to_value('iso', subfmt='date'))
	ax.xaxis.set_major_formatter(formatter)
	plt.savefig('visibility_latest_'+str(k)+'.png')
	plt.clf()
	# exit()
		# print(altaz)
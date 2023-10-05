from glob import glob
# from casatasks import importfits, imstat
import numpy as np
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
from astropy import stats
import os
from scipy.ndimage import binary_dilation  


which = 'gaussians'
# which = 'disks'

dir_sim = f'./../data/{which}_input/'
dir_obs = f'./../data/{which}_observed/'

files_sim = glob(f'{dir_sim}*.fits')
files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix')

files_sim.sort()
files_obs.sort()

sum_sim = ['']*len(files_sim)
sum_obs = ['']*len(files_sim)
rms_arr = ['']*len(files_sim)
conf_arr = ['']*len(files_sim)
wide_arr = ['']*len(files_sim)

for i, file_sim in enumerate(files_sim):

	conf = file_sim.split('/')[-1].split('_')[0]
	wide = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')
	
	for file_obs in files_obs:
		if (conf in file_obs) & (wide in file_obs):
			
			conf_arr[i] = conf
			wide_arr[i] = wide

			print(file_sim.split('/')[-1], file_obs.split('/')[-1])
			
			if not os.path.exists(file_obs.replace('.Jyperpix', '.Jyperpix.fits')):
				print('[INFO] Importing file...') 
				exportfits(file_obs, file_obs.replace('.Jyperpix', '.Jyperpix.fits'))

			file_obs = file_obs.replace('.Jyperpix', '.Jyperpix.fits')

			data_sim = fits.getdata(file_sim)
			data_obs = fits.getdata(file_obs)

			rms_sim = 0 
			rms_obs = stats.mad_std(data_obs, ignore_nan=True)
			rms_obs = stats.mad_std(data_obs[data_obs<rms_obs], ignore_nan=True)
			rms_arr[i] = rms_obs*u.Jy

			# mask_high = data_obs>rms_obs*5
			# mask_low = data_obs>rms_obs*0

			mask_high = data_obs == np.nanmax(data_obs)
			mask_low = data_obs > 0

			mask = binary_dilation(mask_high, mask=mask_low, iterations=-1)

			rms_hdu = fits.PrimaryHDU(mask*1, fits.getheader(file_obs))
			rms_hdu.writeto(file_obs.replace('.Jyperpix.fits', '.Jyperpix.mask.fits'), overwrite=True)

			sum_sim[i] = np.nansum(data_sim)*u.Jy
			sum_obs[i] = np.nansum(data_obs[mask])*u.Jy

for i in range(len(sum_sim)):
	if sum_sim[i] == '':
		sum_sim[i] = np.nan *u.Jy
		sum_obs[i] = np.nan *u.Jy
		rms_arr[i] = np.nan *u.Jy

table = QTable([conf_arr, wide_arr, sum_sim, sum_obs, rms_arr], 
		names=('conf', 'wide', 'sum_sim', 'sum_obs', 'rms_obs'))

conf_ = table['conf'].copy()
wide_ = table['wide'].copy()

for i, (conf, wide) in enumerate(zip(conf_, wide_)):

	if conf == '':
		conf_[i] = np.nan
		wide_[i] = np.nan
	else:
		conf_[i] = float(conf.replace('conf', ''))
		wide_[i] = float(wide.replace('mrs0', ''))

table['conf_'] = np.array(conf_, dtype=float)
table['wide_'] = np.array(wide_, dtype=float)

table['sum_obs/sum_sim'] = table['sum_obs']/table['sum_sim'] # add ratio
table.sort(['conf_', 'wide_'], reverse=False) # sort

ids_nan = np.where(np.isnan(table['conf_']))[0]
table.remove_rows(ids_nan)

table.write(f'table_{which}.fits', overwrite=True)
table.write(f'table_{which}.csv', overwrite=True)
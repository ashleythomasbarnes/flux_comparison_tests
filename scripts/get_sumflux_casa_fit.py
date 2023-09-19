from glob import glob
import numpy as np
from astropy.table import QTable
import astropy.units as u


which = 'gaussians'
# which = 'disks'

dir_sim = f'./../data/{which}_input/'
dir_obs = f'./../data/{which}_observed/'

files_sim = glob(f'{dir_sim}*.fits*')
files_obs = glob(f'{dir_obs}*.image*')

files_sim.sort()
files_obs.sort()

sum_sim = ['']*len(files_sim)
sum_obs = ['']*len(files_sim)
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
			
			if not os.path.exists(file_sim.replace('.fits', '.img')): 
				importfits(file_sim, file_sim.replace('.fits', '.img'), overwrite=True, beam=['0arcsec', '0arcsec', '0deg'])

			file_sim = file_sim.replace('.fits', '.img')

			fit_sim = imfit(file_sim)
			fit_obs = imfit(file_obs)

			sum_sim[i] =  fit_sim['deconvolved']['component0']['flux']['value'][0] *u.Jy
			sum_obs[i] =  fit_obs['deconvolved']['component0']['flux']['value'][0] *u.Jy

table = QTable([conf_arr, wide_arr, sum_sim, sum_obs], 
		names=('conf', 'wide', 'sum_sim', 'sum_obs'))

table['sum_obs/sum_sim'] = table['sum_obs']/table['sum_sim'] # add ratio
table.sort(['conf', 'wide']) # sort

table.write(f'table_fit_{which}.fits', overwrite=True)
table.write(f'table_fit_{which}.csv', overwrite=True)
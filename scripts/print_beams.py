from glob import glob
import numpy as np
from astropy.io import fits

which = 'gaussians'

dir_sim = f'./../data/{which}_input/'
dir_obs = f'./../data/{which}_observed/'

files_sim = glob(f'{dir_sim}*.fits')
files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix.fits')

files_sim.sort()
files_obs.sort()

for i, file_sim in enumerate(files_sim):

    conf = file_sim.split('/')[-1].split('_')[0]
    wide = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')
    
    for file_obs in files_obs:
        if (conf in file_obs) & (wide in file_obs):

            header_sim = fits.getheader(file_sim)
            header_obs = fits.getheader(file_obs)

            bmin = header_obs['BMAJ']
            bmaj = header_obs['BMIN']
            bpa = header_obs['BPA']

            # print('%s %s \t %0.2farcsec %0.2farcsec %0.1fdeg' %(conf, wide, bmin*3600, bmaj*3600, bpa))
            print('%s %s %0.2f %0.2f %0.1f' %(conf.replace('conf', ''), wide.replace('mrs0', ''), bmin*3600, bmaj*3600, bpa))
from glob import glob
import numpy as np
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
from astropy import stats
from scipy.ndimage import binary_dilation  
from astropy.modeling import models, fitting
import numpy as np
from astropy.modeling import models, fitting
from reproject import reproject_interp
import time 
import warnings 
warnings.filterwarnings('ignore')

def fit_2d_gaussian_and_get_sum(image):

    image = np.squeeze(image)
    image[np.isnan(image)] = 0

    # Get the center of the image
    shape_x, shape_y = image.shape
    center_x, center_y = np.array(image.shape) // 2
    
    # Create x, y indices grid for the sub-image
    x, y = np.array(np.mgrid[:shape_x, :shape_y], dtype=np.int32)
    
    # Initialize the Gaussian2D model
    g_init = models.Gaussian2D(amplitude=np.nanmax(image), x_mean=center_x, y_mean=center_y)
    
    # Fit the model to the sub-image
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y, image)
    
    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g(x, y), dtype=np.float32)
    gaussian_sum = np.sum(fitted_data)
    
    return (gaussian_sum*u.Jy, fitted_data)

def remove_padding(data):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return data
    
    # Crop the data array
    cropped_data = data[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    return cropped_data

def remove_padding_2arrays(data1, data2):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data1, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data1, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return data1
    
    # Crop the data array
    cropped_data = data2[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    return cropped_data

which = 'gaussians'

which_time = ''
# which_time = '_6totaltime'
# which_time = '_60totaltime'
# which_time = '_6totaltime_flagged'
# which_time = '_60totaltime_flagged'

dir_sim = f'./../data/{which}_input/'
dir_obs = f'./../data/{which}_observed{which_time}/'

files_sim = glob(f'{dir_sim}*.fits')
files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix.fits')

files_sim.sort()
files_obs.sort()

confs = ['conf1']

start_total = time.time()

for conf in confs: 

    files_sim = glob(f'{dir_sim}{conf}_*.fits')
    files_obs = glob(f'{dir_obs}{conf}_*.pbcor.Jyperpix.fits')
    files_sim.sort()
    files_obs.sort()

    n = 15
    max_sim = ['']*n
    max_obs = ['']*n
    sum_sim = ['']*n
    sum_obs = ['']*n
    sum_mask10_obs = ['']*n
    sum_mask10_sim = ['']*n
    sum_mask50_obs = ['']*n
    sum_mask50_sim = ['']*n
    sum_fit_sim = ['']*n
    sum_fit_obs = ['']*n
    rms_arr = ['']*n
    conf_arr = ['']*n
    wide_arr = ['']*n

    for i, file_sim in enumerate(files_sim):
        
        start = time.time()

        conf_file = file_sim.split('/')[-1].split('_')[0]
        wide_file = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')
        
        for file_obs in files_obs:
            if (conf_file in file_obs) & (wide_file in file_obs):

                if conf_file != conf: 
                    continue

                conf_arr[i] = conf
                wide_arr[i] = wide_file

                print('[INFO] File sim: %s' %file_sim.split('/')[-1])
                print('[INFO] File obs: %s' %file_obs.split('/')[-1])
                hdu_sim = fits.open(file_sim)[0]
                hdu_obs = fits.open(file_obs)[0]
                hdu_sim.data = np.array(hdu_sim.data, dtype=np.float32)
                hdu_obs.data = np.array(hdu_obs.data, dtype=np.float32)

                print('[INFO] Getting squeeze...')
                data_sim = np.array(np.squeeze(hdu_sim.data.copy()), dtype=np.float32)
                data_obs = np.array(np.squeeze(hdu_obs.data.copy()), dtype=np.float32)
                print('[INFO] Getting reproject...')
                data_sim_r, _ = np.array(reproject_interp(hdu_sim, hdu_obs.header), dtype=np.float32)
                # data_sim_r, _ = np.array(reproject_interp(hdu_sim, hdu_obs.header, parallel=True), dtype=np.float32)
                data_sim_r2 = data_sim_r.copy()
                print('[INFO] Time to run reproject: %0.2f [sec]' %((time.time() - start)))

                # print(data_sim.shape)
                # print(data_obs.shape)
                # print(np.squeeze(data_sim_r).shape)
                print('[INFO] Getting cropping...')
                data_sim = remove_padding(data_sim)
                data_obs = remove_padding(data_obs)
                data_sim_r = remove_padding(np.squeeze(data_sim_r))
                print('[INFO] Time to run cropping: %0.2f [sec]' %((time.time() - start)))
                # print(data_sim.shape)
                # print(data_obs.shape)
                # print(np.squeeze(data_sim_r).shape)

                if np.squeeze(data_sim_r).shape[0] < data_obs.shape[0]:
                    print('[INFO] Getting cropping - data to sim...')
                    data_obs = np.array(np.squeeze(hdu_obs.data.copy()), dtype=np.float32)
                    data_obs = remove_padding_2arrays(np.squeeze(data_sim_r), data_obs)
                    print(data_obs.shape)
                    print('[INFO] Time to run data to sim cropping: %0.2f [sec]' %((time.time() - start)))

                print('[INFO] Getting sums...')
                rms_sim = 0  
                rms_obs = stats.mad_std(data_obs, ignore_nan=True)
                rms_obs = stats.mad_std(data_obs[data_obs<rms_obs], ignore_nan=True)
                rms_arr[i] = rms_obs*u.Jy

                mask_high = data_obs == np.nanmax(data_obs)
                mask_low = data_obs > rms_obs*3
                mask = binary_dilation(mask_high, mask=mask_low, iterations=-1)

                # mask_hdu = fits.PrimaryHDU(mask*np.int32(1), fits.getheader(file_obs))
                # mask_hdu.writeto(file_obs.replace('.Jyperpix.fits', '.Jyperpix.mask.fits'), overwrite=True)

                max_sim[i] = np.nanmax(data_sim)*u.Jy
                max_obs[i] = np.nanmax(data_obs)*u.Jy

                sum_sim[i] = np.nansum(data_sim)*u.Jy
                sum_obs[i] = np.nansum(data_obs[mask])*u.Jy
                print('[INFO] Time to run sums: %0.2f [sec]' %((time.time() - start)))

                print('[INFO] Getting fits...')
                sum_fit_sim[i], _ = fit_2d_gaussian_and_get_sum(data_sim)
                print('[INFO] Time to run sim fit: %0.2f [sec]' %((time.time() - start)))
                sum_fit_obs[i], _ = fit_2d_gaussian_and_get_sum(data_obs) 
                print('[INFO] Time to run obs fit: %0.2f [sec]' %((time.time() - start)))

                print('[INFO] Getting masked sums...')
                ### Get masked flux
                mask_sim = hdu_sim.data/np.nanmax(hdu_sim.data) > 0.1
                mask_sim_r = data_sim_r2/np.nanmax(data_sim_r2) > 0.1
                
                sum_mask10_sim[i] = np.nansum(hdu_sim.data[mask_sim])*u.Jy
                sum_mask10_obs[i] = np.nansum(hdu_obs.data[mask_sim_r])*u.Jy

                mask_sim = hdu_sim.data/np.nanmax(hdu_sim.data) > 0.5
                mask_sim_r = data_sim_r2/np.nanmax(data_sim_r2) > 0.5
                
                sum_mask50_sim[i] = np.nansum(hdu_sim.data[mask_sim])*u.Jy
                sum_mask50_obs[i] = np.nansum(hdu_obs.data[mask_sim_r])*u.Jy
                print('[INFO] Time to run masked sums: %0.2f [sec]' %((time.time() - start)))

                print('[INFO] Time to run: %0.2f [sec]' %((time.time() - start)))

    ### Organise into astropy table 
    for i in range(len(sum_sim)):
        if sum_mask10_sim[i] == '':
            sum_fit_sim[i] = np.nan *u.Jy
            sum_fit_obs[i] = np.nan *u.Jy
            sum_mask10_sim[i] = np.nan *u.Jy
            sum_mask10_obs[i] = np.nan *u.Jy
            sum_mask50_sim[i] = np.nan *u.Jy
            sum_mask50_obs[i] = np.nan *u.Jy
            sum_sim[i] = np.nan *u.Jy
            sum_obs[i] = np.nan *u.Jy
            rms_arr[i] = np.nan *u.Jy
            max_obs[i] = np.nan *u.Jy
            max_sim[i] = np.nan *u.Jy

    table = QTable([conf_arr, 
                wide_arr, 
                sum_sim, 
                sum_obs, 
                rms_arr, 
                max_sim,
                max_obs,
                sum_fit_sim, 
                sum_fit_obs,
                sum_mask10_sim, 
                sum_mask10_obs,
                sum_mask50_sim, 
                sum_mask50_obs,
                ], 
            names=('conf', 
                'wide', 
                'sum_sim', 
                'sum_obs', 
                'rms_obs', 
                'max_sim',
                'max_obs',
                'sum_fit_sim', 
                'sum_fit_obs',
                'sum_mask10_sim', 
                'sum_mask10_obs',
                'sum_mask50_sim', 
                'sum_mask50_obs'))

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
    table.sort(['conf_', 'wide_'], reverse=False) # sort

    ids_nan = np.where(np.isnan(table['conf_']))[0]
    table.remove_rows(ids_nan)

    table.write(f'../data/tables/table_fit_{which}{which_time}_{conf}.fits', overwrite=True)
    table.write(f'../data/tables/table_fit_{which}{which_time}_{conf}.csv', overwrite=True)

end_total = time.time()
print('[INFO] Time to run total: %0.2f [min]' %((end_total - start_total)/60))
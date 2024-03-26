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
from scipy.integrate import dblquad
from reproject import reproject_interp
import time 
import warnings 
warnings.filterwarnings('ignore')

def fit_2d_gaussian(image, x_stddev=None, y_stddev=None):

    image = np.squeeze(image)
    image[np.isnan(image)] = 0

    # Get the center of the image
    shape_x, shape_y = image.shape
    center_x, center_y = np.array(image.shape) // 2
    
    # Create x, y indices grid for the sub-image
    x, y = np.array(np.mgrid[:shape_x, :shape_y], dtype=np.int32)
    
    # Initialize the Gaussian2D model
    g_init = models.Gaussian2D(amplitude=np.nanmax(image), 
                               x_mean=center_x, y_mean=center_y, 
                               x_stddev=x_stddev, y_stddev=y_stddev)
    
    # Fit the model to the sub-image
    fit_g = fitting.LevMarLSQFitter(calc_uncertainties=True)
    g = fit_g(g_init, x, y, image)

    # Calculate the uncertainty on the fit using the covariance matrix
    cov_matrix = fit_g.fit_info['param_cov']
    if cov_matrix is not None:
        g_err = np.sqrt(np.diag(cov_matrix))
    else:
        g_err = np.zeros_like(g.parameters)

    g_errl = models.Gaussian2D(amplitude=g.amplitude.value - np.abs(g_err[0]), 
                        x_mean=g.x_mean.value, 
                        y_mean=g.y_mean.value, 
                        x_stddev=g.x_stddev.value -  np.abs(g_err[3]), 
                        y_stddev=g.y_stddev.value -  np.abs(g_err[4]),
                        theta=g.theta.value)

    g_errh = models.Gaussian2D(amplitude=g.amplitude.value + np.abs(g_err[0]), 
                                x_mean=g.x_mean.value, 
                                y_mean=g.y_mean.value, 
                                x_stddev=g.x_stddev.value + np.abs(g_err[3]), 
                                y_stddev=g.y_stddev.value + np.abs(g_err[4]),
                                theta=g.theta.value)
    
    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g(x, y), dtype=np.float32)
    sum_tot = np.nansum(fitted_data)

    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g_errl(x, y), dtype=np.float32)
    sum_errl = np.nansum(fitted_data)

    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g_errh(x, y), dtype=np.float32)
    sum_errh = np.nansum(fitted_data)

    # Compute the reduced chi-squared
    residuals = image - g(x, y)
    chi_squared = np.sum((residuals / image.std())**2)
    degrees_of_freedom = image.size - len(g.parameters)
    reduced_chi_squared = chi_squared / degrees_of_freedom

    return (g, g_err, (sum_tot*u.Jy, sum_errl*u.Jy, sum_errh*u.Jy), reduced_chi_squared)

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

# which_time = ''
which_time = '_6totaltime'
# which_time = '_60totaltime'
# which_time = '_6totaltime_flagged'
# which_time = '_60totaltime_flagged'

# which_times = ['_6totaltime', '_60totaltime', '_6totaltime_flagged', '_60totaltime_flagged']
which_times = ['_6totaltime']

for which_time in which_times:

    dir_sim = f'./../data/{which}_input/'
    dir_obs = f'./../data/{which}_observed{which_time}/'

    files_sim = glob(f'{dir_sim}*.fits')
    files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix.fits')

    files_sim.sort()
    files_obs.sort()

    # confs = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5', 'conf6', 'conf7', 'conf8', 'conf9', 'conf10']
    # confs = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5']
    confs = ['conf1']
    count = 0 

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
        errl_fit_sim = ['']*n
        errh_fit_sim = ['']*n
        errl_fit_obs = ['']*n
        errh_fit_obs = ['']*n
        rms_arr = ['']*n
        conf_arr = ['']*n
        wide_arr = ['']*n
        rchi2_sim = ['']*n
        rchi2_obs = ['']*n

        for i, file_sim in enumerate(files_sim):
            
            start = time.time()

            conf_file = file_sim.split('/')[-1].split('_')[0]
            wide_file = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')
            
            for file_obs in files_obs:
                if (conf_file in file_obs) & (wide_file in file_obs):

                    if conf_file != conf: 
                        continue

                    # count+=1
                    # if count!=12:
                    #     continue

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

                    print('[INFO] Getting cropping...')
                    data_sim = remove_padding(data_sim)
                    data_obs = remove_padding(data_obs)

                    hdu_sim.data = np.squeeze(hdu_sim.data)
                    hdu_obs.data = data_obs

                    hdu_obs.header['CRPIX1'] = hdu_obs.data.shape[1]/2
                    hdu_obs.header['CRPIX2'] = hdu_obs.data.shape[0]/2
                    del hdu_sim.header['*3*']
                    del hdu_sim.header['*4*']
                    del hdu_obs.header['*3*']
                    del hdu_obs.header['*4*']
                    del hdu_obs.header['*PC*']
                    del hdu_sim.header['*PC*']

                    data_sim_r, _ = np.array(reproject_interp(hdu_sim, hdu_obs.header), dtype=np.float32)
                    data_sim_rp = remove_padding(data_sim_r.copy())
                    if np.squeeze(data_sim_rp).shape[0] < data_obs.shape[0]:
                        print('[INFO] Getting cropping - data to sim...')
                        data_obs = np.array(np.squeeze(hdu_obs.data.copy()), dtype=np.float32)
                        data_obs = remove_padding_2arrays(np.squeeze(data_sim_r.copy()), data_obs)
                        print('[INFO] Time to run data to sim cropping: %0.2f [sec]' %((time.time() - start)))

                    print('[INFO] Getting sums...')
                    rms_sim = 0  
                    rms_obs = stats.mad_std(hdu_obs.data, ignore_nan=True)
                    rms_obs = stats.mad_std(hdu_obs.data[hdu_obs.data<rms_obs], ignore_nan=True)
                    rms_arr[i] = rms_obs*u.Jy

                    mask_high = data_obs == np.nanmax(data_obs)
                    mask_low = data_obs > rms_obs*3
                    mask = binary_dilation(mask_high, mask=mask_low, iterations=-1)
                    
                    max_sim[i] = np.nanmax(data_sim)*u.Jy
                    max_obs[i] = np.nanmax(data_obs)*u.Jy

                    sum_sim[i] = np.nansum(data_sim)*u.Jy
                    sum_obs[i] = np.nansum(data_obs[mask])*u.Jy
                    print('[INFO] Time to run sums: %0.2f [sec]' %((time.time() - start)))

                    print('[INFO] Getting fits...')
                    g, g_err, (sum_fit_sim[i], errl_fit_sim[i], errh_fit_sim[i]), rchi2_sim[i] = fit_2d_gaussian(data_sim)
                    print('[INFO] Time to run sim fit: %0.2f [sec]' %((time.time() - start)))
                    g, g_err, (sum_fit_obs[i], errl_fit_obs[i], errh_fit_obs[i]), rchi2_obs[i] = fit_2d_gaussian(data_obs)
                    print('[INFO] Time to run obs fit: %0.2f [sec]' %((time.time() - start)))

                    print('[INFO] Getting masked sums...')
                    mask_sim = data_sim/np.nanmax(data_sim) > 0.1
                    mask_sim_r = data_sim_rp/np.nanmax(data_sim_rp) > 0.1
                    
                    sum_mask10_sim[i] = np.nansum(data_sim[mask_sim])*u.Jy
                    sum_mask10_obs[i] = np.nansum(data_obs[mask_sim_r])*u.Jy

                    mask_sim = data_sim/np.nanmax(data_sim) > 0.5
                    mask_sim_r = data_sim_rp/np.nanmax(data_sim_rp) > 0.5
                    
                    sum_mask50_sim[i] = np.nansum(data_sim[mask_sim])*u.Jy
                    sum_mask50_obs[i] = np.nansum(data_obs[mask_sim_r])*u.Jy
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
                errl_fit_sim[i] = np.nan *u.Jy
                errh_fit_sim[i] = np.nan *u.Jy
                errl_fit_obs[i] = np.nan *u.Jy
                errh_fit_obs[i] = np.nan *u.Jy
                rms_arr[i] = np.nan *u.Jy
                max_obs[i] = np.nan *u.Jy
                max_sim[i] = np.nan *u.Jy
                rchi2_sim = np.nan *u.Jy/u.Jy
                rchi2_obs = np.nan *u.Jy/u.Jy

        table = QTable([conf_arr, 
                    wide_arr, 
                    sum_sim, 
                    sum_obs, 
                    rms_arr, 
                    max_sim,
                    max_obs,
                    sum_fit_sim, 
                    sum_fit_obs,
                    errl_fit_sim,
                    errh_fit_sim,
                    errl_fit_obs,
                    errh_fit_obs,
                    sum_mask10_sim, 
                    sum_mask10_obs,
                    sum_mask50_sim, 
                    sum_mask50_obs,
                    rchi2_sim,
                    rchi2_obs
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
                    'errl_fit_sim',
                    'errh_fit_sim',
                    'errl_fit_obs',
                    'errh_fit_obs',
                    'sum_mask10_sim', 
                    'sum_mask10_obs',
                    'sum_mask50_sim', 
                    'sum_mask50_obs',
                    'rchi2_sim',
                    'rchi2_obs'))

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
from glob import glob
# from casatasks import importfits, imstat
import numpy as np
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
from astropy import stats
import os
from scipy.ndimage import binary_dilation  
from astropy.modeling import models, fitting
import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

def plot_2d_gaussian(image, fitted_data, mask, conf, wide, outputfile=''):

    image = np.squeeze(image)
    image[np.isnan(image)] = 0

    # Get the center of the image
    shape_x, shape_y = image.shape
    center_x, center_y = np.array(image.shape) // 2
    
    # Create x, y indices grid for the sub-image
    x, y = np.array(np.mgrid[:shape_x, :shape_y], dtype=np.int32)

    # Plotting
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    
    ax[0].imshow(image, origin='lower', cmap='inferno')
    ax[1].imshow(fitted_data, origin='lower', cmap='inferno')
    
    # You might want to overlay contours of the fitted Gaussian over the original image
    # for a direct comparison. This would look something like this:
    ax[0].contour(y, x, fitted_data, colors='C0', levels=[0.5*np.nanmax(fitted_data)], alpha=0.6)
    ax[0].contour(y, x, mask, colors='C1', levels=[0], alpha=0.6)

    ax[1].contour(y, x, fitted_data, colors='C0', levels=[0.5*np.nanmax(fitted_data)], alpha=0.6)
    ax[1].contour(y, x, mask, colors='C1', levels=[0], alpha=0.6)
    ax[1].plot([0,0],[0,0], c='C0', label='fit FHWM')
    ax[1].plot([0,0],[0,0], c='C1', label='thresh mask')

    ax[0].text(0.5, 0.9, 'Original Image', color='white', fontweight='normal', size=12, transform=ax[0].transAxes, ha='center')
    ax[1].text(0.5, 0.9, 'Fitted Gaussian', color='white', fontweight='normal', size=12, transform=ax[1].transAxes, ha='center')
    ax[0].text(0.5, 0.95, '%s %s' %(conf, wide), color='white', fontweight='heavy', size=15, transform=ax[0].transAxes, ha='center')
    ax[1].text(0.5, 0.95, '%s %s' %(conf, wide), color='white', fontweight='heavy', size=15, transform=ax[1].transAxes, ha='center')

    ax[1].legend(loc='upper left')

    for a in ax:
        a.set_xlabel('x')
        a.set_ylabel('y')
    
    plt.tight_layout()
    plt.savefig(outputfile, dpi=300, bbox_inches='tight')

    return()


def remove_nan_padding(data):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return data
    
    # Crop the data array
    cropped_data = data[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
       
    return cropped_data

def binArray(data, axis, binstep=2, binsize=2, func=np.nansum):
    data = np.array(data)
    dims = np.array(data.shape)
    argdims = np.arange(data.ndim)
    argdims[0], argdims[axis]= argdims[axis], argdims[0]
    data = data.transpose(argdims)
    data = [func(np.take(data,np.arange(int(i*binstep),int(i*binstep+binsize)),0),0) for i in np.arange(dims[axis]//binstep)]
    data = np.array(data).transpose(argdims)
    return data

which = 'gaussians'

dir_sim = f'./../data/{which}_input/'
dir_obs = f'./../data/{which}_observed/'

files_sim = glob(f'{dir_sim}*.fits')
files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix')

files_sim.sort()
files_obs.sort()

max_sim = ['']*len(files_sim)
max_obs = ['']*len(files_sim)
sum_sim = ['']*len(files_sim)
sum_obs = ['']*len(files_sim)
sum_fit_sim = ['']*len(files_sim)
sum_fit_obs = ['']*len(files_sim)
rms_arr = ['']*len(files_sim)
conf_arr = ['']*len(files_sim)
wide_arr = ['']*len(files_sim)

for i, file_sim in enumerate(files_sim):

    conf = file_sim.split('/')[-1].split('_')[0]
    wide = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')
    
    for file_obs in files_obs:
        if (conf in file_obs) & (wide in file_obs):
            
            if float(conf.replace('conf', '')) not in [5]: 
                continue

            conf_arr[i] = conf
            wide_arr[i] = wide

            print(file_sim.split('/')[-1], file_obs.split('/')[-1])
            file_obs = file_obs.replace('.Jyperpix', '.Jyperpix.fits')

            data_sim = np.array(np.squeeze(fits.getdata(file_sim)), dtype=np.float32)
            data_obs = np.array(np.squeeze(fits.getdata(file_obs)), dtype=np.float32)

            data_obs = remove_nan_padding(data_obs)

            rms_sim = 0  
            rms_obs = stats.mad_std(data_obs, ignore_nan=True)
            rms_obs = stats.mad_std(data_obs[data_obs<rms_obs], ignore_nan=True)
            rms_arr[i] = rms_obs*u.Jy

            mask_high = data_obs == np.nanmax(data_obs)
            # mask_low = data_obs > 0
            mask_low = data_obs > rms_obs*3
            mask = binary_dilation(mask_high, mask=mask_low, iterations=-1)

            rms_hdu = fits.PrimaryHDU(mask*np.int32(1), fits.getheader(file_obs))
            rms_hdu.writeto(file_obs.replace('.Jyperpix.fits', '.Jyperpix.mask.fits'), overwrite=True)

            max_sim[i] = np.nanmax(data_sim)*u.Jy
            max_obs[i] = np.nanmax(data_obs)*u.Jy

            sum_sim[i] = np.nansum(data_sim)*u.Jy
            sum_obs[i] = np.nansum(data_obs[mask])*u.Jy

            sum_fit_sim[i], _ = fit_2d_gaussian_and_get_sum(data_sim)
            sum_fit_obs[i], fitted_data = fit_2d_gaussian_and_get_sum(data_obs) 

            plot_2d_gaussian(data_obs, fitted_data, mask*np.int32(1), conf, wide,
                             outputfile=file_obs.replace('.Jyperpix.fits', '.Jyperpix.png')) 

for i in range(len(sum_sim)):
    if sum_sim[i] == '':
        sum_fit_sim[i] = np.nan *u.Jy
        sum_fit_obs[i] = np.nan *u.Jy
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
            sum_fit_obs], 
        names=('conf', 
            'wide', 
            'sum_sim', 
            'sum_obs', 
            'rms_obs', 
            'max_sim',
            'max_obs',
            'sum_fit_sim', 
            'sum_fit_obs'))

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
table['sum_fit_obs/sum_fit_sim'] = table['sum_fit_obs']/table['sum_fit_sim'] # add ratio
table.sort(['conf_', 'wide_'], reverse=False) # sort

ids_nan = np.where(np.isnan(table['conf_']))[0]
table.remove_rows(ids_nan)

table.write(f'table_fit_{which}.fits', overwrite=True)
table.write(f'table_fit_{which}.csv', overwrite=True)
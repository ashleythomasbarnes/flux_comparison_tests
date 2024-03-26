from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy
import numpy as np
import astropy.units as u
import time 
import warnings 

plt.style.use('paper')
warnings.filterwarnings('ignore')

def custom_sort(arr, text='mrs0'):

    arr_new = [np.float32(item.replace(text,'')) for item in arr]
    arr_new = np.array(arr_new)

    arr_sort = np.sort(arr_new)
    arr_argsort = np.argsort(arr_new)

    return(arr_sort, arr_argsort)

def remove_padding(data):

    data = np.squeeze(data)
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data, axis=1)!=0)[0]
    
    # Crop the data array
    cropped_data = data[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    return cropped_data

#### 
which = 'gaussians'

# which_time = ''
which_time = '_6totaltime'
# which_time = '_60totaltime'
# which_time = '_6totaltime_flagged'
# which_time = '_60totaltime_flagged'

# which_times = ['_6totaltime', '_60totaltime', '_6totaltime_flagged', '_60totaltime_flagged']
which_times = ['_6totaltime']

bbox = dict(facecolor='whitesmoke',  alpha=0.95, boxstyle='round')
# confs = ['conf1', 'conf2', 'conf3', 'conf4', 'conf5', 'conf6', 'conf7', 'conf8', 'conf9', 'conf10']
confs = ['conf1']
### 

# Sims 
for which_time in which_times:

    dir_sim = f'./../data/{which}_input/'
    dir_obs = f'./../data/{which}_observed{which_time}/'

    files_sim = glob(f'{dir_sim}*.fits')
    files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix.fits')

    files_sim.sort()
    files_obs.sort()

    start_total = time.time()
    for conf in confs: 

        files_sim = glob(f'{dir_sim}{conf}_*.fits')
        files_obs = glob(f'{dir_obs}{conf}_*.pbcor.Jyperpix.fits')
        files_sim.sort()
        files_obs.sort()

        n = 15
        conf_arr = ['']*n
        wide_arr = ['']*n

        for i, file_sim in enumerate(files_sim):
            
            start = time.time()

            conf_file = file_sim.split('/')[-1].split('_')[0]
            wide_file = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')

            conf_arr[i] = conf
            wide_arr[i] = wide_file

        fig1 = plt.figure(figsize=(10, 6.5))
        
        title1 = which_time.replace('_', ' ')[1:]
        title2 = conf_arr[0]
        fig1.suptitle('%s %s' %(title1, title2), fontweight='bold', y=0.95)

        wide_arr_sort, wide_arr_argsort = custom_sort(wide_arr)

        for i, j in enumerate(wide_arr_argsort):
            
            start = time.time()

            file_sim = files_sim[j]

            print('[INFO] File sim: %s' %file_sim.split('/')[-1])

            hdu_sim = fits.open(file_sim)[0]
            hdu_sim.data = np.squeeze(hdu_sim.data)
            del hdu_sim.header['*3*'] 
            del hdu_sim.header['*4*'] 

            ax1 = aplpy.FITSFigure(hdu_sim, subplot=(3, 5, i+1), figure=fig1)

            vmin, vmax = np.nanpercentile(hdu_sim.data, (0.1, 99.9))
            ax1.show_colorscale(vmin=vmin, vmax=vmax, cmap='inferno')

            # ax1.recenter(ra, dec, size)
            ax1.axis_labels.hide()
            ax1.tick_labels.hide()

            label = 'conf'
            ax1.add_label(0.05, 0.95, wide_arr[j],  ha='left', va='top', size=10, bbox = bbox, relative=True)

            ax1.set_nan_color('lightgrey')
            ax_plot = fig1.get_axes()[-1]
            ax_plot.grid(True, alpha=0.3, ls=':', color='white')

        fig1.tight_layout(h_pad=0, w_pad=0)
        fig1.savefig(f'./../figs/maps_{which}{which_time}_{conf_arr[i]}_sim.pdf', dpi=300, bbox_inches='tight', transparent=False)
        fig1.savefig(f'./../figs/maps_{which}{which_time}_{conf_arr[i]}_sim.png', dpi=300, bbox_inches='tight', transparent=False)
        plt.close('all')

# Obs
for which_time in which_times:

    dir_sim = f'./../data/{which}_input/'
    dir_obs = f'./../data/{which}_observed{which_time}/'

    files_sim = glob(f'{dir_sim}*.fits')
    files_obs = glob(f'{dir_obs}*.pbcor.Jyperpix.fits')

    files_sim.sort()
    files_obs.sort()

    start_total = time.time()
    for conf in confs: 

        files_sim = glob(f'{dir_sim}{conf}_*.fits')
        files_obs = glob(f'{dir_obs}{conf}_*.pbcor.Jyperpix.fits')
        files_sim.sort()
        files_obs.sort()

        n = 15
        conf_arr = ['']*n
        wide_arr = ['']*n

        for i, file_sim in enumerate(files_sim):
            
            start = time.time()

            conf_file = file_sim.split('/')[-1].split('_')[0]
            wide_file = file_sim.split('/')[-1].split('_')[-1].replace('.fits', '')

            conf_arr[i] = conf
            wide_arr[i] = wide_file

        fig1 = plt.figure(figsize=(10, 6.5))
        
        title1 = which_time.replace('_', ' ')[1:]
        title2 = conf_arr[0]
        fig1.suptitle('%s %s' %(title1, title2), fontweight='bold', y=0.95)

        wide_arr_sort, wide_arr_argsort = custom_sort(wide_arr)

        for i, j in enumerate(wide_arr_argsort):
            
            start = time.time()

            file_obs = files_obs[j]

            print('[INFO] File sim: %s' %file_obs.split('/')[-1])

            hdu_obs = fits.open(file_obs)[0]
            hdu_obs.data = remove_padding(hdu_obs.data)
            del hdu_obs.header['*3*'] 
            del hdu_obs.header['*4*'] 

            ax1 = aplpy.FITSFigure(hdu_obs, subplot=(3, 5, i+1), figure=fig1)

            vmin, vmax = np.nanpercentile(hdu_obs.data, (0.1, 99.9))
            ax1.show_colorscale(vmin=vmin, vmax=vmax, cmap='inferno')

            # ax1.recenter(ra, dec, size)
            ax1.axis_labels.hide()
            ax1.tick_labels.hide()

            ax1.add_label(0.05, 0.95, wide_arr[j],  ha='left', va='top', size=10, bbox = bbox, relative=True)

            ax1.set_nan_color('lightgrey')
            ax_plot = fig1.get_axes()[-1]
            ax_plot.grid(True, alpha=0.3, ls=':', color='white')

        fig1.tight_layout(h_pad=0, w_pad=0)
        fig1.savefig(f'./../figs/maps_{which}{which_time}_{conf_arr[i]}_obs.pdf', dpi=300, bbox_inches='tight', transparent=False)
        fig1.savefig(f'./../figs/maps_{which}{which_time}_{conf_arr[i]}_obs.png', dpi=300, bbox_inches='tight', transparent=False)
        plt.close('all')
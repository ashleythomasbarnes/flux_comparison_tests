from glob import glob
import numpy as np
from astropy.table import QTable
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm

def custom_sort(arr):
    def key_func(item):
        # Extract number from string and convert to integer
        num = int(item[4:])  # 'conf' has length 4
        # Return a tuple with two items:
        # First item is a boolean indicating if number is 10 (True for 'conf10' and False otherwise)
        # Second item is the number itself
        # This will sort 'conf10' last and others based on their numbers
        return (num == 10, num)
    
    return sorted(arr, key=key_func)

# which = 'gaussians'
which = 'disks' 

table = QTable.read(f'table_fit_{which}.fits')
conf_unique = np.unique(table['conf']) 
conf_unique =  custom_sort(conf_unique)

fig1, ax1 = plt.subplots(2, 5, figsize=(15, 5), sharex=True)
fig2, ax2 = plt.subplots(2, 5, figsize=(15, 5), sharex=True)

ax1 = ax1.flatten()
ax2 = ax2.flatten()

colors = cm.turbo(np.linspace(0, 1, len(ax1)))

for i, conf in enumerate(conf_unique):

	conf_tab = table[np.where(table['conf']==conf)]
	conf_arr = conf_tab['conf']
	wide_arr = [''] * len(conf_arr)

	sum_sim = conf_tab['sum_sim'].value
	sum_obs = conf_tab['sum_obs'].value
	sum_fit_sim = conf_tab['sum_fit_sim'].value
	sum_fit_obs = conf_tab['sum_fit_obs'].value
	ratio_fit_arr = conf_tab['sum_fit_obs'].value/conf_tab['sum_fit_sim'].value 
	ratio_arr = conf_tab['sum_obs'].value/conf_tab['sum_sim'].value 

	for j, wide in enumerate(conf_tab['wide']):
		wide_arr[j] = float(wide.replace('mrs0',''))

	wide_arr = np.array(wide_arr)

	#rati_arr
	ax1[i].scatter(wide_arr, ratio_arr, s=25, ec='C0', fc='white', label='sum')
	ax1[i].scatter(wide_arr, ratio_fit_arr, s=25, ec='C1', fc='black', label='fit')

	ax1[i].legend(loc='lower left', fontsize=8)
	ax1[i].set_xlabel('Size')
	ax1[i].set_ylabel('Flux Density ratio (obs/sim)')

	ax1[i].grid(True, ls=':')
	ax1[i].text(0.6, 0.95, conf, transform=ax1[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')
	# ax1[i].set_ylim([rati_arr.min(), 1])

	#togeather
	ax2[i].scatter(wide_arr, sum_obs, s=25, ec='C0', fc='white', label='obs sum')
	ax2[i].scatter(wide_arr, sum_sim, s=25, ec='C0', fc='white', marker='s', label='sim sum')
	ax2[i].scatter(wide_arr, sum_fit_obs, s=25, ec='C1', fc='black', label='sim fit')
	ax2[i].scatter(wide_arr, sum_fit_sim, s=25, ec='C1', fc='black', marker='s', label='obs fit')

	ax2[i].set_xlabel('Size')
	ax2[i].set_ylabel('Flux Density [Jy]')

	ax2[i].grid(True, ls=':')
	ax2[i].text(0.6, 0.95, conf, transform=ax2[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	ax2[i].legend(loc='upper left', fontsize=8)

	# ax1[i].set_yscale('log')
	ax2[i].set_yscale('log')
	
fig1.tight_layout()
fig2.tight_layout()

fig1.savefig(f'scatter_{which}_ratio.pdf', dpi=300, bbox_inches='tight')
fig1.savefig(f'scatter_{which}_ratio.png', dpi=300, bbox_inches='tight')

fig2.savefig(f'scatter_{which}_obs_sim.pdf', dpi=300, bbox_inches='tight')
fig2.savefig(f'scatter_{which}_obs_sim.png', dpi=300, bbox_inches='tight')

for i, conf in enumerate(conf_unique):
	ax1[i].set_yscale('log')

fig1.tight_layout()
fig1.savefig(f'scatter_{which}_ratio_log.pdf', dpi=300, bbox_inches='tight')
fig1.savefig(f'scatter_{which}_ratio_log.png', dpi=300, bbox_inches='tight')

plt.close('all')

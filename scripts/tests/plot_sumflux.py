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

which = 'gaussians'
# which = 'disks' 

table = QTable.read(f'table_{which}.fits')
conf_unique = np.unique(table['conf']) 
conf_unique =  custom_sort(conf_unique)

fig1, ax1 = plt.subplots(2, 5, figsize=(15, 5))
fig2, ax2 = plt.subplots(2, 5, figsize=(15, 5))
fig3, ax3 = plt.subplots(2, 5, figsize=(15, 5))
fig4, ax4 = plt.subplots(2, 5, figsize=(15, 5))

ax1 = ax1.flatten()
ax2 = ax2.flatten()
ax3 = ax3.flatten()
ax4 = ax4.flatten()

colors = cm.turbo(np.linspace(0, 1, len(ax1)))

for i, conf in enumerate(conf_unique):

	conf_tab = table[np.where(table['conf']==conf)]
	conf_arr = conf_tab['conf']
	wide_arr = [''] * len(conf_arr)
	sum_sim = conf_tab['sum_sim'].value
	sum_obs = conf_tab['sum_obs'].value
	rati_arr = conf_tab['sum_obs/sum_sim'].value

	for j, wide in enumerate(conf_tab['wide']):
		wide_arr[j] = float(wide.replace('mrs0',''))

	wide_arr = np.array(wide_arr)

	# sum_sim
	ax1[i].scatter(wide_arr, sum_sim, ec='none', fc='grey', s=25)
	ax1[i].scatter(wide_arr, sum_sim, ec='none', fc='white', s=12)
	ax1[i].scatter(wide_arr, sum_sim, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax1[i].set_xlabel('Size')
	ax1[i].set_ylabel('Flux Density (obs) [Jy]')

	ax1[i].grid(True, ls=':')
	ax1[i].text(0.6, 0.95, conf, transform=ax1[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	# sum_obs
	ax2[i].scatter(wide_arr, sum_obs, ec='none', fc='grey', s=25)
	ax2[i].scatter(wide_arr, sum_obs, ec='none', fc='white', s=12)
	ax2[i].scatter(wide_arr, sum_obs, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax2[i].set_xlabel('Size')
	ax2[i].set_ylabel('Flux Density (sim) [Jy]')

	ax2[i].grid(True, ls=':')
	ax2[i].text(0.6, 0.95, conf, transform=ax2[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	#rati_arr
	ax3[i].scatter(wide_arr, rati_arr, ec='none', fc='grey', s=25)
	ax3[i].scatter(wide_arr, rati_arr, ec='none', fc='white', s=12)
	ax3[i].scatter(wide_arr, rati_arr, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax3[i].set_xlabel('Size')
	ax3[i].set_ylabel('Flux Density ratio (obs/sim)')

	ax3[i].grid(True, ls=':')
	ax3[i].text(0.6, 0.95, conf, transform=ax3[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	#togeather
	ax4[i].scatter(wide_arr, sum_obs, ec='none', fc='grey', s=25, label='obs')
	ax4[i].scatter(wide_arr, sum_obs, ec='none', fc='white', s=12)
	ax4[i].scatter(wide_arr, sum_obs, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax4[i].scatter(wide_arr, sum_sim, ec='none', fc='black', s=25, label='sim')
	ax4[i].scatter(wide_arr, sum_sim, ec='none', fc='white', s=12)
	ax4[i].scatter(wide_arr, sum_sim, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax4[i].set_xlabel('Size')
	ax4[i].set_ylabel('Flux Density (sim) [Jy]')

	ax4[i].grid(True, ls=':')
	ax4[i].text(0.6, 0.95, conf, transform=ax4[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	ax4[i].legend(loc='upper left')

	ax1[i].set_yscale('log')
	ax2[i].set_yscale('log')
	ax3[i].set_yscale('log')
	ax4[i].set_yscale('log')

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()

fig1.savefig(f'scatter_obs_{which}.pdf', dpi=300, bbox_inches='tight')
fig2.savefig(f'scatter_sim_{which}.pdf', dpi=300, bbox_inches='tight')
fig3.savefig(f'scatter_ratio_{which}.pdf', dpi=300, bbox_inches='tight')
fig4.savefig(f'scatter_obs_sim_{which}.pdf', dpi=300, bbox_inches='tight')

plt.close('all')

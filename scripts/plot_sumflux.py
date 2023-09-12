from glob import glob
import numpy as np
from astropy.table import QTable
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm

which = 'gaussians'
# which = 'disks'

table.read(f'table_{which}.fits')

fig1, ax1 = plt.subplots(2, 5, figsize=(15, 5))
fig2, ax2 = plt.subplots(2, 5, figsize=(15, 5))
fig3, ax3 = plt.subplots(2, 5, figsize=(15, 5))

ax1 = ax1.flatten()
ax2 = ax2.flatten()
ax3 = ax3.flatten()

colors = cm.turbo(np.linspace(0, 1, len(ax)))

for conf in table['conf']:

	conf_tab = table[np.where(table['conf']==conf)]
	conf_arr = conf_tab['conf']
	wide_arr = [''] * len(conf_arr)
	sum_sim = conf_tab['sum_sim'].value
	sum_obs = conf_tab['sum_obs'].value
	rati_arr = conf_tab['sum_obs/sum_sim'].value

	for wide in table['wide']:
		wide_arr[i] = float(wide.replace('mrs0',''))

	wide_arr = np.array(wide_arr)

	# sum_sim
	ax1[i].scatter(wide_arr, sum_sim, ec='none', fc='grey', s=25)
	ax1[i].scatter(wide_arr, sum_sim, ec='none', fc='white', s=12)
	ax1[i].scatter(wide_arr, sum_sim, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax1[i].set_xlabel('Size (?)')
	ax1[i].set_ylabel('Flux Density (obs) [Jy]')

	ax1[i].grid(True, ls=':')
	ax1[i].text(0.5, 0.95, conf, transform=ax1[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	# sum_obs
	ax2[i].scatter(wide_arr, sum_obs, ec='none', fc='grey', s=25)
	ax2[i].scatter(wide_arr, sum_obs, ec='none', fc='white', s=12)
	ax2[i].scatter(wide_arr, sum_obs, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax2[i].set_xlabel('Size (?)')
	ax2[i].set_ylabel('Flux Density (sim) [Jy]')

	ax2[i].grid(True, ls=':')
	ax2[i].text(0.5, 0.95, conf, transform=ax2[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

	#rati_arr
	ax3[i].scatter(wide_arr, rati_arr, ec='none', fc='grey', s=25)
	ax3[i].scatter(wide_arr, rati_arr, ec='none', fc='white', s=12)
	ax3[i].scatter(wide_arr, rati_arr, alpha=0.5, s=25, ec='none', fc=colors[i])

	ax3[i].set_xlabel('Size (?)')
	ax3[i].set_ylabel('Flux Density ratio (obs/sim)')

	ax3[i].grid(True, ls=':')
	ax3[i].text(0.5, 0.95, conf, transform=ax3[i].transAxes, weight='extra bold', fontsize=10, va='top', ha='center')

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

fig1.savefig(f'scatter_obs_{which}.pdf', dpi=300, bbox_inches='tight')
fig2.savefig(f'scatter_sim_{which}.pdf', dpi=300, bbox_inches='tight')
fig3.savefig(f'scatter_ratio_{which}.pdf', dpi=300, bbox_inches='tight')

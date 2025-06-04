import numpy as np
import pickle
import platform
import matplotlib.pyplot as plt
import pandas as pd 
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter, MonthLocator


if platform.system()=='Darwin':
    hd = '/Users/ameghino/'
elif platform.system()=='Linux':
    hd = '/home/calving/'

with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_2.pkl', "rb") as file:
    t2 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_4.pkl', "rb") as file:
    t1 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_-4.pkl', "rb") as file:
    t_1 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_-2.pkl', "rb") as file:
    t_2 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_0.pkl', "rb") as file:
    t0 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/temp_seisitivity_obs.pkl', "rb") as file:
    tobs = pickle.load(file)

# plotting setting
plt.rcParams.update({'font.size': 7})
inch2cm = 1/2.54  # centimeters in inches
plt.close('all')
fig, axs = plt.subplots(1, 1, figsize=(12*inch2cm, 15*inch2cm), num=1)
axs.plot(t_1['T'].T[50], label='-4')
axs.plot(t_2['T'].T[50], label='-2')
axs.plot(t0['T'].T[50], label='0')
axs.plot(t2['T'].T[50], label='+2')
axs.plot(t1['T'].T[50], label='+4')
axs.plot(tobs['T'].T[50], label='obs')

axs.legend()
axs.set_ylabel('T at I/S boundary')

fig, axs = plt.subplots(2, 1, figsize=(6*inch2cm, 10*inch2cm), num=2)
t2_2023 = t2['SI'][t2['SI'].index.year >= 2023]
t1_2023 = t1['SI'][t1['SI'].index.year >= 2023]

axs[0].plot(t_1['SI'].index , t_1['SI'].values - np.nanmin(t_1['SI'].values), '--',label='-4$^{\circ}$C', linewidth=0.8)
axs[0].plot(t_2['SI'].index , t_2['SI'].values - np.nanmin(t_2['SI'].values), '--', label='-2$^{\circ}$C', linewidth=0.8)
axs[0].plot(t0['SI'].index , t0['SI'].values - np.nanmin(t0['SI'].values),'-k', label='0$^{\circ}$C', linewidth=0.8)
axs[0].plot(t2_2023.index , t2_2023.values - np.nanmin(t2_2023.values), '-.', label='+2$^{\circ}$C', linewidth=0.8)
axs[0].plot(t1_2023.index , t1_2023.values - np.nanmin(t1_2023.values), '-.', label='+4$^{\circ}$C', linewidth=0.8)
# axs[0].plot(tobs['SI'].index , tobs['SI'].values - np.nanmin(tobs['SI'].values), '-k', label='obs', linewidth=0.8)
axs[0].legend(handlelength=1, loc='upper left',ncol=2)
axs[0].set_ylabel('SI [cm]', fontweight='bold')
axs[0].set_xlim(pd.to_datetime('2023-06-01'), pd.to_datetime('2023-09-01'))
axs[0].set_ylim(0,15)
# Set major ticks to the 1st of each month
axs[0].xaxis.set_major_locator(MonthLocator(bymonthday=1))  # Every 1st of the month
# Set the formatter to show short month name, e.g., Jan, Feb
axs[0].xaxis.set_major_formatter(DateFormatter("%b-%d"))  # Or "%Y-%m-%d" for full format
axs[0].set_xlabel('Date in 2023', fontweight='bold')


axs[1].plot((t_1['CC']/t_1['CC'][0])*1e2, '--', label='-4',  linewidth=0.8)
axs[1].plot((t_2['CC']/t_2['CC'][0])*1e2, '--', label='-2',  linewidth=0.8)
axs[1].plot((t0['CC']/t0['CC'][0])*1e2,'-k', label='0', linewidth=0.8)
axs[1].plot((t2['CC']/t2['CC'][0])*1e2, '-.', label='+2',  linewidth=0.8)
axs[1].plot((t1['CC']/t1['CC'][0])*1e2, '-.', label='+4',  linewidth=0.8)
# axs[1].plot((tobs['CC']/tobs['CC'][0])*1e2, '-k', label='obs',  linewidth=0.8)
# axs[1].legend()
axs[1].set_ylabel('CC [%]', fontweight='bold')
axs[1].set_xlabel('Days since Sep. 9, 2022', fontweight='bold')
axs[1].set_ylim(95, 200)

panel_labels = ["a", "b", "c", "d", "e"]
for ax, label in zip(axs, panel_labels):    
    # Add panel label outside the top-left corner
    ax.text(-0.1, 1.2, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left')
    ax.grid(linestyle='--', linewidth=0.5, zorder=0)

plt.tight_layout()
fig.align_ylabels()

plt.savefig(hd + '/Nextcloud/QaanaaqIceCap/refreezing/sensitivity_T.png',dpi=300)

plt.show()

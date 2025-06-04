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

with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_1.pkl', "rb") as file:
    snow01 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_0.6.pkl', "rb") as file:
    snow06 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_0.8.pkl', "rb") as file:
    snow08 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_1.2.pkl', "rb") as file:
    snow12 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_1.4.pkl', "rb") as file:
    snow14 = pickle.load(file)
with open(hd + '/Nextcloud/QaanaaqIceCap/refreezing/snow_seisitivity_obs.pkl', "rb") as file:
    snowobs = pickle.load(file)


# plotting setting
plt.rcParams.update({'font.size': 7})
inch2cm = 1/2.54  # centimeters in inches
plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(12*inch2cm, 15*inch2cm), num=3)
axs[0].plot(snow01['T'].T[50], label='0.0m')
axs[0].plot(snow06['T'].T[50], label='0.4m')
axs[0].plot(snow08['T'].T[50], label='0.8m')
axs[0].plot(snow12['T'].T[50], label='1.2m')
axs[0].plot(snow14['T'].T[50], label='1.6m')
axs[0].legend()
axs[0].set_ylabel('T at I/S boundary')

fig, axs = plt.subplots(2, 1, figsize=(6*inch2cm, 10*inch2cm), num=2)

axs[0].plot(snow06['SI'].index , snow06['SI'].values - snow06['SI'].values[1],'--', label='-40%', linewidth=0.8)
axs[0].plot(snow08['SI'].index , snow08['SI'].values - snow08['SI'].values[1],'--', label='-20%', linewidth=0.8)
axs[0].plot(snow01['SI'].index , snow01['SI'].values - snow01['SI'].values[1],'-k', label='0%', linewidth=0.8)
axs[0].plot(snow12['SI'].index , snow12['SI'].values - snow12['SI'].values[1],'-.', label='+20%', linewidth=0.8)
axs[0].plot(snow14['SI'].index , snow14['SI'].values - snow14['SI'].values[1],'-.', label='+40%', linewidth=0.8)
# axs[0].plot(snowobs['SI'].index , snowobs['SI'].values - snowobs['SI'].values[1],'-k', label='obs', linewidth=0.8)

axs[0].legend(handlelength=1, loc='upper left',ncol=2)

axs[0].set_ylabel('SI [cm]', fontweight='bold')
axs[0].set_xlim(pd.to_datetime('2023-06-01'), pd.to_datetime('2023-09-01'))
axs[0].set_ylim(0,15)
axs[0].xaxis.set_major_locator(MonthLocator(bymonthday=1))  # Every 1st of the month
# Set the formatter to show short month name, e.g., Jan, Feb
axs[0].xaxis.set_major_formatter(DateFormatter("%b-%d"))  # Or "%Y-%m-%d" for full format
axs[0].set_xlabel('Date in 2023', fontweight='bold')

axs[1].plot((snow06['CC']/snow06['CC'][0])*1e2, '--',label='0.4 m', linewidth=0.8)
axs[1].plot((snow08['CC']/snow08['CC'][0])*1e2, '--', label='0.8 m', linewidth=0.8)
axs[1].plot((snow01['CC']/snow01['CC'][0])*1e2, '-k', label='0.0 m', linewidth=0.8)
axs[1].plot((snow12['CC']/snow12['CC'][0])*1e2, '-.', label='1.2 m', linewidth=0.8)
axs[1].plot((snow14['CC']/snow14['CC'][0])*1e2, '-.', label='1.6 m', linewidth=0.8)
# axs[1].plot((snowobs['CC']/snowobs['CC'][0])*1e2, '-k', label='obs', linewidth=0.8)

axs[1].set_ylabel('CC [%]', fontweight='bold')
axs[1].set_xlabel('Days since Sep. 9, 2022', fontweight='bold')
axs[1].set_ylim(95, 200)

panel_labels = ["c", "d", "e"]
for ax, label in zip(axs, panel_labels):    
    # Add panel label outside the top-left corner
    ax.text(-0.1, 1.2, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left')
    ax.grid(linestyle='--', linewidth=0.5, zorder=0)

plt.tight_layout()
fig.align_ylabels()

plt.savefig(hd + '/Nextcloud/QaanaaqIceCap/refreezing/sensitivity_S.png',dpi=300)

plt.show()

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pandas as pd
from matplotlib import cm
import platform
from matplotlib.colors import Normalize
# from tqdm import tqdm  # For progress bar
import matplotlib.animation as animation
import pandas as pd 
import matplotlib.dates as mdates
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

if platform.system()=='Darwin':
    hd = '/Users/ameghino/'
elif platform.system()=='Linux':
    hd = '/home/calving/'


#TODO temperature should also given at snow and ice layer?
# add snow layer, which changes varies with SIGMA-B accumulation

def implicit_step(T, kappa, dz, dt):
    """
    Perform one implicit (backward Euler) step for 1D heat diffusion equation.

    Parameters:
    T     - Current temperature profile (1D array)
    kappa - Thermal diffusivity (array, same size as T)
    dz    - Grid spacing
    dt    - Time step

    Returns:
    T_new - Updated temperature profile after one time step
    """
    nz = len(T)
    alpha = kappa * dt / dz**2

    # Build the tridiagonal matrix A
    main_diag = 1 + 2 * alpha
    off_diag = -alpha
    diagonals = [off_diag * np.ones(nz-1), main_diag * np.ones(nz), off_diag * np.ones(nz-1)]
    A = diags(diagonals, offsets=[-1, 0, 1], format='csc')

    # Apply boundary conditions
    A = A.toarray()  # Convert to dense for boundary manipulation
    A[0, :] = 0
    A[0, 0] = 1  # Fix lower boundary
    A[-1, :] = 0
    A[-1, -1] = 1  # Fix upper boundary

    # Right-hand side (previous temperature profile)
    T_rhs = T.copy()

    # Enforce boundary conditions on RHS
    T_rhs[0] = T[0]   # Fixed lower boundary temperature
    T_rhs[-1] = T[-1] # Fixed upper boundary temperature

    # Solve the linear system A T_new = T_rhs
    T_new = spsolve(A, T_rhs)

    return T_new


# ------ replace equation for k
def k_s(den, temp): 

    rho = den # QC stake-6
    rho_transition = 450.0
    a = 0.02 #[m^3/kg]
    RHO_I       = 917.0 
    theta       = 1 / (1 + np.exp(-2*a*(rho - rho_transition)))
    kref_firn   = 2.107 + 0.003618 * (rho - RHO_I)
    kref_snow   = 0.024 - 1.23e-4 * rho + 2.5e-6 * rho**2
    kref_i      = 2.107 # [W/m/K]
    kref_a      = 0.024 # [W/m/K] 
    K_air       = kref_a # use this for now; at some point find equation for T-dependence of air
    phi_0 = 273.15+temp # temperature
    K_ice1 = 9.828 * np.exp(-0.0057 * phi_0)
    K_ice1K_air = 2.107 * 0.024     # TODO should be temperature dependent.  ki(T) and ka(T) are the ice and air conductivity at the temperature T
    kref_ikref_a = 2.107 * 0.024    #  [W/m/K] the ice and air conductivities at the reference temperature Tref of −3°C
    K_firn = (1-theta) * K_ice1K_air/(kref_ikref_a) * kref_snow + theta * K_ice1/kref_i * kref_firn # equation 

    return K_firn # [W/m/K]



# geometry
nz = 18
z_min = -9.3 # bottom
z_max = -1.3 # top
Z = np.linspace(z_min,z_max, nz)
dz = Z[1]-Z[0]

# Material parameter
# Calonne et al., 2019 GRL eq (5)
rho_i = 900                 # ice density kg/m^3
cp_i = 2097                 # heat capacity J/kg/K
k_i = 2.1                   # diffusivity W/m/K = J/s/m/K
kappa_i = k_i/(rho_i*cp_i)  # [s/]
# 
rho_s = 300                 # snow/firn density kg/m^3
cp_s = 2090                 # heat capacity J/kg/K
#  k_s = 0.4                # diffusivity W/m/K = J/s/m/K


Lf_i = 333.6 * rho_i # [kJ/m3] latent heat of fusion ice : [kJ/kg] * [kg/m3] 
Lf_s = 333.6 * 500 # [kJ/m3] latent heat of fusion snow : [kJ/kg] * [kg/m3] 

kappa = kappa_i * np.ones(np.size(Z))
# kappa[np.where(Z>-1.3)] = kappa_s

# daily temperature was calculated by ***
df = pd.read_csv(hd + '/Nextcloud/QaanaaqIceCap/thermocouple/daily_mean_temperatures.csv',
                     skiprows=0, index_col='Date')
t_mat = np.fliplr(df)
d = -np.array([130,230,330,530,930])/100
d = np.flipud(d)

#TODO set start day as Jan 1, 2023
ini_date = '2023-06-24'
end_date = '2023-08-28'

dmy = df.loc[ini_date]
dmy2 = df.loc[end_date]

initial_t = np.flipud(dmy[3:].tolist())
end_t = np.flipud(dmy2[3:].tolist())
end_t[-1] = 0
T = np.interp(Z, d, initial_t, left=None, right=None, period=None)
Te = np.interp(Z, d, end_t, left=None, right=None, period=None)

ini_date = pd.to_datetime(ini_date)
end_date = pd.to_datetime(end_date)


# Time loop
#TODO change it to daily resolution
# dt = (dz**2/kappa)/2  # [s]
dt = 1 * 24 * 60 * 60
no_tstep = 61 # len(t_mat[:,0])
t_res = []
#cm = jet(no_tstep)
color = iter(cm.brg(np.linspace(0, 1, no_tstep)))
cb_time = []
# Create colormap and normalization
cmap = cm.brg
norm = Normalize(vmin=0, vmax=no_tstep)

plt.rcParams.update({'font.size': 7})
inch2cm = 1/2.54  # centimeters in inches
plt.close('all')
fig = plt.figure(1, figsize=(10*inch2cm, 6*inch2cm))
gs = fig.add_gridspec(1, 1)
ax1 = fig.add_subplot(gs[0])
# ax2 = fig.add_subplot(gs[1])
ax1.plot([-15, 1], [-1.3, -1.3], '--k')
# ax2.plot([-14, 1], [-1.3, -1.3], '--k')
ax1.grid()
# ax2.grid()
ax1.set_xlim(-8.5,0.5)
# ax2.set_xlim(-13.5,0.5)
ax1.set_ylim(-10,0)
# ax2.set_ylim(-9.5,0)
# ax2.set_title('(b) Modelled', fontweight="bold", loc='left')
# ax2.set_xlabel('T [$^{\circ}$C]')
ax1.set_ylabel('Depth [m]')
# ax2.set_ylabel('Depth [m]')
frames = []
# Predefine the lines before the loop
# mod_line, = ax1.plot([], [], 'x-b', linewidth=1, label='mod')  # Modeled temperature
ax1.plot(T, Z, 'x-', label='initial')
ax1.plot(Te, Z, 'x-', label='end')

plt.legend()
T_mat = []
Tobs_mat = []
date = []
new_date = ini_date
step = 1
Ind = np.arange(1,nz-1) # indexes without boundary
for tstep in range(0, no_tstep, step):
    T_old = T
    # print(f'{tstep} / {no_tstep}')
    # Apply boundary conditions
    # T[0] = t_mat[tstep, 0]   # Lower boundary

    # calculate kappa for firn for given temperature and density
    id = np.where(Z > -1.3)
    k_firn = k_s(rho_s, T[id])
    kappa[id] = k_firn/(rho_s*cp_s)  # [s/]
    ck = kappa * dt / (dz * dz)
    # Solve the equation (explicit method)
    # T[Ind] = T[Ind] + dt * kappa[Ind] / dz**2 * (T[Ind + 1] - 2 * T[Ind] + T[Ind - 1])
    T[Ind] = 1.0 / (1.0 + ck[Ind]) * ((1.0 - ck[Ind]) * T_old[Ind] \
    + 0.5 * ck[Ind] * T_old[Ind+1] + 0.5 * ck[Ind] * T_old[Ind-1] \
    + 0.5 * ck[Ind] * T[Ind+1] + 0.5 * ck[Ind] * T[Ind-1])
    # Update modeled profile
    # mod_line.set_data(T, Z)
    T_mat.append(T.tolist())
    dmy = df.loc[new_date.strftime('%Y-%m-%d')]
    tdmy = np.flipud(dmy[3:].tolist())
    tdmy[-1]=0
    tobs = np.interp(Z, d, tdmy, left=None, right=None, period=None)
    Tobs_mat.append(tobs.tolist())

    # Upper boundary insert air temperature
    new_date = new_date + pd.Timedelta(days=step)
    date.append(new_date)
    tdmy = df.loc[new_date.strftime('%Y-%m-%d')]
    T[-1] = 0 # add 0 degC
    T[0] = tdmy[-1] # set lower boundary
    color = cmap(norm(tstep))
    ax1.plot(T, Z, '-',color=color,linewidth=0.3)


    # plt.title(df.index[tstep])
    # return mod_line
# Create FuncAnimation
# anm = animation.FuncAnimation(fig, update, frames=no_tstep, interval=50, blit=True)
# gif_filename = hd + "/Nextcloud/QaanaaqIceCap/thermocouple/heat_conduction.gif"
# anm.save(gif_filename, writer='pillow', fps=10)


plt.tight_layout()

# 
plt.savefig(hd + '/Nextcloud/QaanaaqIceCap/refreezing/check_refreezing.png',dpi=300)

### 
# modelled refreezing
T_mat = np.array(T_mat)
dH = np.cumsum( np.sum(dz * (np.diff(T_mat, axis=0)) * rho_i * cp_i, axis=1)) * 1e-6 # (MJ/m2)
RF = (dH /333.6)*100 # cm 

# observed refreezing
Tobs_mat = np.array(Tobs_mat)
dHobs = np.cumsum( np.sum(dz * (np.diff(Tobs_mat, axis=0)) * rho_i * cp_i, axis=1)) * 1e-6 # (MJ/m2)
RFobs = (dHobs /333.6)*100 # cm 


fig,ax = plt.subplots(1,1,num=2, figsize=(8*inch2cm, 6*inch2cm))

ax.plot(date[1:],RF,'--',label='Model',color='darkorange')
ax.plot(date[1:],RFobs,label='Observation',color='darkorange')
ax.legend(loc='center right')
ax.set_ylabel('Superimposed ice [cm]',color='darkorange', fontweight='bold')
ax.set_xlabel('Date in 2023', fontweight='bold')
ax.set_xlim(date[1], date[-1])
ax.grid(linestyle='--',linewidth=0.5)
ax.set_ylim(0,60)

# Format the date labels
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))  # Format as 'Month-Day' (e.g., 'Jul-14')
plt.gcf().autofmt_xdate()  # Rotate dates for better readability

#show height calibrated by Q6 observation
# melt of snow [cm i.e.]
data = {
    "Timestamp": ["2023/07/14 13:59:00", "2023/07/30 15:53:00", "2023/08/04 12:48:00"],
    "Value": [68, 4.5, 19]
}
# Convert to DataFrame
df = pd.DataFrame(data)
# Convert Timestamp column to datetime format
df["Timestamp"] = pd.to_datetime(df["Timestamp"])
# load snow height
path = hd + '/Nextcloud/QaanaaqIceCap/AWS/A20241213-002/DATA/SIGMA_AWS_SiteB_2020-2024_Lv1_1.csv'
data = pd.read_csv(
    path, encoding='utf-8')
data['Date'] = pd.to_datetime(data['Date'])
data['sh'] = data['sh'].replace(-9999, np.nan)
data['sh'] = data['sh'].replace(300, np.nan)
data['sh'] = data['sh'].astype(str).astype(float)
ax2 = ax.twinx()
snow2ice = 450/900
closest_times = pd.merge_asof(
    df,
    data,
    left_on='Timestamp',
    right_on='Date',
    direction='nearest'
)
offset=36.3
ax2.plot(data['Date'], (data['sh']-offset)*snow2ice,label = 'AWS',color='green')
ax2.plot(df['Timestamp'], df['Value']*snow2ice,'o',label='Observation',markeredgecolor='green',markerfacecolor='white')
ax2.set_ylim(0,60)
ax2.set_ylabel('Snow depth [cm i.eq.]',color='green', fontweight='bold')
ax2.legend()

plt.tight_layout()
plt.savefig(hd + '/Nextcloud/QaanaaqIceCap/refreezing/refreezing_ts.png',dpi=300)



plt.show()
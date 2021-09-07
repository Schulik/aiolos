import numpy as np 
import matplotlib.pyplot as plt

from load_aiolos import load_aiolos_snap, load_aiolos_diag 

R = 8.31446261815324e7
sigma_SB = 5.670374419e-5


mu_l = 169  # Olivine weight
mu_c = 30   # Gas weight

Cc = 1/(1.3-1) * R / mu_c
Cl = 1.5 * R / mu_c

L = 3.21e10

def Pvap(T):
    return 6.72e14*np.exp(- L*mu_l / (R*T))

def rho_vap(T):
    return Pvap(T) * mu_c / (R*T)

def make_plot(sim_dir, title=None, snap=-1):

    def snapname(species):
        return "{}/output_dusty_wind_{}_t{}.dat".format(sim_dir, species, snap)

    def add_crit_points(ax):
        ax.axvline(r_s, c='k', ls='-')
        ax.axvline(r_H, c='k', ls='--')

    gas  = load_aiolos_snap(snapname("Gas")) 
    dust = load_aiolos_snap(snapname("Grain0.1")) 
    diag = load_aiolos_diag("{}/diagnostic_dusty_wind_t{}.dat".format(sim_dir, snap))


    # Get the Hill radius
    r_H = gas['x'][np.argmax(gas['potential'][1:-1])+1] 

    idx = gas['x'].searchsorted(2*r_H-gas['x'][0])
    gas  = gas[:idx]
    dust = dust[:idx]
    diag = diag[:idx]
    


    gamma_r = diag['kappa_PTsun0_1']/diag['kappa_P0_1']
    print("gamma_r:", gamma_r[0], gamma_r[-1])
    T_thin = (0.0625*gamma_r*diag['Stot']/sigma_SB)**0.25


    # Compute the location of the sonic point:
    Mach = gas['velocity']/gas['soundspeed']
    r_s = gas['x'][np.abs(Mach[1:]-1).argmin()+1] 

    f, ax = plt.subplots(3,1, sharex=True, figsize=(6, 8))

    T_ad = gas['temperature'][-3] * (gas['density']/gas['density'][-3])**0.285714

    cg = ax[2].plot(gas['x'], gas['temperature'],  label='Gas')[0].get_color()
    cd = ax[2].plot(dust['x'], dust['temperature'],  label='Dust')[0].get_color()
    ax[2].plot(dust['x'], T_ad, 'k--',  label='Adiabatic')
    ax[2].plot(diag['x'], T_thin, 'k:', label='T_irr (dust)')

    ax[2].set_xlabel('Radius [cm]')
    ax[2].set_ylabel('Temperature [K]')
    ax[2].legend(loc='best')
    ax[2].set_ylim(0.95*gas['temperature'][1:-1].min(), 
                   1.05*dust['temperature'][1:-1].max())
   


    add_crit_points(ax[2])


    rho_v = rho_vap(dust['temperature']) 
    rho_v /= np.sqrt(gas['temperature']/dust['temperature'])

    rho_v2 = rho_vap(gas['temperature']) 

    ax[0].plot(gas['x'], gas['density'], c=cg, label='Gas')
    ax[0].plot(gas['x'], rho_v, 'k--', label='Vapor (T_dust)')
    ax[0].plot(gas['x'], rho_v2, 'k:', label='Vapor (T_gas)')

   
    ax[0].set_ylabel('Gas Density [g/cm^3]')
    ax[0].set_ylim(ymin=1e-4*gas['density'].max())
    ax[0].legend(loc='best')
    ax[0].set_yscale('log')

    add_crit_points(ax[0])

    ax[1].plot(dust['x'], dust['density']/gas['density'], label='Dust')
   
    ax[1].set_ylabel('Dust-to-gas ratio')

    ax[1].set_yscale('log')


    add_crit_points(ax[1])


    ax[0].set_xlim(gas['x'].min())

    if title is not None:
        f.suptitle(title)
        
    f.subplots_adjust(
        top=0.935, bottom=0.07,
        left=0.135, right=0.955,
        hspace=0.12)

    mdot_g = np.pi*(gas['x']**2 * gas['momentum'])[-1]   
    mdot_d = np.pi*(dust['x']**2 * dust['momentum'])[-1]

    print('Mass-loss Rates: Gas={:.3g}, Dust={:.3g}'.format(mdot_g, mdot_d))
    


if __name__ == "__main__":
    make_plot('./', snap=-1)

    plt.show()

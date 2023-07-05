import os
import numpy as np

from load_aiolos import (
    load_aiolos_snap,
    load_aiolos_diag,
    load_aiolos_params
)
sigma_rad = 5.67e-5

def Guillot_2band(tau, Tint4, Tirr40, gamma0, Tirr41, gamma1,
                  fH=1, mu_star=1):
    
    term_int = 0.75 * tau + 1./(4*fH)

    x = mu_star/gamma0
    term0 = 0.75*mu_star*(1./(3*fH) + x + (1./(3*x) - x)*np.exp(-tau/x))

    x = mu_star/gamma1
    term1 = 0.75*mu_star*(1./(3*fH) + x + (1./(3*x) - x)*np.exp(-tau/x))


    return (Tint4*term_int + Tirr40*term0 + Tirr41*term1)**0.25


def check_structure(sim, L1_target, make_plots):
    filename = 'diagnostic_' + sim + '_t-1.dat'
    
    data = load_aiolos_diag(filename, species=1)

    Tirr40 = 0.25*data['S0'][-1]/sigma_rad
    gamma0 = (data['kappa_PTsun0_0'] / data['kappa_P0_0'])[-1]

    if 'S1' in data.dtype.names:
        gamma1 = (data['kappa_PTsun1_0'] / data['kappa_P0_0'])[-1]
        Tirr41 = 0.25*data['S1'][-1]/sigma_rad
    else:
        gamma1 = 1
        Tirr41 = 0

    tau = data['tau_radial0']/gamma0

    # Get the planet's internal temperature
    pars = load_aiolos_params(f"{sim}.par")
    Tint = 0
    if 'USE_PLANET_TEMPERATURE' in pars:
        if 'PARI_TINT' in pars:
            Tint = float(pars['PARI_TINT'])
        else: 
            Tint = float(pars.get('PARI_TPLANET',0))

    # Get the flux-limiter factor
    fH = 1/float(pars.get('XI_RAD',1))

    T_guillot = Guillot_2band(tau, Tint**4, Tirr40, gamma0, Tirr41, gamma1,
                              fH=fH)

    idx = (tau < 100) & (tau > 1e-10)
    L1 = np.abs(1 - data['T_gas0'][idx]/T_guillot[idx]).mean()

    if L1 <= L1_target:
        print(f'Irradiation test {sim[-1]} L1 check passed')
    else:
        print(f"Irradiation test {sim[-1]} L1 check failed. " + 
              "L1={}, expected={}".format(L1,L1_target))
        


    if make_plots:
        import matplotlib.pyplot as plt
             
        plt.semilogy(data['T_gas0'], tau, label=f'{sim}')
        plt.loglog(T_guillot, tau, 'k--', label='Guillot (2010), $f_H=1$')

        ylim = plt.ylim()
        plt.ylim(ylim[1], 1e-10)
        plt.xlim(200, 2e4)

        plt.xlabel(r'Temperature [K]')
        plt.ylabel(r'Optical Depth')
        
        plt.legend(loc='best')

        os.makedirs('plots', exist_ok=True)
        plt.savefig(f'plots/{sim}.png')
        plt.clf()
    
    np.testing.assert_array_less(L1, L1_target)


def test_structure(make_plots=False):
    par_files = [
        ("irradiation1", 0.02),
        ("irradiation2", 0.06)
        ]
    
    for f, tol in par_files:
        check_structure(f, tol, make_plots)
    #check_structure('irradiation', 175.0, 1.5e-3, make_plots)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser("Test structure of irradiated atmosphere")
    parser.add_argument("-p", "--make_plots", default=None,
                        action="store_true",
                        help="Make plots of the results")
    args = parser.parse_args()

    try:
        test_structure(args.make_plots)
    except AssertionError:
        pass

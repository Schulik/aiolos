import os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_diag

sigma_rad = 5.67e-5

def Guillot_2band(tau, Tint4, Tirr40, gamma0, Tirr41, gamma1,
                  fH=1, mu_star=1):
    
    term_int = 0.75 * tau + 1./(4*fH)

    x = mu_star/gamma0
    term0 = 0.75*mu_star*(1./(3*fH) + x + (1./(3*x) - x)*np.exp(-tau/x))

    x = mu_star/gamma1
    term1 = 0.75*mu_star*(1./(3*fH) + x + (1./(3*x) - x)*np.exp(-tau/x))


    return (Tint4*term_int + Tirr40*term0 + Tirr41*term1)**0.25


def test_structure(sim, Tint, L1_target, make_plots):
    filename = 'diagnostic_' + sim + '_t5.dat'
    
    data = load_aiolos_diag(filename)


    gamma0 = (data['kappa_PTsun0'] / data['kappa_P0'])[-1]
    gamma1 = (data['kappa_PTsun1'] / data['kappa_P1'])[-1]


    Tirr40 = 0.25*data['S0'][-1]/sigma_rad
    Tirr41 = 0.25*data['S1'][-1]/sigma_rad

    tau = data['tau_radial0']
    T_guillot = Guillot_2band(tau, Tint**4, Tirr40, gamma0, Tirr41, gamma1)

    idx = (tau < 100) & (tau > 1e-10)
    L1 = np.abs(1 - data['T_gas'][idx]/T_guillot[idx]).mean()

    if L1 <= L1_target:
        print('Irradiation test L1 check passed')
    else:
        print("Irradiation test L1 check failed. " + 
              "L1={}, expected={}".format(L1,L1_target))

    if make_plots:
        import matplotlib.pyplot as plt
             
        plt.semilogy(data['T_gas'], tau, label='sim')
        plt.loglog(T_guillot, tau, 'k--', label='Guillot (2010), $f_H=1$')

        ylim = plt.ylim()
        plt.ylim(ylim[1], 1e-10)
        plt.xlim(200, 2e4)

        plt.xlabel(r'Temperature [K]')
        plt.ylabel(r'Optical Depth')
        
        plt.legend(loc='best')

        os.makedirs('plots', exist_ok=True)
        plt.savefig('plots/Irradiation.png')

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser("Test structure of irradiated atmosphere")
    parser.add_argument("-p", "--make_plots", default=None,
                        action="store_true",
                        help="Make plots of the results")
    args = parser.parse_args()

    test_structure('irradiation', 175.0, 1.5e-3, args.make_plots)


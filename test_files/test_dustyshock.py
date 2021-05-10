#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params
from test_shock_tube import Riemann


def setup_riemann_solver(param_file):
    params = load_aiolos_params(param_file)
    
    UL = list(map(float, [params['PARI_INIT_DATA_U1L'],
                          params['PARI_INIT_DATA_U2L'],
                          params['PARI_INIT_DATA_U3L']]))
    UR = list(map(float, [params['PARI_INIT_DATA_U1R'],
                          params['PARI_INIT_DATA_U2R'],
                          params['PARI_INIT_DATA_U3R']]))

    GAMMA= float(params['PARI_GAMMA'])

    L = { 'rho' : 2*UL[0], 'u': UL[1], 'P' : UL[2] }
    R = { 'rho' : 2*UR[0], 'u': UR[1], 'P' : UR[2] }
                                              
    return Riemann(L, R, GAMMA)    

def plot_dusty_shock():
    param_file = 'dusty_shock.par'
    gas_file = 'output_dusty_shock_H2_t-1.dat'
    dust_file = 'output_dusty_shock_Dust_t-1.dat'

    #Get the Riemann problem
    exact = setup_riemann_solver(param_file)

    gas = load_aiolos_snap(gas_file)
    dust = load_aiolos_snap(dust_file)

    params = load_aiolos_params(param_file)
    x0 = float(params['PARI_INIT_SHOCK_MID'])
    t  = float(params['PARI_TIME_TMAX'])

    x = gas['x']
    xa = x - x0
   
    import matplotlib.pyplot as plt

    f, axes = plt.subplots(2,2)
    f.suptitle('Dusty shock')

    axes[0][0].plot(x, gas['density'], '+', c='k', label='gas')
    axes[0][0].plot(x, dust['density'], 'x',label='dust')
    axes[0][0].plot(x, exact.density(xa, t)/2, c='k', label='perfect coupling')
    
    axes[0][0].set_ylabel('density')   
    axes[0][0].xaxis.set_ticklabels([])
    axes[0][0].legend(loc='best', fontsize=8)


    axes[0][1].plot(x, gas['velocity'], '+', c='k')
    axes[0][1].plot(x, dust['velocity'], 'x')
    axes[0][1].plot(x, exact.velocity(xa, t), c='k')
    axes[0][1].set_ylabel('velocity')
    axes[0][1].xaxis.set_ticklabels([])

    axes[1][0].plot(x, gas['pressure'], '+', c='k')
    axes[1][0].plot(x, exact.pressure(xa, t), c='k')
    axes[1][0].set_ylabel('pressure')
    axes[1][0].set_xlabel('x')

    u_g= (gas ['energy'] - 0.5*gas ['momentum']*gas ['velocity'])/gas ['density']
    u_d= (dust['energy'] - 0.5*dust['momentum']*dust['velocity'])/dust['density']
    axes[1][1].plot(x, u_g, '+', c='k')
    axes[1][1].plot(x, u_d, 'x')
    axes[1][1].plot(x, exact.internal_energy(xa, t)*2, c='k')
    axes[1][1].set_ylabel('internal energy')
    axes[1][1].set_xlabel('x')

    f.subplots_adjust(top=0.9, bottom=0.11,
                      left=0.095, right=0.965,
                      hspace=0.05, wspace=0.235)

    f.savefig(os.path.join('plots/dusty_shock.png'))

    return plt


def check_dusty_shock(L1s=None):
    param_file = 'dusty_shock.par'
    gas_file = 'output_dusty_shock_H2_t-1.dat'
    dust_file = 'output_dusty_shock_Dust_t-1.dat'

    #Get the Riemann problem
    exact = setup_riemann_solver(param_file)

    gas = load_aiolos_snap(gas_file)
    dust = load_aiolos_snap(dust_file)

    params = load_aiolos_params(param_file)
    x0 = float(params['PARI_INIT_SHOCK_MID'])
    t  = float(params['PARI_TIME_TMAX'])

    x = gas['x']
    
    L1 = [
        np.abs(gas ['density']  - exact.density(x-x0, t)).mean(),
        np.abs(dust['density']  - exact.density(x-x0, t)).mean(),
        ]

    name = "dusty_shock"
    if L1s is not None:
        passed = True
        for l1, target in zip(L1, L1s):
            passed &= l1 <= target

        if passed:
            print('Test {} L1 check passed'.format(name))
        else:
            print('Test {} L1 checked failed:'.format(name))
            print('\tL1={}, target={}'.format(L1, L1s))
    else:
         print('Test {} L1 values:'.format(name))
         print('\tL1={}'.format(L1))
    
if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser("Check shock tube test results")
    parser.add_argument("-p", "--make_plots", default=None,
                        action="store_true",
                        help="Make plots of the results")
    args = parser.parse_args()

    check_dusty_shock([0.541, 0.540])
    
    if args.make_plots:
        f = plot_dusty_shock()
        #f.show()

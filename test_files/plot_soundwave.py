#!/usr/bin/python3
import sys, os
import numpy as np
import matplotlib.pyplot as plt

from load_aiolos import load_aiolos_snap, load_aiolos_params


def get_data(DIR):
    Ns = [32, 64, 128, 256, 512]
    L1s = []
    for N in Ns:
         IC_file = DIR + 'soundwave_{}_H2'.format(N) + '_t0.dat'
         IC = load_aiolos_snap(IC_file)

         try:
             snap_file = DIR + 'soundwave_{}_H2'.format(N) + '_t-1.dat'
             final = load_aiolos_snap(snap_file)
         except OSError:
             snap_file = DIR + 'soundwave_{}_H2'.format(N) + '_t-666.dat'
             final = load_aiolos_snap(snap_file)
             
         L1 = [
             np.abs(final['density'] - IC['density']).mean(),
             np.abs(final['velocity'] - IC['velocity']).mean(),        
             np.abs(final['pressure'] - IC['pressure']).mean(),
         ]
         L1s.append(L1)
         c = plt.plot(final['x'], final['velocity'], '+', label=str(N))
         plt.plot(IC['x'], IC['velocity'], 'x', c=c[0].get_color())
    plt.xlabel('x')
    plt.ylabel('v')        
         
    return np.array(Ns), np.array(L1s)


if __name__ == "__main__":
    import sys
    try:
        DIR = sys.argv[1] + '/'
    except IndexError:
        DIR = './'
        
    Ns, L1s = get_data(DIR)

    plt.figure()
    for L1, l in zip(L1s.T, [r'$\rho$', '$v$', '$p$']):
        plt.loglog(Ns, L1, label=l)
    plt.plot(Ns, 4e-4/Ns**2, 'k--')
    plt.xlabel('N')
    plt.ylabel('L1 error')
    plt.legend()
    plt.show()

#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params


def check_steady_state(problem, L1_target=None):

    IC_file = problem + '_t0.dat'
    snap_file = problem + '_t-1.dat'

    IC = load_aiolos_snap(IC_file)
    final = load_aiolos_snap(snap_file)

    L1 = [
        np.abs(final['density'] - IC['density']).mean(),
        np.abs(final['momentum'] - IC['momentum']).mean(),        
        np.abs(final['pressure'] - IC['pressure']).mean(),
    ]

    if L1_target is not None:
        passed = True
        for l1, target in zip(L1, L1_target):
            passed &= l1 <= target

        if passed:
            print('Test {} L1 check passed'.format(problem))
        else:
            print('Test {} L1 checked failed:'.format(problem))
            print('\tL1={}, target={}'.format(L1, L1_target))
    else:
        print('Test {} L1 values:'.format(problem))
        print('\tL1={}'.format(L1))


if __name__ == "__main__":
    check_steady_state("shock_tube6_H2", [0,0,0])
    check_steady_state("shock_tube7_H2", [4.9e-4, 4.9e-4, 0])

    check_steady_state("soundwave_32_H2", [8.72e-8,]*3)
    check_steady_state("soundwave_64_H2", [2.95e-8,]*3)
    check_steady_state("soundwave_128_H2",[9.07e-9,]*3)
    check_steady_state("soundwave_256_H2",[2.51e-9,]*3)
    check_steady_state("soundwave_512_H2",[6.69e-10,]*3)

    check_steady_state("planet_cartesian_H2", [1.5e-12, 6.2e-14, 7.6e-11])
    check_steady_state("planet_spherical_H2", [2.5e-1, 3.5e+2, 1.4e+2])



    

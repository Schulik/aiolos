#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params


def check_steady_state(problem, L1_target=None):

    IC_file = problem + '_t0.dat'
    snap_file = problem + '_t-666.dat'

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
    check_steady_state("shock_tube6", [0,0,0])
    check_steady_state("shock_tube7", [4.9e-4, 4.9e-4, 0])
    check_steady_state("planet_cartesian", [0.01, 1.7e-10, 0.0])
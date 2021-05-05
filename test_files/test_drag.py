#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params


def check_drag(problem, spc, L1_target=None):

    exact_file = 'output_' + problem + '_an_'+ spc +'_t-1.dat'
    snap_file = 'output_' + problem + '_' + spc + '_t-1.dat'

    exact = load_aiolos_snap(exact_file)
    final = load_aiolos_snap(snap_file)

    L1 = [
        np.abs(final['density'] - 0*exact['density']).mean(),
        np.abs(final['momentum'] - 0*exact['momentum']).mean(),        
        np.abs(final['pressure'] - 0*exact['pressure']).mean(),
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
    check_drag("friction_2spc", 'H0', [3e-16, 3e-14, 2e-14])
    check_drag("friction_2spc_phys", 'H0', [2e-30, 3e-30, 2e-20])

    

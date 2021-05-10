#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import (
    load_aiolos_snap, load_aiolos_params, load_aiolos_species
)

def _solution_2sp(aij_dt, cv1, cv2, T1, T2):
    """Solution to
        cv_i dT_i/dt = aij(Tj - Ti)
    for two species.
    """
    dT = (T2 - T1) * np.exp(-aij_dt * (1/cv1 + 1/cv2)) 

    Tbar = (cv1*T1 + cv2*T2)/(cv1 + cv2)
    
    T1 = Tbar - cv2*dT / (cv1 + cv2)
    T2 = Tbar + cv1*dT / (cv1 + cv2)

    return T1, T2


def check_collisional_heating_2spc(problem, L1_target=None):

    param_file = problem + '.par'
    params = load_aiolos_params(param_file)
    tmax = float(params['PARI_TIME_TMAX'])

    species = load_aiolos_species(problem + '.spc')

    ICs = [load_aiolos_snap('output_' + problem + '_H0_t0.dat'),
           load_aiolos_snap('output_' + problem + '_H1_t0.dat')]

    snaps = [load_aiolos_snap('output_' + problem + '_H0_t-1.dat'),
             load_aiolos_snap('output_' + problem + '_H1_t-1.dat')]


    alpha = float(params['PARI_ALPHA_COLL'])
    cv1 = 0.5 * species[0].dof / species[0].mass
    cv2 = 0.5 * species[1].dof / species[1].mass
    T1 = ICs[0]['temperature'][1:-1].mean()
    T2 = ICs[1]['temperature'][1:-1].mean()

    aij_dt = 3 * alpha * tmax / (species[0].mass + species[1].mass)

    T1_an, T2_an = _solution_2sp(aij_dt, cv1, cv2, T1, T2) 

    T1 = snaps[0]['temperature']
    T2 = snaps[1]['temperature']

    L1 = [np.abs(1-T1[1:-1]/T1_an).sum(), np.abs(1-T2[1:-1]/T2_an).sum()]
 
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
    check_collisional_heating_2spc('collheat_2spc', [1.5e-3, 2.5e-3])

    

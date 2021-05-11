#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import (
    load_aiolos_snap, load_aiolos_params, load_aiolos_species, load_aiolos_diag
)

def check_conservation(problem, species='hydrogen.spc', L1_target=None):
    init = [0,0,0]
    final = [0,0,0]

    params = load_aiolos_params(problem + '.par')
    base_name = "output_" + problem 
    
    try:
        num_spec = int(params['PARI_NUM_SPECIES'])
    except KeyError:
        num_spec = 1
    try:
        num_bound = int(params['PARI_ORDER'])
    except KeyError:
        num_bound = 2


    species = load_aiolos_species(species)
    for i in range(num_spec):

        IC = load_aiolos_snap(base_name + '_'+species[i].name+'_t0.dat')
        init[0] += IC['density'][num_bound:-num_bound].sum()
        init[1] += IC['momentum'][num_bound:-num_bound].sum()
        init[2] += IC['energy'][num_bound:-num_bound].sum()

        snap = load_aiolos_snap(base_name + '_'+species[i].name+'_t-1.dat')
        final[0] += snap['density'][num_bound:-num_bound].sum()
        final[1] += snap['momentum'][num_bound:-num_bound].sum()
        final[2] += snap['energy'][num_bound:-num_bound].sum()


    # Radiative energy:
    if int(params.get('PARI_USE_RADIATION',0)):
        c4pi = 2.99792458e10/(4*np.pi)
        N = int(params.get('NUM_BANDS', 1))
        diag = 'diagnostic_' + problem
        d0 = load_aiolos_diag(diag + '_t0.dat')
        for i in range(N):
            init[2] += d0['J{}'.format(i)][num_bound-1:-num_bound+1].sum()/c4pi

        df = load_aiolos_diag(diag + '_t-1.dat')
        for i in range(N):
            final[2] += df['J{}'.format(i)][num_bound-1:-num_bound+1].sum()/c4pi


    if init[1] != 0:
        L1 = list([abs(1-x2/x1) for x1,x2 in zip(init, final)])
    else:
        L1 = [abs(1-final[0]/init[0]),
              abs(final[1]-init[1]),
              abs(1-final[2]/init[2])]

        
        
    if L1_target is not None:
        passed = True
        for l1, target in zip(L1, L1_target):
            passed &= l1 <= target

        if passed:
            print('Conservation test {} L1 check passed'.format(problem))
        else:
            print('Conservation test {} L1 checked failed:'.format(problem))
            print('\tL1={}, target={}'.format(L1, L1_target))
    else:
        print('Conservation test {} L1 values:'.format(problem))
        print('\tL1={}'.format(L1))


if __name__ == "__main__":
    
    
       
    check_conservation("soundwave_128", L1_target=[5e-16, 2e-11,3e-16]) 
    
    check_conservation("dustywave_stiff", L1_target=[5e-16, 5e-8, 2e-14]) 
    check_conservation("dustywave_nonstiff", L1_target=[4e-16, 2e-10, 2e-14]) 
    
    check_conservation("friction_6spc", "friction_6spc.spc",
                        L1_target=[4e-16, 2e-14, 5e-14])

    check_conservation("friction_2spc", "friction_2spc.spc",
                       L1_target=[4e-16,3e-14,3e-13])
    check_conservation("friction_2spc_phys", "friction_2spc.spc",
                       L1_target=[4e-16,1e-15,3e-14])

    check_conservation("collheat_2spc", "collheat_2spc.spc",
                       L1_target=[4e-16, 4e-14, 2e-13])
    
    check_conservation("collheat_2spc_rad", "collheat_2spc.spc",
                       L1_target=[4e-16, 1e-18, 1e-5])

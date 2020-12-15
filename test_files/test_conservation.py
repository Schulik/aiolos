#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import (
    load_aiolos_snap, load_aiolos_params, load_aiolos_species
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


    L1 = list([abs(1-x2/x1) for x1,x2 in zip(init, final)])

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
    
    
       
    check_conservation("soundwave_128", L1_target=[3e-16, 2e-11,4e-16]) 
    
    check_conservation("dustywave_stiff", L1_target=[0, 3e-9, 1e-14]) 
    check_conservation("dustywave_nonstiff", L1_target=[0, 5e-11, 1e-14]) 
    
    check_conservation("friction_6spc", "friction_6spc.spc",
                        L1_target=[0, 1e-14, 2e-14])

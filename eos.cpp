#include <cmath>

#include "aiolos.h"
#include "eos.h"

/////
///// Adiabatic Equation of state
/////
void Adiabatic_EOS::compute_primitive(const AOS* cons, AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].density = cons[i].u1 ;
        prim[i].speed = cons[i].u2 / cons[i].u1 ;
        prim[i].pres = (_gamma-1) *(cons[i].u3 - 0.5*cons[i].u2*prim[i].speed);
    }
}

void Adiabatic_EOS::compute_conserved(const AOS_prim* prim, AOS* cons, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        cons[i].u1 = prim[i].density ;
        cons[i].u2 = prim[i].density * prim[i].speed ;
        cons[i].u3 = prim[i].pres/(_gamma-1) + 0.5 * cons[i].u2 * prim[i].speed ;
    }
}

void Adiabatic_EOS::compute_auxillary(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].internal_energy = (prim[i].pres/prim[i].density)/(_gamma-1) ;
        prim[i].sound_speed = std::sqrt(_gamma*(_gamma-1)*prim[i].internal_energy) ;
        prim[i].temperature = prim[i].internal_energy / _cv ;
    }
}

/////
///// Polytropic Equation of state
/////
void Polytropic_EOS::compute_primitive(const AOS* cons, AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].density = cons[i].u1 ;
        prim[i].speed = cons[i].u2 / cons[i].u1 ;
        prim[i].pres = _K0 * std::pow(prim[i].density, 1 - 1/_gamma) ;
    }
}

void Polytropic_EOS::compute_conserved(const AOS_prim* prim, AOS* cons, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        cons[i].u1 = prim[i].density ;
        cons[i].u2 = prim[i].density * prim[i].speed ;
        cons[i].u3 = prim[i].pres/(_gamma-1) + 0.5 * cons[i].u2 * prim[i].speed ;
    }
}

void Polytropic_EOS::compute_auxillary(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].internal_energy = prim[i].pres/(prim[i].density*(_gamma-1)) ;
        prim[i].sound_speed = std::sqrt(_gamma*(_gamma-1)*prim[i].internal_energy) ;
        prim[i].temperature = prim[i].internal_energy / _cv ;
    }
}

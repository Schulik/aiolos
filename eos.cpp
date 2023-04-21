/**
 * eos.cpp
 * 
 * This file contains routines implementing different equations of state.
 */
#include <cmath>

#include "aiolos.h"
#include "eos.h"


/**
 * Adiabatic Equation of state, takes conservatives, computes primitives
 * 
 * @param[in] cons pointer to an external conserved variable object
 * @param[out] prim pointer to an external primitive variable object
 * @param[in] num_cells number of cells on which to iterate on both the in&out object
 */
void IdealGas_EOS::compute_primitive(const AOS* cons, AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].density        = cons[i].u1;
        prim[i].speed          = cons[i].u2 / cons[i].u1 ;
        prim[i].pres           = _gamma_m1 *(cons[i].u3 - 0.5*cons[i].u2*prim[i].speed);
    
        if(debug2>3)
            cout<<" In IDEAL GAS EOS i = "<<i<<", primitives = dens/speed/pres = "<<prim[i].density<<"/"<<prim[i].speed<<"/"<<prim[i].pres<<endl;
    }
    
}

/**
 * Adiabatic Equation of state, takes primitives, computes conserved
 * 
 * @param[out] cons pointer to an external conserved variable object
 * @param[in] prim pointer to an external primitive variable object
 * @param[in] num_cells number of cells on which to iterate on both the in&out object
 */
void IdealGas_EOS::compute_conserved(const AOS_prim* prim, AOS* cons, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        cons[i].u1 = prim[i].density;
        cons[i].u2 = prim[i].density * prim[i].speed ;
        cons[i].u3 = prim[i].pres/_gamma_m1 + 0.5 * cons[i].u2 * prim[i].speed ;
        
        if(debug2>3)
            cout<<" In IDEAL GAS EOS i = "<<i<<", conserved = dens/mom/e_tot = "<<cons[i].u1<<"/"<<cons[i].u2<<"/"<<cons[i].u3<<" speed/pres/gamma = "<<prim[i].speed<<"/"<<prim[i].pres<<"/"<<_gamma_m1<<endl;
    }
}

/**
 * Adiabatic Equation of state, takes primitives, computes auxiliaries on those same primitives and writes them into the same objects
 * 
 * @param[in] cons pointer to an external conserved variable object
 * @param[in] num_cells number of cells on which to iterate on both the in&out object
 */
void IdealGas_EOS::compute_auxillary(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].number_density  = prim[i].density / _mass ;
        prim[i].internal_energy = (prim[i].pres/prim[i].density)/_gamma_m1 ;
        prim[i].sound_speed = std::sqrt((_gamma_m1+1)*_gamma_m1*prim[i].internal_energy) ;
        prim[i].temperature = prim[i].internal_energy / _cv ;
        
        if(debug2>3)
            cout<<" In IDEAL GAS EOS i = "<<i<<", aux = number/e_int/c_s/temper = "<<prim[i].number_density<<"/"<<prim[i].internal_energy<<"/"<<prim[i].sound_speed<<"/"<<prim[i].temperature<<" dens/speed/pres = "<<prim[i].density<<"/"<<prim[i].speed<<"/"<<prim[i].pres<<endl;
    }
}


/**
 * Adiabatic Equation of state, assumes the internal energy has changed and updates pressure in the same object to be consistent
 * 
 * @param[in] cons pointer to an external primitive variable object which to update
 * @param[in] num_cells number of cells on which to iterate on both the in&out object
 */
void IdealGas_EOS::update_p_from_eint(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].pres = prim[i].internal_energy * prim[i].density * _gamma_m1;
    }
}

/**
 * Adiabatic Equation of state, assumes the temperature has changed and updates internal energy in the same object to be consistent
 * 
 * @param[in] cons pointer to an external primitive variable object which to update
 * @param[in] num_cells number of cells on which to iterate on both the in&out object
 */
void IdealGas_EOS::update_eint_from_T(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].internal_energy = prim[i].temperature * _cv;
    }
}

/**
 * Adiabatic Equation of state, compute pressure over rho using the equation of state-appropriate function of T (assuming there is no further dependency on rho in this EOS)
 * 
 * @param[in] temperature pointer to an external temperature object. Single value, arrays won't work.
 * @param[out] returnval pointer to external resulting p/rho
 */
void IdealGas_EOS::get_p_over_rho_analytic(const double* temperature, double* returnval) const {
    
    *returnval = *temperature * _gamma_m1 * _cv ;
}

/**
 * Adiabatic Equation of state, compute pressure using the equation of state-appropriate function of T and rho
 * 
 * @param[in] density pointer to an external density object. Single value, arrays won't work.
 * @param[in] temperature pointer to an external temperature object. Single value, arrays won't work.
 * @param[out] pressure pointer to external resulting p
 */
void IdealGas_EOS::get_p_from_rhoT(const double* density,const double* temperature, double* pressure) const {
    
    *pressure = *density * *temperature *  _gamma_m1 * _cv;
    
}

/////
///// Polytropic Equation of state
/////
void Polytropic_EOS::compute_primitive(const AOS* cons, AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].density = cons[i].u1;
        prim[i].speed = cons[i].u2 / cons[i].u1 ;
        prim[i].pres = _K0 * std::pow(cons[i].u1, _gamma) ;
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
        prim[i].number_density  = prim[i].density / _mass ;
        prim[i].internal_energy = prim[i].pres/(prim[i].density*(_gamma-1)) ;
        prim[i].sound_speed     = std::sqrt(_gamma*(_gamma-1)*prim[i].internal_energy) ;
        prim[i].temperature     = prim[i].internal_energy / _cv ;
    }
}


void Polytropic_EOS::update_p_from_eint(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].pres = 0.;
    }
}

void Polytropic_EOS::update_eint_from_T(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].internal_energy = 0.;
    }
}


void Polytropic_EOS::get_p_over_rho_analytic(const double* temperature, double* returnval) const {
    
    *returnval = 0. ;
}

void Polytropic_EOS::get_p_from_rhoT(const double* density,const double* temperature, double* pressure) const {
    
    *pressure = 0.;
    
}

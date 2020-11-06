#include <cmath>

#include "aiolos.h"
#include "eos.h"

/////
///// Adiabatic Equation of state
/////
void IdealGas_EOS::compute_primitive(const AOS* cons, AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].density        = cons[i].u1;
        prim[i].speed          = cons[i].u2 / cons[i].u1 ;
        prim[i].pres           = _gamma_m1 *(cons[i].u3 - 0.5*cons[i].u2*prim[i].speed);
    }
}

void IdealGas_EOS::compute_conserved(const AOS_prim* prim, AOS* cons, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        cons[i].u1 = prim[i].density;
        cons[i].u2 = prim[i].density * prim[i].speed ;
        cons[i].u3 = prim[i].pres/_gamma_m1 + 0.5 * cons[i].u2 * prim[i].speed ;
    }
}

void IdealGas_EOS::compute_auxillary(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].number_density  = prim[i].density / _mass ;
        prim[i].internal_energy = (prim[i].pres/prim[i].density)/_gamma_m1 ;
        prim[i].sound_speed = std::sqrt((_gamma_m1+1)*_gamma_m1*prim[i].internal_energy) ;
        prim[i].temperature = prim[i].internal_energy / _cv ;
    }
}

void update_cons_prim_after_friction(AOS* cons, AOS_prim* prim, double dEkin, double newv, double mass, double gamma, double cv) {
    double temp_u3 = cons->u3;
    
    cons->u3 += dEkin;
    cons->u2 = newv*mass;
    
    prim->speed           = newv;
    prim->pres            = (gamma-1) *(cons->u3 - 0.5*cons->u2*newv);
    prim->internal_energy = prim->pres/(prim->density*(gamma-1)) ;
    prim->sound_speed     = std::sqrt(gamma*(gamma-1)*prim->internal_energy) ;
    prim->temperature     = prim->internal_energy / cv ;
    
    if(prim->pres<0) {
        char a;
        
        cout<<"        Negative pressure in update_after_friction! initial u3 = "<<temp_u3<<" du3 = "<<dEkin<<endl;
        cout<<"        -0.5*u2*speed = "<<(0.5*cons->u2*newv)<<endl;
        cout<<"        newv = "<<newv<<" u1 = "<<cons->u1<<endl; 
        cin>>a;
        
    }
        
    
}

void IdealGas_EOS::update_p_from_eint(AOS_prim* prim, int num_cells) const {
    for (int i=0; i < num_cells; i++) {
        prim[i].pres = prim[i].internal_energy * prim[i].density * _gamma_m1;
    }
}


void IdealGas_EOS::get_p_over_rho_analytic(const double* temperature, double* returnval) const {
    
    *returnval = *temperature * _gamma_m1 * _cv ;
}

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
        prim[i].pres = _K0 * std::pow(cons[i].u1, 1 - 1/_gamma) ;
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


void Polytropic_EOS::get_p_over_rho_analytic(const double* temperature, double* returnval) const {
    
    *returnval = 0. ;
}

void Polytropic_EOS::get_p_from_rhoT(const double* density,const double* temperature, double* pressure) const {
    
    *pressure = 0.;
    
}

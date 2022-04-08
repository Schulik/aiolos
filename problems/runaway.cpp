
#include <stdexcept>

#include "aiolos.h"

void c_Species::user_boundary_left(std::vector<AOS>&) {
    throw std::runtime_error(
        "hydro::user_boundary_left() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}


void c_Species::user_boundary_right(std::vector<AOS>&) {
     double Mdot_max = - 0.01 * 5.97e27 / (3600*24*365.25) ;

    int N = num_cells + 2 - base->num_ghosts ;
    
    AOS ub = u[N-1] ;
    double mdot = 4 * pi * ub.u2 * base->x_i12[N-1] * base->x_i12[N-1] ;
    if (mdot < Mdot_max)
        ub.u2 *= mdot / Mdot_max ;
    
    for (int j=0; j < base->num_ghosts; j++) {
        u[N+j] = ub;
        base->phi[N+j] = base->phi[N-1] ;
    }
}

void c_Species::user_initial_conditions(){

    double T0 = 1000 ; // Take 500K as the initial temperature

    // Set the planet temperature
    base->T_surface = 1000 ;
    base->L_core = 0 ;
    base->use_planetary_temperature = 0 ;
    
    // First get the boundary condition
    AOS_prim p0 ;
    p0.density = 1 ;
    p0.pres = p0.density * Rgas * T0 / mass_amu ;
    p0.speed = 0 ;
    
    // Assume stratification with the gas:
    // Increase Bondi-radius slightly to over-pressure the gas
    double cs2 = Rgas * T0 / 2 ;
    double Rb = 1.0*G*base->planet_mass/cs2 ;
    std::cout << "Bondi Radius:" << Rb/1.01 << "cm\n";
    std::cout << initial_fraction << "\n";

    double Mtot = 0 ;
    for (int j=0; j <= num_cells+1; j++) {
        double f = std::min(1., std::exp(Rb/base->x_i12[j] - Rb/base->x_i12[base->num_ghosts])) ;
        AOS_prim p = p0 ;
        p.density *= f ;
        p.pres *= f ;
        eos->compute_conserved(&p, &u[j], 1) ;
        Mtot += f * base->vol[j] ;
    }
    
    // Normalize the denisty to 0.1 earth mass
    double f = 0.01 * 5.97e27 * initial_fraction / Mtot ;
    for (int j=0; j <= num_cells+1; j++) {
        u[j].u1 *= f ;
        u[j].u2 *= f ;
        u[j].u3 *= f ;
    }
}

void c_Species::user_opacity() {
    // Set constant opacities.
    if (is_dust_like) {
        // Dust
        for (int j = 0; j < num_cells + 2; j++) {
            for (int b = 0; b < num_bands_in; b++) {
                opacity_twotemp(j, b) = 1 ;
            }
            for (int b = 0; b < num_bands_out; b++) {
                opacity(j, b) = 1 ;
                opacity_planck(j, b) = 1 ;
            }
        }
    } else {
        // Gas
        for (int j = 0; j < num_cells + 2; j++) {
            for (int b = 0; b < num_bands_in; b++) {
                opacity_twotemp(j, b) = 0.01 ;
            }
            for (int b = 0; b < num_bands_out; b++) {
                opacity(j, b) = 0.01 ;
                opacity_planck(j, b) = 0.01 ;
            }
        }
    }
}

void c_Sim::user_loop_function() {} ;
void c_Sim::user_heating_function() {} ;
void c_Sim::user_output_function(int) {} ;
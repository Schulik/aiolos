

#include <stdexcept>

#include "aiolos.h"


static double T_eq = 2000 ;

void c_Species::user_boundary_left(std::vector<AOS>& u) {
    
    AOS_prim p ;
    p.density = initial_fraction * mass_amu * amu ;
    p.pres = p.density * Rgas * T_eq / mass_amu ;
    p.speed = std::abs(u[base->num_ghosts].u2/u[base->num_ghosts].u1) ;
    for (int j=0; j < base->num_ghosts; j++) {
        eos->compute_conserved(&p, &u[j], 1) ;
        base->phi[j] = base->phi[base->num_ghosts] ;
    }
}


void c_Species::user_boundary_right(std::vector<AOS>&) {

    int N = num_cells + 2 - base->num_ghosts ;
    for (int j=0; j < base->num_ghosts; j++) {
        u[N+j] = u[N-1] ;
        base->phi[N+j] = base->phi[N-1] ;
    }
}

void c_Species::user_initial_conditions(){

    double n = initial_fraction ;

    std::cout << "Cross-over mass wind setup:\n"
              << "\tT_eq=" << T_eq << "K\n"
              << "\t n=" << n << "g/cm^3\n"; 
    
    // First get the boundary condition
    AOS_prim p0 ;
    p0.density = initial_fraction * mass_amu * amu ;
    p0.pres = p0.density * Rgas * T_eq / mass_amu ;
    p0.speed = 0 ;
    
    // Assume stratification with the gas:
    // Increase Bondi-radius slightly to over-pressure the gas
    double cs2 = Rgas * T_eq / mass_amu ;
    double Rb = 1.01*G*base->planet_mass/cs2 ;
    std::cout << "Bondi Radius:" << Rb/1.01 << "cm\n";

    for (int j=0; j <= num_cells+1; j++) {
        double f = std::min(1., std::exp(Rb/base->x_i12[j] - Rb/base->x_i12[base->num_ghosts])) ;
        AOS_prim p = p0 ;
        p.density *= f ;
        p.pres *= f ;
        eos->compute_conserved(&p, &u[j], 1) ;
    }
}

void c_Species::user_opacity() {
    throw std::runtime_error(
        "hydro::user_opacity() not implemented. You must provide this "
        "method if you want to use user defined initial conditions.") ;
}

void c_Species::user_species_loop_function() {
    for (int j=0; j <= num_cells+1; j++) {
        prim[j].pres = prim[j].density * Rgas * T_eq / mass_amu ;
    }
    eos->compute_conserved(&prim[0], &u[0], num_cells+2);   
    compute_pressure(u);
}

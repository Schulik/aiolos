
#include <stdexcept>

#include "aiolos.h"

// Olivine vapor pressue
static double P_vap(double T) {
    const double m = 169 * amu ;
    const double L = 3.21e10 ; 
    return 6.72e14 *  std::exp(-m*L/(kb*T)) ;
}

static double rho_vap(double T) {
   return 30 * P_vap(T) / (Rgas * T) ;
}

static double T_eq ;

void c_Species::user_boundary_left(std::vector<AOS>& u) {
    
    AOS_prim p ;
    p.density = initial_fraction * rho_vap(T_eq) ;
    p.pres = p.density * Rgas * T_eq / mass_amu ;
    p.speed = std::max(0., u[base->num_ghosts].u2/u[base->num_ghosts].u1) ;
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
    

    // Get the equilibrium temperature
    const double RHO = 1, BETA=1;
    double size = std::pow(3/(4*M_PI)*2.523e12*amu/RHO,1/3.) ;
    
    double k0 = 3/(4*RHO*size) ;
    double T0 = 1/(2*M_PI * 3.475 * size) ;
    
    auto kappa = [BETA, k0, T0](double T) {
        if (BETA == 1) 
            return k0 * T / (T0 + T) ;
        else
            return k0 / (1 + std::pow(T/T0, -BETA)) ;
    } ;

    T_eq = base->T_star ;
    double S = pow(T_eq,4.) * std::pow(base->R_star*rsolar,2.)/std::pow(base->planet_semimajor*au,2.) ;
    double ks = kappa(T_eq) ;
    for (int i=0; i < 100; i++) {
        double k = kappa(T_eq) ;
        T_eq = 0.9*std::pow(S*ks/k/16, 0.25) + 0.1*T_eq ;
    }

    std::cout << "Dusty wind setup: T_eq=" << T_eq << "K\n";
    
    // First get the boundary condition
    AOS_prim p0 ;
    p0.density = initial_fraction * rho_vap(T_eq) ;
    p0.pres = p0.density * Rgas * T_eq / mass_amu ;
    p0.speed = 0 ;
    
    // Assume stratification with the gas:
    // Increase Bondi-radius slightly to over-pressure the gas
    double cs2 = Rgas * T_eq / 30 ;
    double Rb = 1.01*G*base->planet_mass/cs2 ;
    std::cout << "Bondi Radius:" << Rb/1.01 << "cm\n";
    std::cout << initial_fraction << "\n";

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

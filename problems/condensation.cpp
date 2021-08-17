
#include "problems/condensation.h"

int main() {

    double mu_gas = 30 ;
    Condensible silicate(169, 3.21e10, 6.72e14, 0.1) ;

    double Cgas = 1/(1.3-1) * Rgas / mu_gas ;
    double Cl = 1.5 * Rgas / mu_gas ;

    double a_grain = 1e-4 ;
    double rho_grain = 3 ;
    
    SingleGrainCondensation cond(
        silicate, mu_gas, a_grain, rho_grain, Cgas, Cl) ;

    Brent brent ;

    // Initial conditions
    double rho_s = 100 * mu_gas / (Rgas*2200) ;

    SingleGrainCondensation::arr_t T0 = { 2.5e3, 2.0e3} ;
    SingleGrainCondensation::arr_t rho0 = {0.5*rho_s, 0.5*rho_s} ;
    SingleGrainCondensation::arr_t v0 = {1e5, 0e5} ;


    std::vector<double> times = { 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0} ;

    for (auto ti : times) {
        cond.set_state(rho0, T0, v0, ti) ;

        double d0, d1 ; 
        std::tie(d0,d1) = cond.bracket_solution() ;
        
        d0 = brent.solve(d0, d1, cond) ;

        SingleGrainCondensation::arr_t rho, T, v ;
        std::tie(rho, T, v) = cond.update_T_rho_v(d0) ;

        std::cout<< "t=" << ti << "\t"
                 << "T=("<< T[0] << ", " << T[1] << "),\t"
                 << "rho=(" << rho[0] << ", " << rho[1] << "),\t"
                 << "v=(" << v[0] << ", " << v[1] << ")\n" ;

        double E = rho[0]*(Cgas*T[0] + silicate.Lsub)  + rho[1]*Cl*T[1] +
            0.5*rho[0]*v[0]*v[0] + 0.5*rho[1]*v[1]*v[1] ;
        double mom = rho[0]*v[0] + rho[1]*v[1] ;

        std::cout << "\tE_tot=" << E << ", mom_tot="<< mom << "\n" ;
    }

    return 0 ;
} 

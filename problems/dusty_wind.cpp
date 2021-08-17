
#include <stdexcept>

#include "aiolos.h"

#include "problems/condensation.h"



// Setup the condensation parameters 
//    Perez-Becker & Chiang (2013)
Condensible silicate(169, 3.21e10, 6.72e14, 0.1) ;


static double T_eq ;

void c_Species::user_boundary_left(std::vector<AOS>& u) {
    
    AOS_prim p ;
    p.density = initial_fraction * silicate.rho_vap(T_eq) ;
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
    const double RHO = 3, BETA=1;
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
    double S = 0.25 * pow(T_eq,4.) * std::pow(base->R_star*rsolar,2.)/std::pow(base->planet_semimajor*au,2.) ;
    double ks = kappa(T_eq) ;
    double tau0 = 1e300 ;
    for (int i=0; i < 100; i++) {
        double k = kappa(T_eq) ;
        double g = ks/k ;

        double T4 = 0.25*S*g ;//(1 + 3/g + (g - 3/g)*std::exp(-tau0));

        T_eq = 0.9*std::pow(T4, 0.25) + 0.1*T_eq ;
    }

    std::cout << "Dusty wind setup: T_eq=" << T_eq << "K\n";
    
    // First get the boundary condition
    AOS_prim p0 ;
    p0.density = initial_fraction * silicate.rho_vap(T_eq) ;
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

void c_Species::user_species_loop_function() {} ;



void c_Sim::user_heating_function() {

    if (num_species != 2)
        throw std::runtime_error("Condensation requires num_species=2") ;

    double RHO = 3 ;
    double a_grain = std::pow(3/(4*M_PI)*species[1].mass_amu*amu/RHO,1/3.) ;
    
    SingleGrainCondensation cond(
        silicate, species[0].mass_amu, a_grain, RHO, 
        species[0].cv, species[1].cv) ;
        
    for (int b = 0; b < num_bands_in; b++) {
        radial_optical_depth_twotemp(num_cells + 1,b) = 0.;
    }

    RadiationProperties rad(num_bands_in, num_bands_out) ;


    for (int j = num_cells + 1; j > 0; j--) {

        // Get the local conditions:
        std::array<double, 2> rho = 
            { species[0].prim[j].density, species[1].prim[j].density } ;
        std::array<double, 2> T =   
            { species[0].prim[j].temperature, species[1].prim[j].temperature };
        std::array<double, 2> v =   
            { species[0].prim[j].speed, species[1].prim[j].speed } ;

        if (use_rad_fluxes == 1) {

            // Set the stellar radiation properties:
            for (int b=0; b < num_bands_in; b++) {
                rad.kappa_star[b][0] = species[0].opacity_twotemp(j, b) ;
                rad.kappa_star[b][1] = species[1].opacity_twotemp(j, b) ;

                if (j < num_cells + 1)
                    rad.flux_star[b] = 0.25 * solar_heating(b) * 
                        std::exp(-radial_optical_depth_twotemp(j+1,b)) ;
                else
                    rad.flux_star[b] = 0.25 * solar_heating(b) ;
            }

            // Set the thermal radiation properties
            for (int b=0; b < num_bands_out; b++) {
                rad.J_thermal[b] = Jrad_FLD(j, b) ;

                rad.kappa_thermal[b][0] = species[0].opacity_planck(j, b) ;
                rad.kappa_thermal[b][1] = species[1].opacity_planck(j, b) ;

                rad.f_band[b][0] = compute_planck_function_integral3(l_i_out[b], l_i_out[b+1] ,T[0]) ;
                rad.f_band[b][1] = compute_planck_function_integral3(l_i_out[b], l_i_out[b+1] ,T[1]) ;
            }

            // Cell size (for the optical depth)
            rad.dr = dx[j]; 

            // Solve for the new state
            cond.set_state(rho, T, v, rad, dt) ;
        }
        else {
            // No radiation
            cond.set_state(rho, T, v, dt) ;
        }

        double d0, d1 ; 
        std::tie(d0,d1) = cond.bracket_solution() ;
        
        Brent brent(1e-6*rho[1]) ;
        d0 = brent.solve(d0, d1, cond) ;

        std::tie(rho, T, v) = cond.update_T_rho_v(d0) ;

        // Store the new values
        //     Set the new density and velocity
        species[0].prim[j].density = rho[0] ;
        species[1].prim[j].density = rho[1] ;

        species[0].prim[j].speed = v[0] ;
        species[1].prim[j].speed = v[1] ;

        if (T[0] < 0 || T[1] < 0) {
            std::cout << j << " " << dt 
                      << " (" << rho[0] << " " << rho[1] << ")"
                      << " (" << T[0] << " " << T[1] << ")"
                      << " (" << v[0] << " " << v[1] << ")\n" ;
            std::cout << "\t" << species[0].prim[j].temperature << " " << species[1].prim[j].temperature << "\n" ;

            std::array<double,2> heating, cooling ;
            heating = rad.compute_stellar_heating(rho) ;
            cooling = cond.net_heating_rate(rho, T, v) ;

            std::cout << "\t(" << heating[0] << " " << heating[1] << ")"
                      <<  " (" << cooling[0] << " " << cooling[1] << ")\n" ;

            std::cout << "\t(" << (heating[0] + cooling[0])*dt/ (rho[0]*species[0].cv)
                      <<  " " << (heating[1] + cooling[1])*dt / (rho[1]*species[1].cv) << ")\n" ;
        }


        //     Compute the net heating/cooling rate:
        if (use_rad_fluxes) {
            std::array<double,2> heating, cooling ;

            heating = rad.compute_stellar_heating(rho) ;
            cooling = cond.net_heating_rate(rho, T, v) ;

            for (int s=0; s < 2; s++) {
                species[s].dS(j) = heating[s] ;
                species[s].dG(j) = cooling[s] ;
            }

            //     Save the new optical depth etc
            for (int b=0; b < num_bands_in; b++) {

                total_opacity_twotemp(j, b) = rho[0]*rad.kappa_star[b][0] + rho[1]*rad.kappa_star[b][1] ;
                cell_optical_depth_twotemp(j, b) = dx[j]*total_opacity_twotemp(j,b) ;

                if (j < num_cells+1)
                    radial_optical_depth_twotemp(j,b) = cell_optical_depth_twotemp(j, b) +
                        radial_optical_depth_twotemp(j+1, b) ;
                else
                    radial_optical_depth_twotemp(j,b) = cell_optical_depth_twotemp(j, b) ;

                S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j,b)) ;
                if (j == num_cells+1)
                    dS_band(j,b) = 0.25 * solar_heating(b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
                else
                    dS_band(j,b) = 0.25 * S_band(j+1,b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
            }
        }
        else {
            // Without radiation, just set the new temperatues
            for (int s=0; s < 2; s++) 
                species[s].prim[j].temperature = T[s] ;
        }
    }

    // Finally, lets update the conserved quantities
    for (int s = 0; s < num_species; s++) {
        species[s].eos->update_eint_from_T(&species[s].prim[0], num_cells + 2);
        species[s].eos->update_p_from_eint(&species[s].prim[0], num_cells + 2);
        species[s].eos->compute_auxillary(&species[s].prim[0], num_cells + 2);
        species[s].eos->compute_conserved(&species[s].prim[0], &species[s].u[0], num_cells + 2);
    }
}


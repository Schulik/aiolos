
#include <stdexcept>

#include "aiolos.h"

#include "problems/condensation.h"



// Setup the condensation parameters 
//    Perez-Becker & Chiang (2013)
Condensible silicate(169, 3.21e10, 6.72e14, 0.1) ;

static const double RHO = 3 ;
static const double THICKNESS = 10 ;

static double T_core ;
static double f_dust = 0e-20;

void c_Sim::user_output_function(int output_counter) {
    std::ofstream file ;
    std::string filename = workingdir+"user_output.dat" ;

    if (output_counter == 0) {
        file.open(filename) ;
        file << "# snap time[s] T_surf[K] L_surf[K]\n" ;
    } else {
        file.open(filename, std::ios::app) ;
    }
    file << output_counter << " " << globalTime << " " << T_core << " " << 4*pi*R_core*R_core*F_core << "\n" ;
    double T_F = std::pow(std::abs(F_core)/sigma_rad, 0.25) ;
    if (F_core < 0) T_F *= -1 ;
    cout << "T_core, Net thermal flux[K]: " << T_core << ", " << T_F << "K\n";
}

void c_Species::user_boundary_left(std::vector<AOS>& u) {
    
    for (int j=0; j < base->num_ghosts; j++) {
        u[j] = u[2*base->num_ghosts-(j+1)] ;
        u[j].u2 *= -1 ;
        base->phi[j] = base->phi[base->num_ghosts] ;
    }

}


void c_Species::user_boundary_right(std::vector<AOS>& u) {

    int N = num_cells + 2 - base->num_ghosts ;
    for (int j=0; j < base->num_ghosts; j++) {
        u[N+j] = u[N-1] ;
        base->phi[N+j] = base->phi[N-1] ;
    }
}

void c_Species::user_initial_conditions(){
    

    // Get the equilibrium temperature
    const double BETA=1;
    
    double size = std::pow(3/(4*M_PI)*30*7.568e11/3*amu/RHO,1/3.) ;
    
    double k0 = 3/(4*RHO*size) ;
    double T0 = 1/(2*M_PI * 3.475 * size) ;
    
    auto kappa = [BETA, k0, T0](double T) {
        if (BETA == 1) 
            return k0 * T / (T0 + T) ;
        else
            return k0 / (1 + std::pow(T/T0, -BETA)) ;
    } ;

    T_core = base->T_star ;
    double S = 0.25 * pow(T_core,4.) * std::pow(base->R_star*rsolar,2.)/std::pow(base->planet_semimajor*au,2.) ;
    double ks = kappa(T_core) ;
    double tau0 = 0e300 ;
    for (int i=0; i < 100; i++) {
        double k = kappa(T_core) ;
        double g = ks/k ;

        double T4 = 0.25*S*g * (1 + 3/g + (g - 3/g)*std::exp(-tau0));

        T_core = 0.9*std::pow(T4, 0.25) + 0.1*T_core ;
    }
    //T_core = pow(S, 0.25) ;
    // Set the planet temperature
    base->T_surface = T_core ;
    base->L_core = 4*pi*base->R_core*base->R_core*sigma_rad*std::pow(T_core, 4) ;
    base->use_planetary_temperature = 1 ;

    std::cout << "Dusty wind setup:\n\tT_eq=" << T_core << "K\n"
              << "\tGrain size=" << size << "cm\n" ;
    
    // First get the boundary condition
    AOS_prim p0 ;
    p0.pres = initial_fraction *(30/mass_amu)* silicate.P_vap(T_core) ;
    p0.density = p0.pres * mass_amu / (Rgas * T_core) ;
    p0.speed = 0 ;
    
    // Assume stratification with the gas:
    // Increase Bondi-radius slightly to over-pressure the gas
    double cs2 = Rgas * T_core / 30 ;
    double Rb = 1.1*G*base->planet_mass/cs2 ;
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


void c_Sim::user_heating_function() {

    if (num_species != 2)
        throw std::runtime_error("Condensation requires num_species=2") ;

    double a_grain = std::pow(3/(4*M_PI)*species[1].mass_amu*amu/RHO,1/3.) ;
    
    SingleGrainCondensation cond(
        silicate, species[0].mass_amu, a_grain, RHO, 
        species[0].cv, species[1].cv) ;
        
    for (int b = 0; b < num_bands_in; b++) {
        radial_optical_depth_twotemp(num_cells + 1,b) = 0.;
    }

    RadiationProperties rad(num_bands_in, num_bands_out) ;


    for (int j = num_cells + 1; j >= 0; j--) {

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
        
        Brent brent(1e-10*rho[1]) ;
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

            heating = {0, 0 }; // rad.compute_stellar_heating(rho) ;
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
                    radial_optical_depth_twotemp(j,b) = 0 ;

                S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j,b)) ;
                if (j == num_cells+1)
                    dS_band(j,b) = 0.25 * solar_heating(b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
                else
                    dS_band(j,b) = 0.25 * S_band(j+1,b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
            }
            // Same for the thermal bands
            for (int b=0; b < num_bands_out; b++) {
                total_opacity(j, b) = rho[0]*rad.kappa_thermal[b][0] + rho[1]*rad.kappa_thermal[b][1] ;
                cell_optical_depth(j, b) = dx[j]*total_opacity(j,b) ;

                if (j < num_cells+1)
                    radial_optical_depth(j,b) = cell_optical_depth(j, b) + radial_optical_depth(j+1, b) ;
                else
                    radial_optical_depth(j,b) = 0 ;
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

    // Recompute radiative heating etc given new densities.
    update_dS() ;
}

void c_Sim::user_loop_function() {

    if (num_species != 2)
        throw std::runtime_error("Condensation requires num_species=2") ;

    // Compute evaporation and condensation from the surface

    RadiationProperties rad(num_bands_in, 1) ;

    int j = num_ghosts ;
    SurfaceCondensation planet_surf(silicate, species[0].mass_amu, 
                                   species[0].cv, species[1].cv, RHO*THICKNESS,surf[j-1]/vol[j]);
    if (use_rad_fluxes == 1) {
        // Stellar radiation reaching the planet
        for (int b=0; b < num_bands_in; b++) 
            rad.flux_star[b] = 0.25 * solar_heating(b) * 
                    std::exp(-radial_optical_depth_twotemp(j,b)) ;

        // Thermal radiation reaching the planet (F_down = \pi J in diffusion limit)
        rad.J_thermal[0] = - F_core / pi ;
    } else {
        throw std::runtime_error("Radiation must be turned on") ;
    }

    // Energy due to rain out
    double L_rain = 0 ;
    if (species[1].prim[j].speed < 0) {
        AOS_prim& prim = species[1].prim[j] ;

        double E_tot = prim.density*
            (species[1].cv*prim.temperature + 0.5*prim.speed*prim.speed) ;

        L_rain = - E_tot * prim.speed ;

        prim.density += prim.density * prim.speed * dt * surf[j-1]/vol[j] ;
    } 

    // Get the surface temperature and mass-loss rate
    planet_surf.set_state(T_core,
                          species[0].prim[j].temperature, 
                          species[0].prim[j].pres,
                          rad, L_rain, dt) ;

    double T0, T1 ; 
    std::tie(T0,T1) = planet_surf.bracket_solution() ;
    
    Brent brent(1e-10*T_core) ;

    T_core = brent.solve(T0, T1, planet_surf) ;
    
    L_core = 4*pi*R_core*R_core*sigma_rad*std::pow(T_core, 4) ;

    /*
    std::cout 
        << "t, dt, T " << globalTime << " " << dt << " " << T_core << ", " 
        << rad.flux_star[0]<< " " << pi*rad.J_thermal[0] << " "
        << L_core/(4*pi*R_core*R_core) << " " << L_rain << "\n" ;
    */

    double drho = dt*planet_surf.mass_flux(T_core)*surf[j-1]/vol[j] ;
    double heat = dt*planet_surf.gas_heating_rate(T_core)*surf[j-1]/vol[j] ;

    std::array<double, 2> frac = {1-f_dust, f_dust} ;

    for (int s=0; s < 2; s++) {
        // Update the mass, momentum and energy of the corresponding cell
        AOS_prim& prim = species[s].prim[j] ;

        double E_tot = prim.density*
            (species[s].cv*prim.temperature + 0.5*prim.speed*prim.speed) ;

        E_tot += frac[s]*heat;

        prim.speed /= (1 + frac[s]*drho/prim.density) ;
        prim.density += drho ;
        
        E_tot -= 0.5*prim.density*prim.speed*prim.speed ;

        prim.temperature = E_tot / (prim.density*species[s].cv) ;
    }

    for (int s = 0; s < num_species; s++) {
        species[s].eos->update_eint_from_T(&species[s].prim[j], 1);
        species[s].eos->update_p_from_eint(&species[s].prim[j], 1);
        species[s].eos->compute_auxillary(&species[s].prim[j], 1);
        species[s].eos->compute_conserved(&species[s].prim[j], &species[s].u[j], 1);
    }
} ;

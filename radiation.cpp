///////////////////////////////////////////////////////////
//
//
//  radiation.cpp
//
// This file contains routines computing the absorption of solar radiation, radiation transport and updating temperatures accordingly.
//
//
//
///////////////////////////////////////////////////////////

//#define EIGEN_RUNTIME_NO_MALLOC

#include <cassert>
#include "aiolos.h"

//////////////////////////////////////////////////////////////////
//
// void update_opacities():
//      Computation of opacities based on material properties.
//      Furthermore the attenuation of solar radiation is ray-traced,
//      flux loss per cell per band is computed and assigned to individual species based on 
//
//
//
//
//
//
//
//
//////////////////////////////////////////////////////////////////
void c_Sim::reset_dS() {
    for(int j = num_cells + 1; j>0; j--)
        for(int s=0; s<num_species; s++)
                species[s].dS(j)  = 0.;
}

void c_Sim::update_dS() {
    
    //
    // Compute optical depths, solar radiation attenuation and heating function for low-energy bands
    //
    for(int j = num_cells + 1; j>0; j--) {
        
        //////////////////////////////////////////////////////////////////
        //////////////Solar bands
        //////////////////////////////////////////////////////////////////
        for(int b=0; b<num_bands_in; b++) {
            
            if(!BAND_IS_HIGHENERGY[b]) {
            
                total_opacity_twotemp(j,b)        = 0 ;
                radial_optical_depth_twotemp(j,b) = 0.;
                
                for(int s=0; s<num_species; s++) 
                    total_opacity_twotemp(j,b)  += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 ;
                
                cell_optical_depth_twotemp(j,b) = total_opacity_twotemp(j, b) * dx[j] ;
                
                if(j==num_cells+1)
                    radial_optical_depth_twotemp(j,b) = 0.; //cell_optical_depth(j,b);
                else
                    radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + cell_optical_depth_twotemp(j,b);
                
                //
                // After the total optical depth per band is known, we assign the fractional optical depths
                // Maybe merge with previous loop for optimization
                //
                
                for(int s=0; s<num_species; s++) 
                    species[s].fraction_total_solar_opacity(j,b) = species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth_twotemp(j,b);
                    //species[s].fraction_total_opacity(j,b) = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth(j,b);
                
                //
                // Now compute the attenuation of solar radiation and then assign the lost energy back to individual species in a manner that conserves energy
                //
                
                S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j,b));
                
                //if(steps == 0)
                //    cout<<"IN SOLAR HEATING, const_opacity_solar_factor = "<<const_opacity_solar_factor<<endl;
                
                if(debug >= 1)
                    cout<<" in cell ["<<j<<"] band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<" tau_rad(j,b) = "<<radial_optical_depth(j,b)<<" exp(tau) = "<<std::exp(-radial_optical_depth(j,b))<<endl;
                
                // Regular dS computation
                if(steps > -1) {
                            
                    if(j<num_cells+1)
                        dS_band(j,b) = 0.25 * solar_heating(b) * (-exp(-radial_optical_depth_twotemp(j+1,b)) * expm1(-cell_optical_depth_twotemp(j,b)) )/dx[j] ;
                    else
                        // Use optically thin limit  
                        dS_band(j,b) = 0.25 * solar_heating(b)*total_opacity_twotemp(j,b);
                    
                    dS_band_zero(j,b) = dS_band(j,b);
                
                    
                }// Irregular dS computation, in case we want to fix the solar heating function to its initial value
                else 
                    dS_band(j,b) = dS_band_zero(j,b);
                //
                // In low-energy bands, individual species are heated according to their contribution to the total cell optical depth
                //
                for(int s=0; s<num_species; s++)
                    species[s].dS(j)  += dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
            }
        }
        
        //////////////////////////////////////////////////////////////////
        ///////////////////Thermal bands
        //////////////////////////////////////////////////////////////////
        
        for(int b=0; b<num_bands_out; b++) {
            
                cell_optical_depth(j,b)           = 0.;
                total_opacity(j,b)                = 0 ;
                radial_optical_depth(j,b)         = 0.;
                
                for(int s=0; s<num_species; s++) {
                    total_opacity(j,b)          += species[s].opacity(j,b)         * species[s].u[j].u1 ; //TODO: Replace this with sensible addition rule for Rosseland opacities!
                }
                cell_optical_depth(j,b)         = total_opacity(j, b)         * dx[j] ;
                
                if(j==num_cells+1) 
                    radial_optical_depth(j,b)         = 0.;
                else
                    radial_optical_depth(j,b)         = radial_optical_depth(j+1,b)         + cell_optical_depth(j,b);
        }
        
    }
    
    //Initialize optically thin regions with low energy density
    if((steps==0) && (init_J_factor > 0.)) {
        
        cout<<" Setting J to custom values. init_J_factor = "<<init_J_factor<<" and init_T_temp ="<<init_T_temp<<" const_opacity_solar_factor = "<<const_opacity_solar_factor<<endl;
        
        for(int j = num_cells + 1; j>0; j--) {
            for(int b=0; b<num_bands_out; b++) {
                
                if(radial_optical_depth(j,b) < 1.) 
                    Jrad_FLD(j,b) *= init_J_factor/pow(x_i12[j]/x_i12[34] ,2.);
                
                //for(int s=0; s<num_species; s++)
                  //  species[s].prim[j].temperature = init_T_temp;
                
            }
        }
    }
    
    
    for(int s=0; s<num_species; s++) 
        if (species[s].const_T_space > 0) {
            species[s].prim[num_cells].temperature   = species[s].const_T_space;
            species[s].prim[num_cells+1].temperature = species[s].const_T_space;
        }


    //Compute surface temperature from core flux
    double J_remnant = 0;
    double T_target;
    
    /*
    for(int b = 0; b < num_bands_in; b++) {
        J_remnant += 0.25 * S_band(num_ghosts-1,b);
    }
    */
    
    J_remnant += L_core/(4*pi*R_core*R_core) ;

    T_target  = pow(std::abs(J_remnant) / sigma_rad, 0.25);
    if (J_remnant < 0) 
       T_target *= -1 ; // Use T < 0 to represent negative flux

    T_surface = T_target;
    
    //Set planetary temperature (computed in update_opacities)
    if(false && use_planetary_temperature == 1) {
        
        for(int s=0; s<num_species; s++) {
            species[s].prim[1].temperature = T_surface;
            species[s].prim[0].temperature = T_surface;
        }
        for(int b=0; b<num_bands_out; b++) {
            Jrad_FLD(1, b) = sigma_rad * pow(T_surface,4.) * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], std::abs(T_surface)) ;
            Jrad_FLD(0, b) = sigma_rad * pow(T_surface,4.) * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], std::abs(T_surface)) ;
        }
        //cout<<"t ="<<globalTime<<" Jrad(0,b) = "<<Jrad_FLD(0, 0)<<" Tsurface = "<<T_surface<<" Tcore = "<<T_core<<" T_rad_remnant = "<<T_rad_remnant<<endl;
    }

    
    if(debug >= 1) {
        
        for(int b=0; b<num_bands_in;b++) {
           cout<<" in band_in ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<endl;
        }
        cout<<"    in update opa "<<endl;
        for(int j=num_cells; j>0; j--) {
            
            for(int b=0; b<num_bands_in;b++) {
                cout<<"    band ["<<b<<"][j="<<j<<"] = "<<S_band(j,b)<<" dS = "<<dS_band(j,b)<<" tau = "<<const_opacity_solar_factor*radial_optical_depth_twotemp(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fract["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_solar_opacity(j,b))<<" rho/T = "<<species[s].prim[j].density<<"/"<<species[s].prim[j].temperature<<" x = "<<x_i12[j];
                    //" opa = "<<species[s].opacity(j,b)<<" total opa(j+1)"<<total_opacity(j+1,b)<<
                }
                //cout<<endl;
            }
            cout<<endl;
        }
        
        if(debug > 3) {
            char stop;
            cin>>stop;    
        }
        
    }
}


bool c_Sim::update_fluxes_FLD(double dt_step) {
    
    // If test, then set J to initial values before computing on
    if(radiation_matter_equilibrium_test == 1) {
        
        for (int j=0; j < num_cells+2; j++) {
            for(int b=0; b<num_bands_out; b++) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
        }
    }
    
    // Diffusion coefficient for the core:
    std::vector<double> D_core(num_bands_out, 0)  ;

    int num_vars = num_bands_out + num_species;
    int stride = num_vars * num_vars ;
    int size_r = (num_cells + 2) * num_vars ;
    int size_M = (num_cells + 2) * stride ;

    std::vector<double> 
        l(size_M, 0.), d(size_M, 0.), u(size_M, 0.), r(size_r, 0.) ;
    
    // Step 1: setup transport terms (J)
    for(int b=0; b<num_bands_out; b++) {
        for (int j=0; j < num_cells+2; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;

            // Time dependent terms:
            d[idx] +=  vol[j] / (c_light * dt_step) ;
            r[idx_r] += (vol[j] / (c_light * dt_step)) * Jrad_FLD(j, b) ;

            // Flux across right boundary
            if (j > 0 && j < num_cells + 1) {
                double dx      = (x_i12[j+1]-x_i12[j]) ;
                double rhokr   = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx );
                       rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);
                double tau_inv = 1 / (dx * rhokr) ;
                //double tau_inv = 0.5 / (dx * (total_opacity(j,b) + total_opacity(j+1,b))) ;
                double R       = 2 * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300) ;
                double D       = no_rad_trans * surf[j] * flux_limiter(R) * tau_inv;

                // divergence terms
                u[idx] = -D ;
                d[idx] += D ;
                d[idx+stride] = D ;
                l[idx+stride] = -D ;
                
                if(debug > 1)
                    cout<<" radiation part 0. t,j,b="<<steps<<","<<j<<","<<b<<" tau_inv/R/D = "<<tau_inv<<"/"<<R<<"/"<<D<<" J/J/dJ = "<<Jrad_FLD(j+1,b)<<"/"<<Jrad_FLD(j,b)<<"/"<<(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<" flux = "<<D*(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<endl;
            }
        }
        
        // Boundaries:
        // Left boundary:
        //    Reflecting / no flux or planetary temperature
        for (int j=0; j < num_ghosts; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;

            // Compute heating due to radiation reaching surface
            if (j == num_ghosts-1 && use_planetary_temperature) {
                double S ;
                S = sigma_rad * pow(T_surface,4.) * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], std::abs(T_surface)) ;  
                if (T_surface < 0) 
                    S *= -1 ; 

                double dx      = (x_i12[j+1]-x_i12[j]) ;
                double rhokr   = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx );
                       rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);
                double tau_inv = 1 / (dx * rhokr) ;
                double R       = 2 * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300) ;

                D_core[b]    = no_rad_trans * flux_limiter(R) * tau_inv;
                double Chi   = 0.125*(1 + 3*eddington_coeff(R)) ;
                
                l[idx] = 0 ;
                d[idx] =  0.5*(D_core[b] + Chi) * surf[j] ;
                u[idx] = -0.5*(D_core[b] - Chi) * surf[j] ;
                r[idx_r] = surf[j] * S / (4*pi) ;
            } else {
                l[idx] = 0 ;
                u[idx] = -d[idx] ;
                r[idx_r] = 0 ;
            }
        }
        
        //   Right boundary: reflective?
        //if(geometry == Geometry::cartesian) {
        if(closed_radiative_boundaries) {

            int Ncell = num_cells - 2*(num_ghosts - 1);
            for (int j=0; j < num_ghosts; j++) {
                int i = Ncell + num_ghosts + j ;

                int idx = i*stride + b*(num_vars + 1) ;
                int idx_r = i*num_vars + b ;    

                l[idx] = -d[idx] ;
                u[idx] = 0 ;
                r[idx_r] = 0 ;
            }
        }
        else {//   Right boundary: free stream, no emission / absorbtion.

            // Only need to set the very last cell.
            //  Assume F = J and \div(F) = const

            int idx = (num_cells+1)*stride + b*(num_vars + 1) ;
            int idx_r = (num_cells+1)*num_vars + b ;  
            
            double f = x_i12[num_cells]/x_i12[num_cells+1] ;
            switch (geometry) {
                case Geometry::cartesian:
                    f = 1 ;
                    break;
                case Geometry::cylindrical:
                    break;
                case Geometry::spherical:
                    f *= f;
                    break;
            }
            
            l[idx] = -f*d[idx] ;
            u[idx] = 0;
            r[idx_r] = 0 ;
        }
        
        if(debug >= 1) {
            
            for(int index=0; index < num_cells+2; index++) {    
//                 /int index = (num_cells/2+1);
                
                cout<<" radiation part1, t = "<<steps<<" band["<<b<<"] cell["<<index<<"] l/d/u/r = "<<l[index]<<"/"<<d[index]<<"/"<<u[index]<<"/"<<r[index];
                cout<<" temps = ";
                for(int si = 0; si<num_species; si++) {
                        cout<<species[si].prim[2].temperature<<" ";
                }
                cout<<endl;
            }
            
            if(debug > 3) {
                char a;
                cin>>a;
                
            }
        }
    }

    // Step 2: Energy exchange terms kappa*rho*(J-B) + dS + Pi + Lambda
    
    if(radiation_matter_equilibrium_test <= 2) { //radtests 3 and 4 are delta-radiation peaks without energy-matter coupling
        
        for (int j=0; j < num_cells+2; j++) {
            for (int s=0; s < num_species; s++) {
                
                if(debug > 1) {
                        cout<<" Going into radpart 2. ";
                        cout<<" opa    = "<<species[s].opacity_planck(j, 0);
                        cout<<" temper = "<<species[s].prim[j].temperature;
                        cout<<" dens   = "<<species[s].prim[j].density<<endl;
                }
                
                int idx_s  = j*stride + (s + num_bands_out) * (num_vars+1) ;
                int idx_rs = j*num_vars + (s + num_bands_out) ;
                
                
                double Ts = species[s].prim[j].temperature ;
                double rhos = species[s].prim[j].density ;

                d[idx_s ] = 1 / dt_step ;
                r[idx_rs] = Ts / dt_step ;
                r[idx_rs] += (species[s].dS(j) + species[s].dG(j)) / species[s].u[j].u1 / species[s].cv;
                
                for(int b=0; b<num_bands_out; b++) {
                    int idx_b  = j*stride + b * (num_vars+1) ;
                    int idx_bs = j*stride + b * num_vars + (s + num_bands_out) ; 
                    int idx_sb = j*stride + (s + num_bands_out) * num_vars + b ;
                    int idx_rb = j*num_vars + b ; 
                        
                    double fac = no_rad_trans * species[s].opacity_planck(j, b) * sigma_rad*Ts*Ts*Ts / pi * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], species[s].prim[j].temperature);
                    
                    d[idx_s ] += 16 * pi * fac / species[s].cv ;
                    d[idx_sb] = - 4 * pi * species[s].opacity_planck(j, b)  * no_rad_trans / species[s].cv ;
                    r[idx_rs] += 12 * pi * fac * Ts / species[s].cv ;
                    
                    if (j < num_ghosts || j >= num_cells + 2-num_ghosts) continue ;
                    //cout<<" in rad matrix["<<j<<"]: term1 = "<<12 * pi * fac * Ts / species[s].cv<<" term2 = "<<species[s].dS(j)  / species[s].u[j].u1 / species[s].cv<<endl;
                    
                    d[idx_b ] += vol[j] * rhos * species[s].opacity_planck(j, b) * no_rad_trans;
                    d[idx_bs] = - 4 * vol[j] * rhos * fac ;
                    r[idx_rb] -=  3 * vol[j] * rhos * fac * Ts ;
                    
                    if(debug > 1)
                        cout<<" radiation part2, t = "<<steps<<" b["<<b<<"] i["<<j<<"] s["<<s<<"] l/d/u/r_rs/rs_rb = "<<d[idx_s]<<"/"<<d[idx_sb]<<"/"<<d[idx_b]<<"/"<<d[idx_bs]<<"/"<<r[idx_rs]<<"/"<<r[idx_rb]<<" opac = "<<vol[j] * rhos * species[s].opacity(j, b)<<endl;
                }
            }
            if (use_collisional_heating && num_species > 1) {
                fill_alpha_basis_arrays(j);
                compute_collisional_heat_exchange_matrix(j);

                for (int si=0; si < num_species; si++) { 
                    int idx  = j*stride + (si + num_bands_out) * num_vars + num_bands_out;
                    for (int sj=0; sj < num_species; sj++) {
                        d[idx+sj] -= friction_coefficients(si,sj) ;
                    }
                }
            }
        }
    }
    if(debug > 4) {
        char stepstop;
        cin>>stepstop;
    }

    tridiag.factor_matrix(&l[0], &d[0], &u[0]) ;
    tridiag.solve(&r[0], &r[0]) ; 

    // Check whether the solution worked (J, T >= 0)
    //    If not, repeat with a smaller time-step
    bool passed = true  ;
    double Jmin = 0, Tmin = 0 ;
    int jmin, bmin ;
    for (int j=0; j < num_cells+2; j++) 
        for(int b=0; b<num_bands_out+num_species; b++)
            if (r[j*num_vars + b] < 0) {
                passed = false ;
                if (b < num_bands_out) {
                    if (r[j*num_vars + b] < Jmin) {
                        Jmin = r[j*num_vars + b] ;
                        jmin = j ;
                        bmin = b ;
                    }
                }
                else {
                    if (r[j*num_vars + b] < Tmin) {
                        Tmin = r[j*num_vars + b] ;
                        jmin = j ;
                        bmin = b - num_bands_out ;
                    }
                }
            }
    
    if (not passed) {
        if (dt_step*1024 <= dt) { 
            std::cout << "Radiation step failed, reached too small timestep fraction\n" ;
            return false ;
        }
        else {
            std::cout << "Radiation failed, halving step:\n\t"
                      << "dt_rad=" << dt_step << ", ratio=" << dt_step/dt
                      << ", Tmin=" << Tmin << ", Jmin=" << Jmin << "\n"
                      << "\tjmin=" << jmin << ", bmin=" << bmin 
                      << std::endl ;
            
            if (dt_step == dt)
                F_core = 0 ;
                
            // If the sub-steps fail we will just use the single full-step, so
            // pass the failure back up
            passed = update_fluxes_FLD(dt_step/2) ;
            if (passed) 
                passed = update_fluxes_FLD(dt_step/2) ;

            // Exit early if we've passed because the sub-steps have updated the
            // radiation and temperature
            if (passed)            
                return true ;
        }
    }



    // Store the result
    for (int j=0; j < num_cells+2; j++) {
        for(int b=0; b<num_bands_out; b++) {
            Jrad_FLD(j, b) = r[j*num_vars + b] ;
            
            if(radiation_matter_equilibrium_test == 1) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
                    
          //  std::cout << j << "\t" <<  Jrad_FLD(j, b) << "\n" ;
        }

        for(int s=0; s<num_species; s++) {
            species[s].prim[j].temperature = r[j*num_vars + (s + num_bands_out)] ;
        }
    }
    
    if (use_planetary_temperature) {
        if (dt == dt_step) F_core = 0 ;
        for(int b=0; b<num_bands_out; b++) 
            F_core -= 4*pi * D_core[b] * (dt_step / dt) * (Jrad_FLD(num_ghosts, b) - Jrad_FLD(num_ghosts-1, b)) ;
    }
    
    // Update energies. 
    // TODO: We should add cv * (Tf - Ti) to u to conserve energy properly.
    for(int si=0; si<num_species; si++) {
        species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }

    return passed ;
}


void c_Sim::update_temperatures_simple() {
    
    //Compute change in energy
    for (int j=0; j < num_cells+2; j++) {
        for(int s=0; s<num_species; s++) {
            species[s].prim[j].temperature      += (species[s].dS(j) + species[s].dG(j) ) * dt / species[s].u[j].u1 / species[s].cv ;
            species[s].prim[j].internal_energy  += (species[s].dS(j) + species[s].dG(j) ) * dt / species[s].u[j].u1;
            
            if(steps > 605 && j == 25 && s == 2)
                cout<<" In Updte_T_simple, T, j, steps, s =  "<<species[s].prim[j].temperature<<", "<<j<<", "<<steps<<", "<<s<<" dS, dG = "<<species[s].dS(j)<<", "<<species[s].dG(j)<<" dT = "<<(species[s].dS(j) + species[s].dG(j) ) * dt  / species[s].u[j].u1 / species[s].cv<<endl;
                
            if(species[s].prim[j].temperature < 0) {
                
                cout<<" In Updte_T_simple, T<0!!!, j, steps, s =  "<<j<<", "<<steps<<", "<<s<<" dS, dG = "<<species[s].dS(j)<<", "<<species[s].dG(j)<<" dT = "<<(species[s].dS(j) + species[s].dG(j) ) * dt  / species[s].u[j].u1 / species[s].cv<<endl;
            }
        }
            
    }
    
    // Update Temperatures
    for(int si=0; si<num_species; si++) {
            species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
            species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
            species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);       
            
            //For self-consistency
            species[si].compute_pressure(species[si].u);
    }
    
    //Debug
    for (int j=0; j < num_cells+2; j++) {
        for(int s=0; s<num_species; s++) {
            if(steps > 605 && j == 25 && s == 2)
                cout<<" AFTER Updte_T_simple, T =  "<<species[s].prim[j].temperature<<", "<<" p = "<<species[s].prim[j].pres<<" E = "<<species[s].u[j].u3<<" e = "<<species[2].prim[25].internal_energy<<endl;
        }
    }
    
}

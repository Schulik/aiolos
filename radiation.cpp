/** 
 * radiation.cpp
 * 
 * This file contains routines computing the absorption of solar radiation, radiation transport and updating temperatures accordingly.
 */

#define EIGEN_RUNTIME_NO_MALLOC

#include <cassert>
#include "aiolos.h"

/**
 * Resets the heating and cooling functions and applies a time-dependent ramp up to the top-of-the atmosphere radiation spectrum.
 */
void c_Sim::reset_dS() {
    for(int j = num_cells + 1; j>0; j--) {
        for(int s=0; s<num_species; s++) {
            species[s].dS(j)    = 0.;
            species[s].dG(j)    = 0.;
            species[s].dGdT(j)  = 0.;
        }
    }
        
    if(globalTime < radiation_rampup_time) {
        for(int b=0; b<num_bands_in; b++) {
        
            solar_heating(b) = solar_heating_final(b) * ( init_radiation_factor + (1.-init_radiation_factor) * globalTime/radiation_rampup_time );
        }
    } else {
        for(int b=0; b<num_bands_in; b++) {
            solar_heating(b) = solar_heating_final(b);
        }        
    }
                
}

/**
 * Integrate the opacity to stellar irradiation in band b over cell j for all contributions from species s. Update radial optical depth. Called before update_dS_jb().
 */
void c_Sim::update_tau_s_jb(int j, int b) {
    
    total_opacity_twotemp(j,b)        = 0 ;
    radial_optical_depth_twotemp(j,b) = 0.;
                
    for(int s=0; s<num_species; s++) 
        total_opacity_twotemp(j,b)  += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 ;
                
    cell_optical_depth_twotemp(j,b) = total_opacity_twotemp(j, b) * dx[j] ;
                
    if(j==num_cells+1)
        radial_optical_depth_twotemp(j,b) = 0.; //cell_optical_depth(j,b);
    else
        radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + cell_optical_depth_twotemp(j,b);
    
}

/**
 * Update heating function. Keep track of highenergy- and thermal opacities.
 */
void c_Sim::update_dS_jb(int j, int b) {
    
                //
                // After the total optical depth per band is known, we assign the fractional optical depths
                // Maybe merge with previous loop for optimization
                //
                
                for(int s=0; s<num_species; s++) 
                    species[s].fraction_total_solar_opacity(j,b) = highenergy_switch(s,b) * species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth_twotemp(j,b);
                    //species[s].fraction_total_opacity(j,b) = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth(j,b);
                
                //
                // Now compute the attenuation of solar radiation and then assign the lost energy back to individual species in a manner that conserves energy
                //
                if(j>=num_cells+1)
                    S_band(j,b) = solar_heating(b);
                else
                    S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b));
                //S_band(j,b) = solar_heating(b) * fastexp2(-radial_optical_depth_twotemp(j,b));
                
                if(debug >= 1)
                    cout<<" in cell ["<<j<<"] band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<" tau_rad(j,b) = "<<radial_optical_depth(j,b)<<" exp(tau) = "<<std::exp(-radial_optical_depth(j,b))<<endl;
                
                //
                // Compute the heating dS as the difference between flux in and out of the cell
                // Contains the geometric averaging factor 1/4 as well as a crude approximation of albedo effects, as we do not solve for scattering.
                // Different approximations of exp and expm1 can improve speed in cases when many bands are integrated. The current setting is fast and accurate enough.
                //
                // Highenergy optical depths are computed separately in chemistry.cpp lines ~706 and used in line 336. This is because the heating rates need to be split onto the
                // ionisation products, not the absorbing reactants, so there the heating rates need to be split.
                //
                if(steps > -1) {
                            
                    if(j<num_cells+1){
                        //dS_band(j,b) = 0.25 * solar_heating(b) * fastexp2(-radial_optical_depth_twotemp(j+1,b))  /dx[j] * (1.-bond_albedo);
                        //dS_band(j,b) = 0.25 * solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b))  /dx[j] * (1.-bond_albedo);
                        dS_band(j,b) = 0.25 * S_band(j,b)  /dx[j] * (1.-bond_albedo);
                        double dtau_tot = -(radial_optical_depth_twotemp(j+1,b) - radial_optical_depth_twotemp(j,b));
                        double dS_he_temp = dS_band(j,b);
                        
                        if(dtau_tot > 1e-3)
                            dS_band(j,b) *= (-expm1(-dtau_tot));
                        else
                            dS_band(j,b) *= (-fastexpm1_2(- dtau_tot));
                            //dS_band(j,b) *= (-expm1(-const_opacity_solar_factor * dtau_tot));
                            
                        if( cell_optical_depth_highenergy(j,b) > 1e-3)
                            dS_he_temp *=  -expm1(-cell_optical_depth_highenergy(j,b)) ;
                        else 
                            dS_he_temp *=  -fastexpm1_2(-cell_optical_depth_highenergy(j,b)) ;
                            //dS_he_temp *=  -expm1(-const_opacity_solar_factor * cell_optical_depth_highenergy(j,b)) ;
                        
                        //cout<<" j/b = "<<j<<"/"<<b<<" dS_band(j,b) "<<dS_band(j,b)<<" dS_he_temp "<<dS_he_temp<<endl;
                        dS_band(j,b) -= dS_he_temp; //Substract highenergy heating again, as this has already been added in update_dS_jb_photochem() around line 700 in chemistry.cpp
                        dS_band(j,b) = std::max(dS_band(j,b), 0.); //Safeguard against highenergy shenanigans
                        if(BAND_IS_HIGHENERGY[b] && photochemistry_level==1)
                            dS_band(j,b) = 0.;
                        
                        if(b==999) {
                            cout<<" in thermal heating: total heating ~ "<<dS_band(j,b)<<" of which highenergy is "<<dS_he_temp;
                            cout<<" resulting in dS_band("<<j<<","<<b<<") = "<<dS_band(j,b)<<" F = "<<solar_heating(b)<< " F*e-tau = "<<solar_heating(b) * exp(-radial_optical_depth_twotemp(j+1,b))<<" tau "<<dtau_tot<<" ";
                            for(int s=0; s<num_species; s++) {
                                cout<<species[s].u[j].u1<<" ";
                            }
                            cout<<endl;
                            
                        }
                        
                    }
                    else
                        // Use optically thin limit  
                        dS_band(j,b) = dS_band(j-1,b); //0.25 * solar_heating(b)*total_opacity_twotemp(j,b)*(1-bond_albedo);
                    
                    dS_band_zero(j,b) = dS_band(j,b);
                
                    //
                    // Planetary heating 2
                    //
                    if(use_planetary_temperature == 1){
                        double lum = 1.0 * sigma_rad * T_core*T_core*T_core*T_core * 0.5; //Factor 0.5 is not understood. Use proper radiation boundary conditions for avoiding this.
                        //Spread the luminosity for fewer crashes
                        //dS_band(2,b) += 3./6. * lum / (dx[2]); 
                        //dS_band(3,b) += 2./6. * lum / (dx[3]); 
                        //dS_band(4,b) += 1./6. * lum / (dx[4]); 
                        
                        dS_band(2,b) += lum * surf[2]/ (vol[2]); 
                    }
                    
                }// Irregular dS computation, in case we want to fix the solar heating function to its initial value
                else {
                    dS_band(j,b) = dS_band_zero(j,b);
                }
                    
                //
                // In low-energy bands, individual species are heated according to their contribution to the total cell optical depth
                //
                for(int s=0; s<num_species; s++) {
                    //if lowenergy or photochem < 2
                    
                    if(photochemistry_level <= 2) {
                        species[s].dS(j)  += highenergy_switch(s,b) * dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                        if(species[s].dS(j) < 1e-50)
                            species[s].dS(j) = 0.;
                        
                        if(species[s].dG(j) > -1e-50)
                            species[s].dG(j) = 0.;
                    }
                        
                    else {
                        
                        //if(species.s participates in photoreaction)
                                    //take educt.dS(j,b) and assign it to products, weighed with (1-E_lim/hv) and mass ratios
                        
                        //for all reacs:
                            //for educt:
                        //         Equant = dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                        //    for products:
                        //         dS += Equant * mass_weight;
                        
                    }
                    
                    
                    /*
                    else {
                        
                        //Assuming one photoreaction only acts on one reactant and photon band
                        for(c_photochem_reaction& reaction : photoreactions) {
                            threshold_energy = 13.6; //TODO: replace with reaction-specific number
                            heating_mask_educ  = {0.};
                            heating_mask_prod  = {0.,1.};
                            
                            double energy_available = dS_band(j,b) * (1 - threshold_energy * ev_to_K * kb / photon_energies[b]);
                            double cooling = nX[2]*HOnly_cooling(nX, Tx(2));
                            
                            for(int& ej : reaction.educts ) {
                                
                                species[ej].dS(j)  += heating_mask[ej] * species[ej].fraction_total_solar_opacity(j,b) * energy_available;
                                species[ej].dG(j)  += heating_mask[ej]
                            }
                            for(int& pj : reaction.products ) {
                                
                                species[pj].dS(j)  += heating_mask[pj] * species[pj].fraction_total_solar_opacity(j,b) * energy_available;
                            }                            
                        
                        }
                    
                    }
                    */
                }
                    
    
    
    
    
}

/**
 * Wrapper function for computing the heating function, without solving for the new temperatures just yet.
 * By the time this is called from the main loop, photo/thermo chemistry is already done and the photoionisation rates 
 * and photoheating rates (contributions of dS_band due to cell_optical_depth_highenergy(cell, b) are already added to dS_band.
 */
void c_Sim::update_dS() {
    
    if(debug > 1)
        cout<<"Updating dS..."<<endl;
    
    //
    // Compute optical depths, solar radiation attenuation and heating function for low-energy bands
    //
    for(int j = num_cells + 1; j>0; j--) {
        
        //////////////////////////////////////////////////////////////////
        //////////////Stellar irradiation bands
        //////////////////////////////////////////////////////////////////
        for(int b=0; b<num_bands_in; b++) {
            
            update_tau_s_jb(j, b); //See above in this file
            update_dS_jb(j, b);   //See above in this file
        }
        
        //NOTE: update_dS_jb_photochem(j,b) is now moved into chemistry.cpp, due to algorithmic requirements.
        
        //////////////////////////////////////////////////////////////////
        ///////////////////Thermal bands
        //////////////////////////////////////////////////////////////////
        
        for(int b=0; b<num_bands_out; b++) {
            
                cell_optical_depth(j,b)           = 0.;
                total_opacity(j,b)                = 0 ;
                radial_optical_depth(j,b)         = 0.;
                
                for(int s=0; s<num_species; s++) {
                    //total_opacity(j,b)          += species[s].opacity(j,b)         * species[s].u[j].u1 ; //TODO: Replace this with sensible addition rule for Rosseland opacities!
                    total_opacity(j,b)          = std::fmax(species[s].opacity(j,b) * species[s].u[j].u1, total_opacity(j,b)) ; //TODO: Replace this with sensible addition rule for Rosseland opacities!
                }
                cell_optical_depth(j,b)         = total_opacity(j, b)         * dx[j] ;
                
                if(j==num_cells+1) 
                    radial_optical_depth(j,b)         = 0.;
                else
                    radial_optical_depth(j,b)         = radial_optical_depth(j+1,b)         + cell_optical_depth(j,b);
        }
        
    }
    
    //For tests of the radiation solver when initialising the radiation field out of equilibrium.
    if((steps==0) && (init_J_factor > 1e10)) {
        
        cout<<" Setting J to custom values. init_J_factor = "<<init_J_factor<<" and init_T_temp ="<<init_T_temp<<" const_opacity_solar_factor = "<<const_opacity_solar_factor<<endl;
        
        for(int j = num_cells + 1; j>0; j--) {
            for(int b=0; b<num_bands_out; b++) {
                
                if(radial_optical_depth(j,b) < 1.) 
                    Jrad_FLD(j,b) *= init_J_factor; //*std::pow(x_i[j]/x_i[34] , -2.);
                
                for(int s=0; s<num_species; s++)
                    species[s].prim[j].temperature = init_T_temp;
                
            }
        }
    }
    
    if(debug >= 1) {
        
        for(int b=0; b<num_bands_in;b++) {
           cout<<" in band_in ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<endl;
        }
        cout<<"    in update opa "<<endl;
        for(int j=num_cells+1; j>=0; j--) {
            
            for(int b=0; b<num_bands_in; b++) {
                cout<<"    band ["<<b<<"][j="<<j<<"] = "<<S_band(j,b)<<" dS = "<<dS_band(j,b)<<" tau_s = "<<radial_optical_depth_twotemp(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fract["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_solar_opacity(j,b))<<" opas = "<<species[s].opacity_twotemp(j,b)<<" rho = "<<species[s].prim[j].density;
                     
                     //cout<<" dS_fract["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_solar_opacity(j,b))<<" rho/T/Etot-Ekin = "<<species[s].prim[j].density<<"/"<<species[s].prim[j].temperature<<"/"<<(species[s].u[j].u3 - 0.5*std::pow(species[s].u[j].u2,2.)/species[s].u[j].u1)<<" x = "<<x_i12[j];
                    //" opa = "<<species[s].opacity(j,b)<<" total opa(j+1)"<<total_opacity(j+1,b)<<
                }
                //cout<<endl;
            }
            cout<<endl;
        }
        
        for(int j=num_cells+1; j>=0; j--) {
            
            for(int b=0; b<num_bands_out; b++) {
                cout<<"    band ["<<b<<"][j="<<j<<"] = "<<S_band(j,b)<<" dS = "<<dS_band(j,b)<<" tau_s = "<<radial_optical_depth_twotemp(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" opap["<<species[s].speciesname<<"] = "<<(species[s].opacity_planck(j,b))<<" opa = "<<species[s].opacity(j,b);
                     
                }
            }
            cout<<endl;
        }
        
        if(debug > 3) {
            char stop;
            cin>>stop;    
        }
        
    }
    
    //End of update_dS
}

/**
 * Routine containing the actual radiation transport solver.
 */
void c_Sim::update_fluxes_FLD() {
    
    if(debug > 1)
        cout<<"Starting update_fluxes_FLD.."<<endl;
    
    // If test, then set J to initial values before computing on
    if(radiation_matter_equilibrium_test == 1) {
        
        for (int j=0; j < num_cells+2; j++) {
            for(int b=0; b<num_bands_out; b++) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
        }
    }
    
    auto flux_limiter = [](double R) {
        if (R <= 2)
            return 2 / (3 + std::sqrt(9 + 10*R*R)) ;
        else 
            return 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
    } ;
    //A few other possible flux limiters. Use with caution.
    /* 
    auto flux_limiter = [](double R) {
        if (R <= 1.5)
            return 2 / (3 + std::sqrt(9 + 12*R*R)) ;
        else 
            return 1. / (1. + R + std::sqrt(1. + 2.*R)) ;
    } ;*/ 
    
    /*
     auto flux_limiter = [](double R) {
         
         return (2.+ R) / (6. + 3*R + R*R) ;
     } ;*/
    
     /*
     auto flux_limiter = [](double R) {
         
         return 1. / (3.+ R) ;
     } ;*/
     
    int num_vars = num_bands_out + num_species;
    int stride = num_vars * num_vars ;
    int size_r = (num_cells + 2) * num_vars ;
    int size_M = (num_cells + 2) * stride ;

    std::vector<double> 
        l(size_M, 0.), d(size_M, 0.), u(size_M, 0.), r(size_r, 0.) ;
    //
    
    // Step 1: setup transport terms (J)
    for(int b=0; b<num_bands_out; b++) {
        for (int j=0; j < num_cells+1; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;

            // Time dependent terms:
            d[idx] +=  vol[j] / (c_light * dt) ;
            r[idx_r] += (vol[j] / (c_light * dt)) * Jrad_FLD(j, b) ;

            // Flux across right boundary
            if (j > 0 && j < num_cells + 1) {
                double dx      = (x_i12[j+1]-x_i12[j]) ;                
                //double rhokr   = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx ); //From Ramsey Dullemond 2015
                //double rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);
                double rhokr = std::sqrt(total_opacity(j,b)*total_opacity(j+1,b));

                double tau_inv = 1. / (dx * rhokr) ;
                double R       = xi_rad * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j, b) + 1e-300) ; // Put in 1.0 as prefactor to get correct rad shock
                double D       = 1. * tau_inv * 1. * surf[j] * flux_limiter(R);
                
                //if(1./tau_inv < 1.e-13)
                //    D = 0.;
                
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
        //if(use_planetary_temperature == 0) {
            for (int j=0; j < num_ghosts; j++) {
                int idx = j*stride + b*(num_vars + 1) ;
                int idx_r = j*num_vars + b ;
                
                l[idx] = 0 ;
                u[idx] = -d[idx] ;
                //u[idx] = 0. ;
                //d[idx] =  ;
                r[idx_r] = 0 ; 
            }
        //}
        
        //   Right boundary: reflective?
        //if(geometry == Geometry::cartesian) {
        if(closed_radiative_boundaries) {
        //if(true) {

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
                //cout<<endl;
            }
            
            if(debug > 3) {
                char a;
                cin>>a;
                
            }
        }
    }

    // Step 2: Energy exchange terms kappa*rho*(J-B) + dS + Pi + Lambda
    
    if(radiation_matter_equilibrium_test <= 2) { //radtests 3 and 4 are delta-radiation peaks without energy-matter coupling
        
        for (int j=0; j < num_cells+1; j++) {
            for (int s=0; s < num_species; s++) {
                
                if(debug > 1) {
                        cout<<" Going into radpart 2. i["<<j<<"]";
                        cout<<" opas/p/r    = "<<species[s].opacity_twotemp(j, 0)<<"/"<<species[s].opacity_planck(j, 0)<<"/"<<species[s].opacity(j, 0);
                        cout<<" temper = "<<species[s].prim[j].temperature;
                        cout<<" dens   = "<<species[s].prim[j].density<<" dS = "<<species[s].dS(j)<<endl;
                }
                
                int idx_s  = j*stride + (s + num_bands_out) * (num_vars+1) ;
                int idx_rs = j*num_vars + (s + num_bands_out) ;
                
                double Ts = species[s].prim[j].temperature ;
                double rhos = species[s].prim[j].density ;

                d[idx_s ] = 1 / dt ;
                r[idx_rs] = Ts / dt ;
                r[idx_rs] +=  (species[s].dS(j) + species[s].dG(j)) / species[s].u[j].u1 / species[s].cv; //Misc heating terms that are not directly related to self-consistent temperature are just added to the rhs
                
                for(int b=0; b<num_bands_out; b++) {
                    int idx_b  = j*stride + b * (num_vars+1) ;
                    int idx_bs = j*stride + b * num_vars + (s + num_bands_out) ; 
                    int idx_sb = j*stride + (s + num_bands_out) * num_vars + b ;
                    int idx_rb = j*num_vars + b ; 
                        
                    double fac = 1.* no_rad_trans * species[s].opacity_planck(j, b) * sigma_rad*Ts*Ts*Ts / pi * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], species[s].prim[j].temperature);
                    
                    d[idx_s ] += 16 * pi * fac / species[s].cv ;
                    d[idx_sb] = - 4 * pi * species[s].opacity_planck(j, b) * no_rad_trans / species[s].cv ;
                    r[idx_rs] += 12 * pi * fac * Ts / species[s].cv ;
                    
                    if (j < num_ghosts || j >= num_cells + 2-num_ghosts) continue ;
                    //cout<<" in rad matrix["<<j<<"]: term1 = "<<12 * pi * fac * Ts / species[s].cv<<" term2 = "<<species[s].dS(j)  / species[s].u[j].u1 / species[s].cv<<endl;
                    
                    d[idx_b ] += vol[j] * rhos * species[s].opacity_planck(j, b) * no_rad_trans;
                    d[idx_bs] = - 4 * vol[j] * rhos * fac ;
                    r[idx_rb] -=  3 * vol[j] * rhos * fac * Ts ;
                    
                    if(debug >= 1)
                        cout<<" radiation part2, t = "<<steps<<" b["<<b<<"] i["<<j<<"] s["<<s<<"] l/d/u/r_rs/rs_rb = "<<d[idx_s]<<"/"<<d[idx_sb]<<"/"<<d[idx_b]<<"/"<<d[idx_bs]<<"/"<<r[idx_rs]<<"/"<<r[idx_rb]<<" opac = "<<vol[j] * rhos * species[s].opacity(j, b)<<" fac = "<<fac<<endl;
                }
            }
            if (use_collisional_heating && num_species > 1) {
                
                double tau = total_opacity(j,0) * (x_i12[j+1]-x_i12[j]);
                
                 //Heuristic fix for boundary heat bug at large optical depth
                    fill_alpha_basis_arrays(j);
                    compute_collisional_heat_exchange_matrix(j);

                    for (int si=0; si < num_species; si++) { 
                        int idx  = j*stride + (si + num_bands_out) * num_vars + num_bands_out;
                        for (int sj=0; sj < num_species; sj++) {
                            
                            if (tau < 1e3) {
                            
                            d[idx+sj] -= friction_coefficients(si,sj) ;
                            
                            } else
                                d[idx+sj] -= friction_coefficients(si,sj) * 1.e-4 ;
                            //cout<<" si/sj = "<<si<<"/"<<sj<<" coeff = "<<friction_coefficients(si,sj);
                        }
                    }
                
                
                //char a;
                //cin>>a;
            }
        }
    }
    
    
        // Making space for the convective energy transport, following Tajima & Nakagawa 1997
        // Lconv = 2pir^2 c_p dT**3/2 std::sqrt(rho g Lmabda \partial rho/\partial T_P=const )
        //         dT     = Lambda (dT'-dT)/2
        //         Lambda = P/dP
        
        //
        // Step3: Transport terms for convective fluxes in the T-equation
        //
        if(use_convective_fluxes) {
            
           auto smooth = [](double dT, double dT_ad) {
                
                double dTrel = (dT - dT_ad)/dT_ad;
                return (dT-dT_ad) *1e-20 * (-std::expm1(-dTrel*dTrel));//; (1.-std::exp(-dTrel*dTrel));
                
                if(dT - dT_ad > 0.)
                   return (dT-dT_ad) *0.01 * (1.-std::exp(-dTrel*dTrel));
                else
                    return 0.;
            } ;
        
           for (int j=2; j < num_cells-1; j++) {
                for (int s=0; s < num_species; s++) {
                    
                    int idx      = j*stride   + (s + num_bands_out) * (num_vars+1) ;
                    int idx_r    = j*num_vars + (s + num_bands_out) ;
                    
                    double dx = (x_i12[j+1]-x_i12[j]) ;
                    double rhoavg = std::sqrt(species[s].prim[j].density * species[s].prim[j+1].density) ; //Averages as in Kutter&Sparks 1972  
                    double Pavg   = std::sqrt(species[s].prim[j].pres * species[s].prim[j+1].pres);
                    double Tavg   = std::sqrt(species[s].prim[j].temperature * species[s].prim[j+1].temperature);
                    double dP     = (species[s].prim[j].pres - species[s].prim[j+1].pres)/Pavg;
                    double glocal = -get_phi_grav(x_i[j], enclosed_mass[j])/x_i[j];
                    
                    double dT     = (species[s].prim[j].temperature - species[s].prim[j+1].temperature)/Tavg / dP;
                    double nabla_ad = 1.-1./species[s].gamma_adiabat;
                    
                    double lam = Pavg / (species[s].prim[j].pres - species[s].prim[j+1].pres);  // The mixing length
                    double DT = smooth(dT, nabla_ad);//(dT > nabla_ad ? dT - nabla_ad : 0.); //smooth(dT, nabla_ad); //   // Gradient comparison and switch for Lconv
                           DT = (dx * total_opacity(j,0)) > 2./3. ? DT : DT*(dx * total_opacity(j,0)); //Guardian to not use convection in optically thin areas
                           DT = max(0.,DT);
                    //double DT = dT - nabla_ad;    // Gradient comparison and switch for Lconv
                    
                    double alphaconv = 0.5 * species[s].cv * lam * lam * std::sqrt(DT) * rhoavg * std::sqrt(glocal/Tavg); //Prefactor
                    
                    double Lconv = alphaconv;    //Convection
                    //double kappa_conductive = 23.e2 * pow(Tavg/273., 0.7); //From engineering toolbox and Watson 1981 
                    //double Lconv =  kappa_conductive * dTabs / species[s].cv;  //Conduction
                    
                    if(species[s].is_dust_like)
                        Lconv = 0;
                    
                    u[idx] += -Lconv ;
                    d[idx] += Lconv ;
                    d[idx+stride] += Lconv ;
                    l[idx+stride] += -Lconv ;
                    //r[idx_r]      -= Lconv * nabla_ad;
                    
                    species[s].lconvect[j] = Lconv * dT;
                    
                    if(globalTime > 1e100) {
                    //if(dT > nabla_ad) {
                        cout<<" step "<<steps<<" j = "<<j<<"DT = "<<DT<<" dT = "<<dT<<" nabla= "<<nabla_ad<<" dt= "<<dT*dP<<" dTrel= "<<(dT - nabla_ad)/nabla_ad<<" convective flux nonzero = "<<alphaconv<<endl;
                    }
                        
                    if(globalTime > 1e16) {
                        
                        cout<<j<<" nabla_ad = "<<nabla_ad<<" nabla_actual = "<<dT<<" DT = "<<DT<<" deltaT = "<<(species[s].prim[j].temperature - species[s].prim[j+1].temperature)<<" Lconv = "<<Lconv<<endl;
                        
                        
                    }
                        
                }
            }
            
            if(globalTime > 1e100) {
                char a;
                cin>>a;
            }
            
            // Convective boundaries
            //
            for (int j=0; j < num_ghosts; j++) {
                
                for (int s=0; s < num_species; s++) {
                    int idx      = j*stride   + (s + num_bands_out) * (num_vars+1) ;
                    
                    double Lconv = 0.*species[s].lconvect[2]*1.5;
                    
                    u[idx] += -Lconv ;
                    d[idx] += Lconv ;
                    d[idx+stride] += Lconv ;
                    l[idx+stride] += -Lconv ;
                }
                
            }
            
            if(globalTime > 1e16) {
                
                char a;
                cin>>a;
            }
            
        }
        
    
    if(debug >= 3) {
        
        cout<<"L ="<<endl;
        for(int i = 0; i < size_M; i++) {
            
            if(i%num_vars == 0)
                cout<<endl;
            if(i%stride == 0)
                cout<<endl;
            
            cout<<l.at(i)<<" ";
            
        }
        
        cout<<"D ="<<endl;
        for(int i = 0; i < size_M; i++) {
            
            if(i%num_vars == 0)
                cout<<endl;
            if(i%stride == 0)
                cout<<endl;
            
            cout<<d.at(i)<<" ";
            
        }
        
        cout<<"u ="<<endl;
        for(int i = 0; i < size_M; i++) {
            
            if(i%num_vars == 0)
                cout<<endl;
            if(i%stride == 0)
                cout<<endl;
            
            cout<<u.at(i)<<" ";
            
        }
        
        char stepstop;
        cin>>stepstop;
    }
    if(debug >= 3 ) {
        cout<<" i  l d  u"<<endl;
        for(int i = 0; i < size_M; i++) {
            cout<<i<<" "<<l.at(i)<<" "<<d.at(i)<<" "<<u.at(i)<<" "<<endl;
        }
        char stepstop;
        cin>>stepstop;
    }
    
    //
    // Solve!
    //
    tridiag.factor_matrix(&l[0], &d[0], &u[0]) ;
    tridiag.solve(&r[0], &r[0]) ; // Solve in place (check it works)

    // Store the result
    for (int j=0; j < num_cells+2; j++) {
        for(int b=0; b<num_bands_out; b++) {
            Jrad_FLD(j, b) = r[j*num_vars + b] ;
            
            if(radiation_matter_equilibrium_test == 1) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
                
        }

        for(int s=0; s<num_species; s++) {
            double tt = r[j*num_vars + (s + num_bands_out)];
            
            if(tt<temperature_floor)
                tt=temperature_floor;
            
            if(tt>max_temperature)
                tt=max_temperature;
            
            species[s].prim[j].temperature = tt ;
        }
    }
    
    
    // Update energies. 
    // TODO: We should add cv * (Tf - Ti) to u to conserve energy properly.
    for(int si=0; si<num_species; si++) {
        species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }

}


/**
 *  Explicit temperature update. Used for tests only, as not very stable.
 */
void c_Sim::update_temperatures_simple() {
    
    //Compute change in energy
    for (int j=0; j < num_cells+2; j++) {
        for(int s=0; s<num_species; s++) {
            species[s].prim[j].temperature      += (species[s].dS(j) + species[s].dG(j) ) * dt / species[s].u[j].u1 / species[s].cv ;
            species[s].prim[j].internal_energy  += (species[s].dS(j) + species[s].dG(j) ) * dt / species[s].u[j].u1;
            
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

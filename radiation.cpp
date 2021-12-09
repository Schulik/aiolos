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

#define EIGEN_RUNTIME_NO_MALLOC

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
    for(int j = num_cells + 1; j>0; j--) {
        for(int s=0; s<num_species; s++) {
            species[s].dS(j)  = 0.;
            species[s].dG(j)  = 0.;
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
    
       
        if(steps==1107) {
            
            cout<<" in reset_dS Steps == 1107, dS = "<<dS_band(60,0)<<" S = "<<S_band(60,0)<<" F = "<<solar_heating(0)<<" F0 = "<<solar_heating_final(0)<<endl;
            
        }   
                
}

void c_Sim::update_dS() {
    
    if(debug > 1)
        cout<<"Updating dS..."<<endl;
    
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
                        dS_band(j,b) = 0.25 * solar_heating(b) * (-exp(-radial_optical_depth_twotemp(j+1,b)) * expm1( const_opacity_solar_factor* ( radial_optical_depth_twotemp(j+1,b) - radial_optical_depth_twotemp(j,b))) ) /dx[j] * (1.-bond_albedo);
                    else
                        // Use optically thin limit  
                        dS_band(j,b) = 0.25 * solar_heating(b)*total_opacity_twotemp(j,b)*(1-bond_albedo);
                    
                    dS_band_zero(j,b) = dS_band(j,b);
                
                    
                    //
                    // Planetary heating 2
                    //
                    if(use_planetary_temperature == 1){
                        
                        //Spread the luminosity for fewer crashes
                        dS_band(2,b) += 3./6. *0.5* pow(T_core,4.)*sigma_rad / (dx[2]); 
                        dS_band(3,b) += 2./6. *0.5* pow(T_core,4.)*sigma_rad / (dx[3]); 
                        dS_band(4,b) += 1./6. *0.5* pow(T_core,4.)*sigma_rad / (dx[4]); 
                        
                    }
                    
                }// Irregular dS computation, in case we want to fix the solar heating function to its initial value
                else {
                    dS_band(j,b) = dS_band_zero(j,b);
                }
                    
                //
                // In low-energy bands, individual species are heated according to their contribution to the total cell optical depth
                //
                for(int s=0; s<num_species; s++)
                    species[s].dS(j)  += dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
            
                if(steps==1107 && j==60) {
            
                        cout<<" in compute dS, Steps == 1107, dS = "<<dS_band(60,0)<<" S = "<<S_band(60,0)<<" F = "<<solar_heating(0)<<" F0 = "<<solar_heating_final(0)<<endl;
                        cout<<" nominal ds = "<<0.25 * solar_heating(b)*total_opacity_twotemp(j,b)*(1-bond_albedo)<<endl;
                        
                }
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
    
    //
    //  bad boundaries
    //
    //for(int s=0; s<num_species; s++) {
    //    species[s].dS(num_cells+2)  = species[s].dS(num_cells+1);   
   // }
    
    //Initialize optically thin regions with low energy density
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
    
    //Compute surface temperature from remaining solar flux and core flux
    double cooling_time  = 1e9 / core_cv;
    double J_remnant = 0;
    double T_rad_remnant;
    double T_target;
    
    for(int b = 0; b < num_bands_in; b++) {
        J_remnant += S_band(1,b);
    }
    
    //
    // To dampen oscillating outer temperature fluctuations
    //
    for(int s=0; s<num_species; s++)  {
        //species[s].prim[num_cells].temperature   = species[s].prim[num_cells-1].temperature; // species[s].const_T_space;
        //species[s].prim[num_cells+1].temperature = species[s].prim[num_cells-1].temperature; //species[s].const_T_space;
    }
    
    for(int b=0; b<num_bands_out; b++)  {
        //Jrad_FLD(num_cells, b)   = 0.9*Jrad_FLD(num_cells-1, b);
        //Jrad_FLD(num_cells+1,b)  = 0.8*Jrad_FLD(num_cells-1, b);
    }
        //if (species[s].const_T_space > 0) {
            
        //}
    
    
    /*T_rad_remnant = pow(pi * J_remnant / sigma_rad, 0.25);
    T_target      = pow(pow(T_core,4.) + pow(T_rad_remnant,4.), 0.25);
    T_surface     = T_target; //+ (T_surface - T_target)*std::exp(-dt/cooling_time);
    
    //Set planetary temperature (computed in update_opacities)
    if(use_planetary_temperature == 1) {
        for(int s=0; s<num_species; s++) {
            species[s].prim[1].temperature = T_surface;
            species[s].prim[0].temperature = T_surface;
        }
        for(int b=0; b<num_bands_out; b++) {
            Jrad_FLD(1, b)            = sigma_rad * pow(T_surface,4.) * pow(R_core,2.) * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], T_surface) ;
            Jrad_FLD(0, b)            = sigma_rad * pow(T_surface,4.) * pow(R_core,2.) * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], T_surface) ;
        }
        //cout<<"t ="<<globalTime<<" Jrad(0,b) = "<<Jrad_FLD(0, 0)<<" Tsurface = "<<T_surface<<" Tcore = "<<T_core<<" T_rad_remnant = "<<T_rad_remnant<<endl;
    }*/
    
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
    
    if(steps==1107) {
            
            cout<<" in update_dS Steps == 1107, dS = "<<dS_band(60,0)<<" S = "<<S_band(60,0)<<" F = "<<solar_heating(0)<<endl;
            
        }
}


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
    
    /*auto flux_limiter = [](double R) {
        if (R <= 2)
            return 2 / (3 + std::sqrt(9 + 10*R*R)) ;
        else 
            return 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
    } ;*/
    
    auto flux_limiter = [](double R) {
        if (R <= 1.5)
            return 2 / (3 + std::sqrt(9 + 12*R*R)) ;
        else 
            return 1. / (1. + R + std::sqrt(1. + 2.*R)) ;
    } ; 
    
    /*
     auto flux_limiter = [](double R) {
         
         return (2.+ R) / (6. + 3*R + R*R) ;
     } ;
    */
     /*
     auto flux_limiter = [](double R) {
         
         return 1. / (3.+ R) ;
     } ;
     */
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
                double rhokr   = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx );
                       rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);

                double tau_inv = 1. / (dx * rhokr) ;
                double R       = 4. * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300) ; // Put in 1.0 as prefactor to get correct rad shock
                double D       = 1. * tau_inv * no_rad_trans * surf[j] * flux_limiter(R);
                
                // restarting iteration with 4 0.25 0.25 (everything was working)
                //titan single species tables works with 4 0.5 0.25
                
                //correct static G10 model prefactors: 4/2/ 0.25.
                //fake-correct supercritical shock prefactors: 4/4/ 3.141/8  or 16/2/ pi/16
                //correct supercrit shock prefactors: 2 0.5 0.5
                
                // divergence terms
                u[idx] = -D ;
                d[idx] += D ;
                d[idx+stride] = D ;
                l[idx+stride] = -D ;
                
                if(debug > 1)
                    cout<<" radiation part 0. t,j,b="<<steps<<","<<j<<","<<b<<" tau_inv/R/D = "<<tau_inv<<"/"<<R<<"/"<<D<<" J/J/dJ = "<<Jrad_FLD(j+1,b)<<"/"<<Jrad_FLD(j,b)<<"/"<<(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<" flux = "<<D*(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<endl;
                
                /*if(j==20 || j==21 || j==22) {
                    if( (base->steps == 3205600) || (base->steps == 3206000) || (base->steps == 3205700) || (base->steps == 3205625) || (base->steps == 3205650) || (base->steps == 3205675) || (base->steps == 3205605) || (base->steps == 3205610) || (base->steps == 3205615) || (base->steps == 3205620) ) {
                    
                    
                    //cout<<"Radiation 1,  R  "<<   
                    
                    
                    }
                }*/
                    
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
            
            int Ncell   = num_cells - 2*(num_ghosts - 1) ;
            for (int j=0; j < num_ghosts-1; j++) {
                int i       = Ncell + num_ghosts + j;

                int idx = i*stride + b*(num_vars + 1) ;
                int idxim1 = (i-1)*stride + b*(num_vars + 1) ;
                int idx_r = i*num_vars + b ;  

                if (j == 0) // FLD flux at the edge of the last cell
                    l[idx] = -surf[i-1] / (x_i12[i]-x_i12[i-1]) ; //u[i-1] ;u[idxim1];//;
                else // Free-stream 
                    l[idx] = -surf[i-1] / (x_i12[i]-x_i12[i-1]) ;
                
                //l[idx] = -surf[i-1] / (x_i12[i]-x_i12[i-1]);// (x_i12[i]-x_i12[i-1]) ;
                
                
                d[idx] = -0.*l[idx] + surf[ i ] / (x_i12[i+1]-x_i12[i]); // (x_i12[i+1]-x_i12[i]) ;//;// 
                u[idx] = 0 ;
                r[idx_r] = 0 ;
            }

            int idx = (num_cells+1)*stride + b*(num_vars + 1) ;
            int idx_r = (num_cells+1)*num_vars + b ;  
            l[idx] = -1;
            d[idx] = 1;
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
                        
                    double fac = no_rad_trans * species[s].opacity_planck(j, b) * sigma_rad*Ts*Ts*Ts / pi * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], species[s].prim[j].temperature);
                    
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
                fill_alpha_basis_arrays(j);
                compute_collisional_heat_exchange_matrix(j);

                for (int si=0; si < num_species; si++) { 
                    int idx  = j*stride + (si + num_bands_out) * num_vars + num_bands_out;
                    for (int sj=0; sj < num_species; sj++) {
                        d[idx+sj] -= friction_coefficients(si,sj) ;
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
                    double rhoavg = (species[s].prim[j].density + species[s].prim[j+1].density) * 0.5;
                    double Pavg   = (species[s].prim[j].pres + species[s].prim[j+1].pres) * 0.5;
                    double Tavg   = (species[s].prim[j].temperature + species[s].prim[j+1].temperature) * 0.5;
                    double dP     = (species[s].prim[j].pres - species[s].prim[j+1].pres)/Pavg;
                    double dT     = (species[s].prim[j].temperature - species[s].prim[j+1].temperature)/Tavg / dP;
                    double glocal = -get_phi_grav(x_i[j], enclosed_mass[j])/x_i[j];
                    
                    double nabla_ad = 1.-1./species[s].gamma_adiabat;
                    double lam = Pavg / (species[s].prim[j].pres - species[s].prim[j+1].pres);  // The mixing length
                    double DT =  (dT > nabla_ad ? dT - nabla_ad : 0.); //smooth(dT, nabla_ad); //   // Gradient comparison and switch for Lconv
                           DT = (dx * total_opacity(j,0)) > 2./3. ? DT : 0.; //Guardian to not use convection in optically thin areas
                    //double DT = dT - nabla_ad;    // Gradient comparison and switch for Lconv
                    
                    double alphaconv = 0.5 * species[s].cv * lam * lam * DT * rhoavg * std::sqrt(glocal/Tavg); //Prefactor
                    double Lconv = alphaconv;
                    double Lconvsymm = 0;
                    
                    u[idx] += -Lconv ;
                    d[idx] += Lconv ;
                    d[idx+stride] += Lconv ;
                    l[idx+stride] += -Lconv ;
                    //r[idx_r]      += Lconv * nabla_ad;
                    
                    species[s].lconvect[j] = Lconv * dT;
                    
                    if(false) {
                    //if(dT > nabla_ad) {
                        cout<<" step "<<steps<<" j = "<<j<<" DT = "<<DT<<" convective flux nonzero = "<<alphaconv<<endl;
                        
                        char a;
                        cin>>a;
                    }
                        
                    if(globalTime > 1e16) {
                        
                        cout<<j<<" nabla_ad = "<<nabla_ad<<" nabla_actual = "<<dT<<" DT = "<<DT<<" deltaT = "<<(species[s].prim[j].temperature - species[s].prim[j+1].temperature)<<" Lconv = "<<Lconv<<endl;
                        
                        
                    }
                        
                }
            }
            
            // Convective boundaries
            //
            for (int j=0; j < num_ghosts; j++) {
                
                for (int s=0; s < num_species; s++) {
                    int idx      = j*stride   + (s + num_bands_out) * (num_vars+1) ;
                    int idx_r    = j*num_vars + (s + num_bands_out) ;
                    
                    double Lconv = species[s].lconvect[2]*1.5;
                    
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
        //cout<<"L ="<<endl<<l<<endl;
        //cout<<"D ="<<endl<<d<<endl;
        //cout<<"U ="<<endl<<u<<endl;
        
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
                    
          //  std::cout << j << "\t" <<  Jrad_FLD(j, b) << "\n" ;
        }

        for(int s=0; s<num_species; s++) {
            double tt = r[j*num_vars + (s + num_bands_out)];
            
            if(tt<temperature_floor)
                tt=temperature_floor;
            
            species[s].prim[j].temperature = tt ;
            //if(j==num_cells)
            //    cout<<" num_cells = "<<num_cells<<" temp = "<<species[s].prim[j].temperature<<endl;
        }
    }
    
    //
    // Emergency BC
    //
    
    for(int s=0; s<num_species; s++)  {
        //species[s].prim[num_cells].temperature   = species[s].prim[num_cells-1].temperature; // species[s].const_T_space;
        //species[s].prim[num_cells+1].temperature = species[s].prim[num_cells-1].temperature; //species[s].const_T_space;
        //species[s].prim[num_cells+2].temperature = species[s].prim[num_cells-1].temperature; //species[s].const_T_space;
    }
    for(int b=0; b<num_bands_out; b++)  {
        //Jrad_FLD(num_cells, b)   = 0.99*Jrad_FLD(num_cells-1, b);
        //Jrad_FLD(num_cells+1,b)  = 0.98*Jrad_FLD(num_cells-1, b);
        //Jrad_FLD(num_cells+2,b)  = 0.97*Jrad_FLD(num_cells-1, b);
    }
    
    
    // Update energies. 
    // TODO: We should add cv * (Tf - Ti) to u to conserve energy properly.
    for(int si=0; si<num_species; si++) {
        species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }

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

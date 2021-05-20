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

        ///////////////////////////////////////////////////////////
        //
        // Chapter on radiative fluxes
        //
        //
        ///////////////////////////////////////////////////////////
        
        //void c_Sim::update_radiation(std::vector<AOS>& u)
void c_Sim::transport_radiation() {
        
        int iters = 0;
        
        if(steps == 1)
            cout<<"   in transport_radiation, before loop, steps = "<<steps<<endl;
        
        //while( ( iters < rad_solver_max_iter) && (epsilon_rad > epsilon_rad_min) ) {
        while( iters < 1) {

            if(steps == 1)
                cout<<"   in transport_radiation, in loop, steps = "<<steps<<endl;
            
            update_opacities(); //Contains solar absorption 
            
            if(debug >= 1)
                cout<<"   in transport_radiation, after update opacities."<<endl;
            
            if(radiation_solver == 0)
                update_fluxes_FLD(); //Contains radiation transport and sourcing
            else 
                update_fluxes_simple();
            
            if(steps == -1 || steps == -2) {
                cout<<"   in transport_radiation, after update fluxes."<<endl;
                cout<<"Jrad_FLD2 and T_FLD2= "<<endl<<Jrad_FLD2<<endl<<endl<<T_FLD2<<endl;
            }
            
            iters++;
        }
}

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

void c_Sim::update_opacities() {
    
    //Compute kappa
    for(int s=0; s<num_species; s++)
        species[s].update_opacities(); //now to be found in opacities.cpp
    
    //Compute dtau = dx * kappa * dens
    
    for(int j = num_cells + 1; j>0; j--) {
        
        for(int s=0; s<num_species; s++)
                species[s].dS(j)  = 0.;
        
        for(int b=0; b<num_bands; b++) {
            
            cell_optical_depth(j,b)           = 0.;
            total_opacity(j,b)                = 0 ;
            total_opacity_twotemp(j,b)        = 0 ;
            radial_optical_depth(j,b)         = 0.;
            radial_optical_depth_twotemp(j,b) = 0.;
            
            for(int s=0; s<num_species; s++) {
                total_opacity(j,b)          += species[s].opacity(j,b)         * species[s].u[j].u1 ;
                total_opacity_twotemp(j,b)  += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 ;
                //cell_optical_depth(j,b)  = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j];
            }
            cell_optical_depth(j,b)         = total_opacity(j, b)         * dx[j] ;
            cell_optical_depth_twotemp(j,b) = total_opacity_twotemp(j, b) * dx[j] ;
            
            if(j==num_cells+1) {
                radial_optical_depth(j,b)         = 0.; //cell_optical_depth(j,b);
                radial_optical_depth_twotemp(j,b) = 0.; //cell_optical_depth(j,b);
            }
            else{
                radial_optical_depth(j,b)         = radial_optical_depth(j+1,b)         + cell_optical_depth(j,b);
                radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + cell_optical_depth_twotemp(j,b);
            }
            
            //
            // After the total optical depth per band is known, we assign the fractional optical depths
            // Maybe merge with previous loop for optimization
            //
            
            for(int s=0; s<num_species; s++) 
                species[s].fraction_total_opacity(j,b) = species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth_twotemp(j,b);
                //species[s].fraction_total_opacity(j,b) = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth(j,b);
            
            //
            // Now compute the attenuation of solar radiation and then assign the lost energy back to individual species in a manner that conserves energy
            //
            
            S_band(j,b) = solar_heating(b) * std::exp(-const_opacity_solar_factor*radial_optical_depth_twotemp(j,b));
            
            //if(steps == 0)
            //    cout<<"IN SOLAR HEATING, const_opacity_solar_factor = "<<const_opacity_solar_factor<<endl;
            
            if(debug >= 1)
                cout<<" in cell ["<<j<<"] band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<" tau_rad(j,b) = "<<radial_optical_depth(j,b)<<" exp(tau) = "<<std::exp(-radial_optical_depth(j,b))<<endl;
            
            if(steps > -1) {
                        
                if(j<num_cells+1)
                    //dS_band(j,b) = (surf[j+1] * S_band(j+1,b) - surf[j] * S_band(j,b))/vol[j]/4.; //       4pi J = c Erad
                    //dS_band(j,b) = (S_band(j+1,b) - S_band(j,b))/dx[j]/4. * (1. - bond_albedo); //       4pi J = c Erad
                    dS_band(j,b) = solar_heating(b) * (-exp(-const_opacity_solar_factor*radial_optical_depth_twotemp(j+1,b)) * expm1( const_opacity_solar_factor* ( radial_optical_depth_twotemp(j+1,b) - radial_optical_depth_twotemp(j,b))) ) /dx[j]/4. * (1.-bond_albedo);
                else
                    // Use optically thin limit  
                    dS_band(j,b) = solar_heating(b)*const_opacity_solar_factor*total_opacity_twotemp(j,b)/4;
                
                dS_band_zero(j,b) = dS_band(j,b);
            }
            else 
                dS_band(j,b) = dS_band_zero(j,b);
            
            
            
            
            //if(const_opacity_solar_factor*radial_optical_depth_twotemp(j,b))
            
            for(int s=0; s<num_species; s++)
                species[s].dS(j)  += dS_band(j,b) * species[s].fraction_total_opacity(j,b);
        }
    }
    
    //Initialize optically thin regions with low energy density
    if((steps==0) && (init_J_factor > 0.)) {
        
        cout<<" Setting J to custom values. init_J_factor = "<<init_J_factor<<" and init_T_temp ="<<init_T_temp<<" const_opacity_solar_factor = "<<const_opacity_solar_factor<<endl;
        
        for(int j = num_cells + 1; j>0; j--) {
            for(int b=0; b<num_bands; b++) {
                
                if(radial_optical_depth(j,b) < 1.) 
                    Jrad_FLD(j,b) *= init_J_factor/pow(x_i[j]/x_i[34] ,2.);
                
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
    
    for(int b = 0; b < num_bands; b++) {
        J_remnant += S_band(1,b);
    }
    
    for(int s=0; s<num_species; s++) 
        if (species[s].const_T_space > 0) {
            species[s].prim[num_cells].temperature   = species[s].const_T_space;
            species[s].prim[num_cells+1].temperature = species[s].const_T_space;
        }
    
    T_rad_remnant = pow(pi * J_remnant / sigma_rad, 0.25);
    
    T_target      = pow(pow(T_core,4.) + pow(T_rad_remnant,4.), 0.25);
    T_surface     = T_target; //+ (T_surface - T_target)*std::exp(-dt/cooling_time);
    
    //Set planetary temperature (computed in update_opacities)
    if(use_planetary_temperature == 1) {
        
        for(int s=0; s<num_species; s++) {
            species[s].prim[1].temperature = T_surface;
            species[s].prim[0].temperature = T_surface;
        }
        for(int b=0; b<num_bands; b++) {
            Jrad_FLD(1, b)            = sigma_rad * pow(T_surface,4.) * pow(R_core,2.) * compute_planck_function_integral3(l_i[b], l_i[b+1], T_surface) ;
            Jrad_FLD(0, b)            = sigma_rad * pow(T_surface,4.) * pow(R_core,2.) * compute_planck_function_integral3(l_i[b], l_i[b+1], T_surface) ;
        }
        //cout<<"t ="<<globalTime<<" Jrad(0,b) = "<<Jrad_FLD(0, 0)<<" Tsurface = "<<T_surface<<" Tcore = "<<T_core<<" T_rad_remnant = "<<T_rad_remnant<<endl;
    }
    
    if(debug >= 1) {
        
        for(int b=0; b<num_bands;b++) {
           cout<<" in band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<endl;
        }
        cout<<"    in update opa "<<endl;
        for(int j=num_cells; j>0; j--) {
            
            for(int b=0; b<num_bands;b++) {
                cout<<"    band ["<<b<<"][j="<<j<<"] = "<<S_band(j,b)<<" dS = "<<dS_band(j,b)<<" tau = "<<const_opacity_solar_factor*radial_optical_depth_twotemp(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fract["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_opacity(j,b))<<" opa = "<<species[s].opacity(j,b)<<" total opa(j+1)"<<total_opacity(j+1,b)<<" rho/T = "<<species[s].prim[j].density<<"/"<<species[s].prim[j].temperature<<" x = "<<x_i12[j];
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


void c_Sim::update_fluxes_FLD() {
    
    // If test, then set J to initial values before computing on
    if(radiation_matter_equilibrium_test == 1) {
        
        for (int j=0; j < num_cells+2; j++) {
            for(int b=0; b<num_bands; b++) {
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

    int num_vars = num_bands + num_species;
    int stride = num_vars * num_vars ;
    int size_r = (num_cells + 2) * num_vars ;
    int size_M = (num_cells + 2) * stride ;

    std::vector<double> 
        l(size_M, 0.), d(size_M, 0.), u(size_M, 0.), r(size_r, 0.) ;
    
    // Step 1: setup transport terms (J)
    for(int b=0; b<num_bands; b++) {
        for (int j=0; j < num_cells+2; j++) {
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
                double tau_inv = 0.5 / (dx * rhokr) ;
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
        //if(use_planetary_temperature == 0) {
            for (int j=0; j < num_ghosts; j++) {
                int idx = j*stride + b*(num_vars + 1) ;
                int idx_r = j*num_vars + b ;
                
                l[idx] = 0 ;
                u[idx] = -d[idx] ;
                r[idx_r] = 0 ; 
            }
        //}
        
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
            
            int Ncell   = num_cells - 2*(num_ghosts - 1) ;
            for (int j=0; j < num_ghosts-1; j++) {
                int i       = Ncell + num_ghosts ;

                int idx = i*stride + b*(num_vars + 1) ;
                int idx_r = i*num_vars + b ;  

                if (j == 0) // FLD flux at the edge of the last cell
                    l[idx] = u[i-1] ;
                else // Free-stream 
                    l[i] = -surf[i-1] / (x_i12[i]-x_i12[i-1]) ;

                d[idx] = -l[idx] + surf[ i ] / (x_i12[i+1]-x_i12[i]) ;
                u[idx] = 0 ;
                r[idx_r] = 0 ;
            }

            int idx = (num_cells+1)*stride + b*(num_vars + 1) ;
            int idx_r = (num_cells+1)*num_vars + b ;  
            l[idx] = -d[idx];
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

    // Step 2: Energy exchange terms
    
    if(radiation_matter_equilibrium_test <= 2) { //radtests 3 and 4 are delta-radiation peaks without energy-matter coupling
        
        for (int j=0; j < num_cells+2; j++) {
            for (int s=0; s < num_species; s++) {
                
                if(debug > 1) {
                        cout<<" Going into radpart 2. ";
                        cout<<" opa    = "<<species[s].opacity_planck(j, 0);
                        cout<<" temper = "<<species[s].prim[j].temperature;
                        cout<<" dens   = "<<species[s].prim[j].density<<endl;
                }
                
                int idx_s  = j*stride + (s + num_bands) * (num_vars+1) ;
                int idx_rs = j*num_vars + (s + num_bands) ;
                
                
                double Ts = species[s].prim[j].temperature ;
                double rhos = species[s].prim[j].density ;

                d[idx_s ] = 1 / dt ;
                r[idx_rs] = Ts / dt ;
                r[idx_rs] += species[s].dS(j) / species[s].u[j].u1 / species[s].cv; //* 4. * pi /c_light
                
                //if (j < num_ghosts || j >= num_cells + 2-num_ghosts) continue ;
                //if (j < num_ghosts ) continue ;

                
                for(int b=0; b<num_bands; b++) {
                    int idx_b  = j*stride + b * (num_vars+1) ;
                    int idx_bs = j*stride + b * num_vars + (s + num_bands) ; 
                    int idx_sb = j*stride + (s + num_bands) * num_vars + b ;
                    int idx_rb = j*num_vars + b ; 
                        
                    double fac = no_rad_trans * species[s].opacity_planck(j, b) * sigma_rad*Ts*Ts*Ts / pi * compute_planck_function_integral3(l_i[b], l_i[b+1], species[s].prim[j].temperature);
                    
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
                    int idx  = j*stride + (si + num_bands) * num_vars + num_bands;
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
    tridiag.solve(&r[0], &r[0]) ; // Solve in place (check it works)

    // Store the result
    for (int j=0; j < num_cells+2; j++) {
        for(int b=0; b<num_bands; b++) {
            Jrad_FLD(j, b) = r[j*num_vars + b] ;
            
            if(radiation_matter_equilibrium_test == 1) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
                    
          //  std::cout << j << "\t" <<  Jrad_FLD(j, b) << "\n" ;
        }

        for(int s=0; s<num_species; s++) {
            species[s].prim[j].temperature = r[j*num_vars + (s + num_bands)] ;
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


void c_Sim::update_fluxes_simple() {
    
    //Compute change in energy
    for (int j=0; j < num_cells+2; j++) {
        for(int s=0; s<num_species; s++)
            species[s].prim[j].temperature += species[s].dS(j) * dt  / species[s].u[j].u1 / species[s].cv ;
    }
    
    
    // Update Temperatures
    for(int si=0; si<num_species; si++) {
            species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
            species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
            species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }
    
    
}

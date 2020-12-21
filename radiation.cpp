///////////////////////////////////////////////////////////
//
//
//  radiation.cpp
//
// This file contains routines doing radiation transport and updating temperatures accordingly.
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
        double epsilon_rad = 1.;
        
        if(steps == 1)
            cout<<"   in transport_radiation, before loop, steps = "<<steps<<endl;
        
        //while( ( iters < rad_solver_max_iter) && (epsilon_rad > epsilon_rad_min) ) {
        while( iters < 1) {

            if(steps == 1)
                cout<<"   in transport_radiation, in loop, steps = "<<steps<<endl;
            
            update_opacities();
            
            if(debug >= 1)
                cout<<"   in transport_radiation, after update opacities."<<endl;
            
            update_fluxes_FLD();
            
            //update_fluxes_FLD(1.0*dt,   Jrad_FLD, Jrad_FLD2, T_FLD); //J_in, J_out, T_in
            //update_temperatures(1.0*dt, Jrad_FLD2, T_FLD, T_FLD2); // J_in,  T_in, T_out
            
            if(steps == -1 || steps == -2) {
                cout<<"   in transport_radiation, after update fluxes."<<endl;
                cout<<"Jrad_FLD2 and T_FLD2= "<<endl<<Jrad_FLD2<<endl<<endl<<T_FLD2<<endl;
            }
            
            iters++;
        }
}



//
// Computation of opacities based on material properties
//

void c_Sim::update_opacities() {
    
    //Compute kappa
    for(int s=0; s<num_species; s++)
        species[s].update_opacities();
    
    //Compute dtau = dx * kappa * dens
    
    for(int j=0; j<num_cells+2; j++) {
        
        for(int s=0; s<num_species; s++)
                species[s].dS(j)  = 0.;
        
        for(int b=0; b<num_bands; b++) {
            
            cell_optical_depth(j,b) = 0.;
            total_opacity(j,b) = 0 ;
            radial_optical_depth(j,b) = 0.;
            
            for(int s=0; s<num_species; s++) {
                total_opacity(j,b)      += species[s].opacity(j,b) * species[s].u[j].u1 ;
                //cell_optical_depth(j,b)  = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j];
            }
            cell_optical_depth(j,b) = total_opacity(j, b) * dx[j] ;
            
            if(j==num_cells)
                radial_optical_depth(j,b) = cell_optical_depth(j,b);
            else
                radial_optical_depth(j,b) = radial_optical_depth(j+1,b) + cell_optical_depth(j,b);
                
            //
            // Afyter the total optical depth per band is known, we assign the fractional optical depths
            // Maybe merge with previous loop for optimization
            //
            
            for(int s=0; s<num_species; s++) 
                    species[s].fraction_total_opacity(j,b) = species[s].opacity(j,b) * species[s].u[j].u1 * dx[j] / cell_optical_depth(j,b);
            
            //
            // Now compute the attentiation of solar radiation and then assign the lost energy back to individual species in a manner that conserves energy
            //
            
            S_band(j,b)  = solar_heating(b) * std::exp(-radial_optical_depth(j,b));
            
            if(debug >= 1)
                cout<<" in cell ["<<j<<"] band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<" tau_rad(j,b) = "<<radial_optical_depth(j,b)<<" exp(tau) = "<<std::exp(-radial_optical_depth(j,b))<<endl;
            
            if(j<num_cells)
                dS_band(j,b) = (surf[j+1] * S_band(j+1,b) - surf[j] * S_band(j,b))/vol[j];
            else
                dS_band(j,b) = 0;
            
            for(int s=0; s<num_species; s++)
                species[s].dS(j)  += dS_band(j,b) * species[s].fraction_total_opacity(j,b);
        }
    }
    
    if(debug >= 1) {
        
        for(int b=0; b<num_bands;b++) {
           cout<<" in band ["<<b<<"] top-of-the-atmosphere heating = "<<solar_heating(b)<<endl;
        }
        
        for(int j=num_cells; j>0; j--) {
            
            cout<<"    in update opa "<<endl;
            
            for(int b=0; b<num_bands;b++) {
                cout<<"    band ["<<b<<"][j="<<j<<"] = "<<S_band(j,b)<<" banddS = "<<dS_band(j,b)<<" tau = "<<radial_optical_depth(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fraction["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_opacity(j,b))<<" opa = "<<species[s].opacity(j,b)<<" total opa(j+1)"<<total_opacity(j+1,b)<<" x = "<<x_i12[j];
                }
                cout<<endl;
            }
            cout<<endl;
        }
        
        if(debug > 3) {
            char stop;
            cin>>stop;    
        }
        
    }
}

void c_Species::update_opacities() {
    
    for(int j=0; j< num_cells+2; j++) {
        for(int b=0; b<num_bands; b++) {
            
            if(base->radiation_matter_equilibrium_test <= 3) {
                opacity(j,b)        = const_opacity; // TODO: Replace with some_complex_tabulated_opacity_function() or link with the opacity_data matrix
                opacity_planck(j,b) = const_opacity; // TODO: Replace with some_complex_tabulated_opacity_function() or link with the opacity_data matrix   
            }
            else if(base->radiation_matter_equilibrium_test == 4){ //The nonlinear commercon radiative diffusion test
                opacity(j,b)        = 1e11/this->u[j].u1*pow(base->Jrad_FLD(j,b)/c_light*4.*pi, -1.5 );
                opacity_planck(j,b) = 1e11/this->u[j].u1*pow(base->Jrad_FLD(j,b)/c_light*4.*pi, -1.5 );                
            }
            
            
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
                double tau_inv = 0.5 / (dx * (total_opacity(j,b) + total_opacity(j+1,b))) ;
                double R       = 2 * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300) ;
                double D       = 1. * surf[j] * flux_limiter(R) * tau_inv;

                // divergence terms
                u[idx] = -D ;
                d[idx] += D ;
                d[idx+stride] = D ;
                l[idx+stride] = -D ;
                
                if(debug > 1)
                    cout<<" radiation part 0. t,j,b="<<steps<<","<<j<<","<<b<<" tau_inv/R/D = "<<tau_inv<<"/"<<R<<"/"<<D<<" J/J = "<<Jrad_FLD(j+1,b)<<"/"<<Jrad_FLD(j,b)<<endl;
            }
            

        }
        
        // Boundaries:
        
        for (int j=0l; j < num_ghosts; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;
            // Left boundary:
            //    Reflecting / no flux
            l[idx] = 0 ;
            d[idx] = +1 ;
            u[idx] = -1 ;
            r[idx_r] = 0 ;
        }
        
        //   Right boundary: reflective?
        //if(geometry == Geometry::cartesian) {
        if(false) {

            int Ncell = num_cells - 2*(num_ghosts - 1);
            for (int j=0; j < num_ghosts; j++) {
                int i = Ncell + num_ghosts + j ;

                int idx = i*stride + b*(num_vars + 1) ;
                int idx_r = i*num_vars + b ;    

                l[idx] = -1;
                d[idx] = +1 ;
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
        
        for (int j=0; j < num_cells+2; j++)
            for (int s=0; s < num_species; s++) {
                
                int idx_s  = j*stride + (s + num_bands) * (num_vars+1) ;
                int idx_rs = j*num_vars + (s + num_bands) ;

                double Ts = species[s].prim[j].temperature ;
                double rhos = species[s].prim[j].density ;

                d[idx_s ] = 1 / dt ;
                r[idx_rs] = Ts / dt ;

                if (j < num_ghosts || j >= num_cells + 2-num_ghosts) continue ;

                for(int b=0; b<num_bands; b++) {
                    int idx_b  = j*stride + b * (num_vars+1) ;
                    int idx_bs = j*stride + b * num_vars + (s + num_bands) ; 
                    int idx_sb = j*stride + (s + num_bands) * num_vars + b ;
                    int idx_rb = j*num_vars + b ; 

                    // TODO multiply by fraction in band.
                    double fac = species[s].opacity(j, b) * sigma_rad*Ts*Ts*Ts / pi; // * compute_planck_function_integral3(l_i[b], l_i[b+1], species[s].prim[j].temperature);
                    
                    d[idx_s ] += 16 * pi * fac / species[s].cv ;
                    d[idx_sb] = - 4 * pi * species[s].opacity(j, b) / species[s].cv ;
                    r[idx_rs] += 12 * pi * fac * Ts / species[s].cv ;
                    r[idx_rs] += species[s].dS(j)  / species[s].u[j].u1 / species[s].cv; //* 4. * pi /c_light

                    //cout<<" in rad matrix["<<j<<"]: term1 = "<<12 * pi * fac * Ts / species[s].cv<<" term2 = "<<species[s].dS(j)  / species[s].u[j].u1 / species[s].cv<<endl;
                    
                    d[idx_b ] += vol[j] * rhos * species[s].opacity(j, b) ;
                    d[idx_bs] = - 4 * vol[j] * rhos * fac ;
                    r[idx_rb] -=  3 * vol[j] * rhos * fac * Ts ;
                    
                    if(debug > 1)
                        cout<<" radiation part2, t = "<<steps<<" b["<<b<<"] i["<<j<<"] s["<<s<<"] l/d/u/r_rs/rs_rb = "<<d[idx_s]<<"/"<<d[idx_sb]<<"/"<<d[idx_b]<<"/"<<d[idx_bs]<<"/"<<r[idx_rs]<<"/"<<r[idx_rb]<<" opac = "<<vol[j] * rhos * species[s].opacity(j, b)<<endl;
                }
            }
            
    }
    if(debug > 4) {
        char stepstop;
        cin>>stepstop;
    }
    
    // TODO:
    //   Add collision terms here

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
       //     std::cout <<  "\t" << species[s].prim[j].temperature << "\n" ;
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

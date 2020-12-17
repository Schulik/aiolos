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
            
            cout<<"    in update opa ";
            
            for(int b=0; b<num_bands;b++) {
                cout<<"    band ["<<b<<"] = "<<S_band(j,b)<<" banddS = "<<dS_band(j,b)<<" tau = "<<radial_optical_depth(j,b);
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fraction["<<species[s].speciesname<<"] = "<<(species[s].fraction_total_opacity(j,b))<<" opa = "<<species[s].opacity(j,b)<<" total opa(j+1)"<<total_opacity(j+1,b)<<" x = "<<x_i12[j];
                }
                cout<<endl;
            }
            cout<<endl;
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
                double R       = 2 * tau_inv * std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b)) / (Jrad_FLD(j+1,b) + Jrad_FLD(j, b)) ;
                double D       = 1. * surf[j] * flux_limiter(R) * tau_inv;

                // divergence terms
                u[idx] = -D ;
                d[idx] += D ;
                d[idx+stride] = D ;
                l[idx+stride] = -D ;
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
                
                cout<<" timestep "<<steps<<" band["<<b<<"] cell["<<index<<"] l/d/u/r = "<<l[index]<<"/"<<d[index]<<"/"<<u[index]<<"/"<<r[index];
                cout<<" temps = ";
                for(int si = 0; si<num_species; si++) {
                        cout<<species[si].prim[2].temperature<<" ";
                }
                cout<<endl;
            }
            
            if(steps <= 2) {
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
                }
            }
            
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

void c_Sim::fill_rad_basis_arrays(int j, double timestep, Eigen::MatrixXd &J_in, Eigen::MatrixXd &T_in ) { //Called in compute_radiation() in source.cpp
    
        if(debug >= 2) cout<<" IN FILL RAD BASIS, pos0.;"<<endl;
        
        for(int si=0; si<num_species; si++) {
            double J_bandsum = 0.;
            
            dens_vector(si)        = species[si].u[j].u1; 
            //numdens_vector(si)     = species[si].prim[j].number_density; 
            //mass_vector(si)        = species[si].mass_amu;
            
            if(debug >= 2) cout<<" IN FILL RAD BASIS, species = "<<si<<" j ="<<j<<" x[j]="<<x_i[j]<<endl;
            
            //radiation_vec_input(si)= dt/species[si].cv*(12.*species[si].opacity_planck(j)*sigma_rad*pow(species[si].prim[j].temperature,4.) + species[si].dS(j)/dens_vector(si) + c_light * Erad_FLD_total(j) ); //TODO: kappa* c*Erad needs to be a band sum! -> c sum_b kappa_b Erad_b
            
            for(int b=0;b<num_bands; b++)
                J_bandsum += species[si].opacity(j,b) * J_in(j,b);
                //J_bandsum += species[si].opacity(j,b) * Jrad_FLD(j,b);
            
            radiation_T3_vector(si)= timestep/species[si].cv * 16 * species[si].opacity_planck(j) * sigma_rad * pow(T_in(j,si),3.);
            
            radiation_vec_input(si) = T_in(j,si) * (1. + timestep/species[si].cv*12.*species[si].opacity_planck(j)*sigma_rad*pow(T_in(j,si),3.)  ) ;
            radiation_vec_input(si) += timestep/species[si].cv * ( species[si].dS(j)/dens_vector(si) + 4.*pi*J_bandsum ); 
            
            if(debug >= 2) cout<<" IN FILL RAD BASIS, rhs, T3-term debug, T = "<<species[si].prim[j].temperature<< " dt/cv = "<<timestep/species[si].cv<<" kappa_planck = "<<species[si].opacity_planck(j)<<" sT3 = "<<sigma_rad*pow(T_in(j,si),3.)<<" T = "<<T_in(j,si)<<endl;
            
            if(debug >= 2) cout<<" IN FILL RAD BASIS, rhs, T3-term = "<<(timestep/species[si].cv*(12.*species[si].opacity_planck(j)*sigma_rad*pow(T_in(j,si),3.) ))<<" solar heating = "<<(dt/species[si].cv * ( species[si].dS(j)/dens_vector(si)))<<" cooling "<<(timestep/species[si].cv*4.*J_bandsum )<<endl;
            
            radiation_cv_vector(si)= timestep/species[si].cv;
            
            //TODO:Specify S properly and band_integral over Erad_FLD
        }
        
        if(debug >= 2) cout<<" IN FILL RAD BASIS, posLast.;"<<endl;
}

//
// Implicit energy update based on fluxes
//
void c_Sim::update_temperatures(double timestep, Eigen::MatrixXd &J_in, Eigen::MatrixXd &T_in, Eigen::MatrixXd &T_out) {
    
    Eigen::internal::set_is_malloc_allowed(false) ;
         
    if(debug >= 1) cout<<" IN UPDATE TEMPERATURES, dense friction, num_species = "<<num_species<<endl;
            
    double beta_local;
    double coll_b;

    for(int j=0; j < num_cells+2; j++){

        fill_rad_basis_arrays(j, timestep, J_in, T_in);
        compute_alpha_matrix(j, 1);
        
        // 
        // Loop for determining friction coefficients is not needed here, as those coeffs have already been computed in friction();
        //
        
        if(debug >= 2) cout<<"    RHS = "<<endl<<radiation_vec_input<<endl;
        
        radiation_matrix_M.noalias() = radiation_cv_vector.asDiagonal() * friction_coefficients;
        
        radiation_matrix_T = identity_matrix;
        radiation_matrix_T.diagonal().noalias() += radiation_T3_vector;
        radiation_matrix_T.diagonal().noalias() +=  (friction_coefficients * unity_vector) * radiation_cv_vector;
        //radiation_matrix_T.noalias() += radiation_matrix_M;
        
        LU.compute(radiation_matrix_T);
        radiation_vec_output.noalias() = LU.solve(radiation_vec_input);
        
        //if(debug >= 1) cout<<"    T_matrix = "<<endl<<radiation_matrix_T<<endl;
        //if(debug >= 1 && j==7 && steps == 1) cout<<"    T = "<<radiation_matrix_T<<endl;
        //if(debug >= 1 && j==7 && steps == 1) cout<<"    inp_vector = "<<radiation_vec_input<<endl;
        //if(debug >= 1 && j==7 && steps == 1) cout<<"    T_out = "<<radiation_vec_output<<endl;
        
        if(debug >= 1 && j==num_cells) {
            cout<<"    T_init = "<<endl;
            for(int si=0; si<num_species; si++)
                cout<<species[si].prim[j].temperature<<" "<<endl;
            
            cout<<endl;
            cout<<"    T3_vector = "<<endl<<radiation_T3_vector<<endl<<endl;
            cout<<"    cv_vector = "<<endl<<radiation_cv_vector<<endl<<endl;
            cout<<"    T_out = "<<endl<<radiation_vec_output<<endl<<endl;
            cout<<"Radiation matrix: "<<endl<<radiation_matrix_T<<endl<<endl;
            cout<<"alpha coefficients matrix: "<<endl<<friction_coefficients<<endl<<endl;
            
            //cout<<"    T_out = "<<radiation_vec_output;
            //cout<<"    inp_vector = "<<radiation_vec_input<<endl;
            
        }
        
        //
        // Update new temperature and internal energy
        //
        //if(debug >= 1) cout<<"UPDATE TEMPERATURES:"<<endl<<" T before: ";
        //for(int si=0;si<num_species; si++) {
        //    cout<<" s["<<species[si].name<<"].T = "<<species[si].prim[j].temperature;
        //}
        //cout<<endl;
        
        
        for(int si=0; si<num_species; si++)
            T_out(j,si) = radiation_vec_output(si);
        
    }
    /*
    for (int j=0; j < num_cells+1; j++) {
            if(globalTime > 1.) {
                        cout<<" time = "<<globalTime<<" netflux in cell "<<j<<" = "<<(sigma_rad * pow(T_in(j,si),4.)-J_in(j,0))<<" r = "<<r[j]<<endl;
                    
            }
        }*/
    
    if(debug >= 1)
    {            
        char chr;
        cin>>chr;
    }
    
    for(int si=0; si<num_species; si++) {
            species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
            species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
            species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
        }

}


//
// Computation of fluxes based on radiative transport theory
//
void c_Sim::update_fluxes(double timestep) {
/*    
    for(int b=0; b<num_bands; b++) {
        F_up(0,b) = 0.; //Some value that makes sense with T_internal
        F_down(num_cells,b) = 0.; //Some value that makes sense with T_irradiated + Nonthermal irradiation
    }
    
    for(int j=1; j<=num_cells; j++) {
        for(int b=0; b<num_bands; b++) {
            
            F_up(j,b) = F_up(j-1,b) - cell_optical_depth(j,b) *1.; 
            //TODO: Replace 1. with appropriate geometric factors for (non)-plane parallel emission
            //TODO: Add appropriate term for source function
        }
    }
    
    
    for(int j=num_cells-1; j<=0; j--) {
        for(int b=0; b<num_bands; b++) {
            
            F_down(j,b) = F_up(j+1,b) - cell_optical_depth(j,b) *1.; 
            //T ODO: Replace 1. with appropriate geometric factors for (non)-plane parallel emission
            //TODO: Add appropriate term for source function
            
        }
        
    }
    
    F_plus  = F_up + F_down;
    F_minus = F_up - F_down;*/
}

void c_Sim::update_fluxes_FLD2(double timestep, Eigen::MatrixXd &J_in, Eigen::MatrixXd &J_target, Eigen::MatrixXd &T_in) {
/*
    auto flux_limiter = [](double R) {
        if (R <= 2)
            return 2 / (3 + std::sqrt(9 + 10*R*R)) ;
        else 
            return 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
    } ;

    std::vector<double> l(num_cells+2), d(num_cells+2), u(num_cells+2), r(num_cells+2) ;
    
    //double saveR, saveL, savetau, saveJ1, saveJ2;
    
    for(int b=0; b<num_bands; b++) {
   
        for (int j=0; j < num_cells+1; j++) {
            double dx      = (x_i12[j+1]-x_i12[j]) ;
            double tau_inv = 0.5 / (dx * (total_opacity(j,b) + total_opacity(j+1,b))) ;
            double R       = 2 * tau_inv * std::abs(J_in(j+1,b) - J_in(j,b)) / (J_in(j+1,b) + J_in(j, b)) ;
            //safeguard against 0/0 at low intensities
            if((J_in(j+1,b) + J_in(j, b)) < 1e-20)
                R = 2 * tau_inv * std::abs(J_in(j+1,b) - J_in(j,b)) / (2e-20) ;
            double D       =c_light* timestep * surf[j] * flux_limiter(R) * tau_inv;
            
            // divergence terms
            u[j] = -D ;
            d[j] += D ;
            d[j+1] = D ;
            l[j+1] = -D ;

            // source terms
            d[j] +=c_light* timestep * vol[j] * total_opacity(j,b) ;
            
            //Thermal emission field from all species in this band
            r[j] = 0;
            if (num_bands == 1) {
                for(int s=0; s<num_species; s++) {
                    r[j] +=c_light* timestep * vol[j] * species[s].prim[j].density * species[s].opacity_planck(j) * 
                        sigma_rad*pow(T_in(j,s),4) / pi;
                }
            } else {
                for(int s=0; s<num_species; s++) {
                    r[j] +=c_light* timestep * vol[j] * species[s].prim[j].density * species[s].opacity(j,b) * 
                        compute_planck_function_integral(l_i[b], l_i[b+1], T_in(j,s)) ;
                }
            }
            
            // Time dependent terms:
            d[j] += vol[j]  ;
            r[j] += vol[j] * J_in(j, b) ;
        }
        
        // Boundaries:
        
        for (int j=0; j < num_ghosts; j++) {
            // Left boundary:
            //    Reflecting / no flux
            l[j] = r[j] = 0 ;
            d[j] = +1 *c_light* timestep ;
            u[j] = -1 *c_light* timestep;

        }
        
        //   Right boundary: reflective?
        //if(geometry == Geometry::cartesian) {
        if(true) {

            int Ncell = num_cells - 2*(num_ghosts - 1);
            for (int j=0; j < num_ghosts; j++) {
                int i = Ncell + num_ghosts + j ;
                
                l[i] = -1*c_light* timestep;
                d[i] = +1*c_light* timestep ;
                u[i] = 0 ;
                r[i] = 0 ;
            }
        }
        else {//   Right boundary: free stream, no emission / absorbtion.
            
            int Ncell   = num_cells - 2*(num_ghosts - 1) ;
            for (int j=0; j < num_ghosts-1; j++) {
                int i       = Ncell + num_ghosts ;
            
                if (j == 0) // FLD flux at the edge of the last cell
                    l[i] = u[i-1] * c_light * timestep ;
                else // Free-stream 
                    l[i] = -surf[i-1] / (x_i12[i]-x_i12[i-1])  * c_light * timestep;

                d[i] = -l[i] + surf[ i ] / (x_i12[i+1]-x_i12[i]) * c_light * timestep ;
                u[i] = r[i] = 0 ;
            }
            l[num_cells+1] = -1 * c_light*timestep;
            d[num_cells+1] =  1 * c_light*timestep;
            u[num_cells+1] = 0;
            r[num_cells+1] = 0;
        }
        if(debug >= 1) {
            cout<<" band["<<b<<"] l/d/u/r = "<<l[2]<<"/"<<d[2]<<"/"<<u[2]<<"/"<<r[2];
            cout<<" l+d+r = "<<l[2]+d[2]+r[2];
            //<<" l("<<saveR<<") = "<<saveL;
            //cout<<" tau_inv = "<<savetau<<" J1 = "<<saveJ1<<" J2 = "<<saveJ2<<" opa = "<<total_opacity(0,0);
            cout<<" T_in = ";
            for(int si = 0; si<num_species; si++) {
                    cout<<T_in(2,si)<<" ";
            }
            
            double invcdt = 1./c_light/timestep;
            double f1 = c_light*timestep * species[0].prim[2].density * species[0].opacity_planck(2);
            double f2 = f1 * species[0].opacity_planck(2)*sigma_rad*pow(T_in(2,0),4) / pi;
            double result = (J_in(2,0) + f2)/(1.+f1);
            
            cout<<"    in update flux2 J_old = "<<J_in(2,0)*invcdt<<" cdtkrho = "<<(1.+f1)<<" cdtkrhosigmaT^4 = "<<(f2)<< " result = "<<result;
            
            for(int b=0; b<num_bands;b++) {
                cout<<"    band ["<<b<<"] ";
                for(int s=0; s<num_species; s++) {
                     cout<<" dS_fraction["<<species[s].name<<"] = "<<(species[s].fraction_total_opacity(2,b))<<" opa = "<<total_opacity(2,b);
                }
                cout<<endl;
            }
            
        }
        
        tridiag.factor_matrix(&l[0], &d[0], &u[0]) ;
        tridiag.solve(&r[0], &r[0]) ; // Solve in place (check it works)
        
        if(debug >= 1) {
                cout<< " r_final = "<<r[2]<<endl;
        }
        
        // Store result
        for (int j=0; j < num_cells+2; j++) 
            J_target(j, b) = r[j];
            //Jrad_FLD(j, b) = r[j];
    }*/
}



void c_Sim::correct_totalenergy(double, Eigen::MatrixXd &Etot) {
    
    
    //
    // Get initial energy
    //
    double etot_now = 0.;
    
    if(steps == 0) {
        
        for(int i = 0; i < num_cells+2; i++) {
            
            Etot(i) = 0; 
            
            for(int b=0; b<num_bands; b++) {
                Etot(i) += 0.*Jrad_FLD(i,b)*4*pi/c_light; 
            }
            
            for(int s=0; s<num_species; s++) {
                Etot(i) += 0.*species[s].prim[i].density*species[s].cv*species[s].prim[i].temperature;
            }
        }
        
        
    }
    
    //
    // Use analytic formula to correct total energy
    //
    else{
        
        for(int i = 0; i < num_cells+2; i++) {
            double A = dt * c_light * species[0].opacity(i,0) * species[0].prim[i].density;
            double Q = dt / species[0].cv * species[0].opacity(i,0) * sigma_rad * pow(species[0].prim[i].temperature,3.);
            double J_minus_B = +(Jrad_FLD(i,0)*4*pi/c_light - species[0].prim[i].density*species[0].cv*species[0].prim[i].temperature);
            
            double correction_term = 64.*pi*A*Q/c_light/(1.+16.*Q)/(1.+A)*J_minus_B;
            
            if(correction_term > 0)
                Etot(i) += 1e4*c_light*correction_term;
            else
                Etot(i) += 0.;
        }
    }
}

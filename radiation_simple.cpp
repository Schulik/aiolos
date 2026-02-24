
/**
 * radiation_simple.cpp
 * 
 * This file contains routines computing the transport of radiation in the one outgoing-band approximation. 
 * This allows to decouple solving for the radiation field and the species' temperatures. Great for debugging more complex problems involving thermal radiation.
 */

#define EIGEN_RUNTIME_NO_MALLOC

#include <cassert>
#include "aiolos.h"

/**
 * Solves first for J at the advanced timestep, i.e. J^{n+1} and then plugs that into the formula for all the T_s^{n+1}, which is solved independently in another timestep.
 * Analoguous to Bitsch+2013, in their Appendix to the radiation solver.
 * 
 * @param[in] ddt timestep
 */
void c_Sim::update_fluxes_FLD_simple(double ddt) {

    for(int si=0; si<num_species; si++) {
        for(int j =0; j<= num_cells; j++) {
            AOS      tmp  = species[si].u[j];
            AOS_prim tmpp = species[si].prim[j];

            double tt = tmp.u3 - 0.5 * tmp.u2*tmp.u2/tmp.u1;
            tt /= (tmp.u1*species[si].cv);

            if(std::isnan(tmpp.temperature) || std::isnan(tt) )
                cout<<" before T_solver "<<steps<<" "<<j<<" "<<si<<" "<<tt<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tmpp.temperature<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;                
        }
    }
    //
    // Initial checks
    //
    
    if(debug > 2)
        cout<<"Starting update_fluxes_FLD_simple.."<<endl;
   
    auto flux_limiter = [](double R) {
        if (R <= 2)
            return 2 / (3 + std::sqrt(9 + 10*R*R)) ;
        else 
            return 10 / (10*R + 9 + std::sqrt(81 + 180*R));
    };     
    int num_vars = num_bands_out; // + num_species
    int stride = num_vars * num_vars ;
    int size_r = (num_cells + 2) * num_vars ;
    int size_M = (num_cells + 2) * stride ;
    int size_etas = (num_cells + 2) * num_species;

    std::vector<double> 
        l(size_M, 0.), d(size_M, 0.), u(size_M, 0.), r(size_r, 0.), eta1(size_etas, 0.), eta2(size_etas, 0.), denoms(size_etas, 0.) ;
     std::vector<double> arr_rhokr(size_r, 0.), arr_R(size_r, 0.), arr_D(size_r, 0.), arr_JDJ(size_r, 0.);
     std::vector<double> exchange_d_sums(num_cells+2, 0.), exchange_r_sums(num_cells+2, 0.);

    //std::fill(lhs_sc.begin(), lhs_sc.end(), 0); //zero out Subcycling arrays
    rhs_sc = Vector_t::Zero( (num_cells+2)*num_species );//  std::fill(rhs_sc.begin(), rhs_sc.end(), 0);
    //std::fill(denoms_sc.begin(), denoms_sc.end(), 0);
    
    int numcells_offset = 2; //Nominally 1
    // Step 1: setup transport terms (J)
    for(int b=0; b<num_bands_out; b++) {
        for (int j=0; j < num_cells + numcells_offset; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;

            // Time dependent terms:
            d[idx] +=  vol[j] / (c_light * ddt) ;
            r[idx_r] += (vol[j] / (c_light * ddt)) * Jrad_FLD(j, b) ;

            // Flux across right boundary
            if (j > 0 && j < num_cells + numcells_offset-1) {
                double dx      = (x_i12[j+1]-x_i12[j]) ;                
                double rhokr   = 0;
                
                rhokr = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx );
                rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);
                
                double tau_inv = 1. / (dx * rhokr) ;
                double dJ = std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b))/(Jrad_FLD(j, b) + 1e-300);
                    
                double R       = 1.* xi_rad * tau_inv *  dJ ; // Put in 1.0 as prefactor to get correct rad shock
                double D       = 1.* tau_inv * surf[j] * flux_limiter(R) * no_rad_trans;
                
                arr_rhokr[j] = rhokr;
                arr_R[j]        = R;
                arr_D[j]        = D;
                arr_JDJ[j]      =  std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b));///(4.*pi*(Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300));
                
                // divergence terms
                u[idx] = -D ;
                d[idx] += D ;
                d[idx+stride] = D ;
                l[idx+stride] = -D ;
                
                if(debug > 33)
                    cout<<" radiation part 0. t,j,b="<<steps<<","<<j<<","<<b<<" tau_inv/R/D = "<<tau_inv<<"/"<<R<<"/"<<D<<" J/J/dJ = "<<Jrad_FLD(j+1,b)<<"/"<<Jrad_FLD(j,b)<<"/"<<(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<" flux = "<<D*(Jrad_FLD(j+1,b)-Jrad_FLD(j,b))<<endl;
            }
        }
        
        
        // Boundaries:
        // Left boundary:
        //    Reflecting / no flux or planetary temperature
            for (int j=0; j < num_ghosts; j++) {
                int idx = j*stride + b*(num_vars + 1) ;
                int idx_r = j*num_vars + b ;
                
                l[idx] = 0 ;
                u[idx] = -d[idx] ;
                r[idx_r] = 0 ; 
            }
        
        //   Right boundary: reflective?
        if(closed_radiative_boundaries) {
        

            int Ncell = num_cells - 2*(num_ghosts - 1);
            for (int j=0; j < num_ghosts; j++) {
                int i = Ncell + num_ghosts + j ;

                int idx = i*stride + b*(num_vars + 1) ;
                int idx_r = i*num_vars + b ;    
                
                l[idx] = -l[idx] ;
                u[idx] = 0 ;
                r[idx_r] = 0 ;
            }
        }
        else {//   Right boundary: free stream, no emission / absorbtion.
            
            //  Assume F = J and \div(F) = const

            int idx = (num_cells)*stride ;
            int idx_r = (num_cells)*num_vars;  
            
            double f = 1./xi_rad * x_i12[num_cells]/x_i12[num_cells+1] ;
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
        
        if(debug >= 33) {
            
            for(int index=0; index < num_cells+2; index++) {    
                //                 /int index = (num_cells/2+1);
                
                cout<<" radiation part1, t = "<<steps<<" band["<<b<<"] cell["<<index<<"] l/d/u/r = "<<l[index]<<"/"<<d[index]<<"/"<<u[index]<<"/"<<r[index];
                cout<<" temps = ";
                for(int si = 0; si<num_species; si++) {
                        cout<<species[si].prim[2].temperature<<" ";
                }
                //cout<<endl;
            }
            
            if(debug > 33) {
                char a;
                cin>>a;
                
            }
        }
    }

    // Step 2: Energy exchange terms kappa*rho*(J-B) + dS + Pi + Lambda
    
    if(radiation_matter_equilibrium_test <= 2) { //radtests 3 and 4 are delta-radiation peaks without energy-matter coupling
        
        for (int j=0; j < num_cells+numcells_offset; j++) {
            
            double exchange_d_sum = 0.;
            double exchange_r_sum = 0.;

            //Compute etas
            for (int s=0; s < num_species; s++) {
                
                int idx_s = j * (num_species) + s;
                double Ts = species[s].prim[j].temperature ;
                double Ts3 = Ts*Ts*Ts;
                double rhos = species[s].prim[j].density ;
                double kappa = species[s].opacity_planck(j, 0);

                double fac = 1. * ddt * no_rad_trans * kappa / species[s].cv * sigma_rad * Ts3;
                double moredenom = - photocooling_expansion * species[s].dGdT(j) * ddt / ( species[s].cv * species[s].u[j].u1);
                double denom = 1. + 16. * fac + moredenom ;
                
                double tempeta = 0;
                
                //Non-subcycling heating function
                denoms[idx_s] = denom;
                eta1[idx_s] += Ts * ( 1. + 12. * fac);
                eta1[idx_s] += 1. * ddt * (species[s].dS(j) - species[s].dG(j) - photocooling_expansion * species[s].dGdT(j)*Ts     ) / species[s].u[j].u1 / species[s].cv;
                eta2[idx_s] += 4.* pi * ddt * kappa * no_rad_trans / species[s].cv;

                //Sub-cycling heating function
                if(couple_J_into_T)
                        rhs_sc(idx_s) +=  4.* pi * kappa * Jrad_FLD(j, 0);
                
                exchange_d_sum += no_rad_trans * rhos * kappa * (1 - 4.*sigma_rad*Ts3/pi * eta2[idx_s]/denoms[idx_s]);
                exchange_r_sum += no_rad_trans * rhos * kappa * sigma_rad * Ts3/pi * (4 * eta1[idx_s]/denoms[idx_s] - 3 * Ts );
            }
            
            int idx_b  = j*stride;
            int idx_rb = j*num_vars;
            
            d[idx_b]  += vol[j] * exchange_d_sum;
            r[idx_rb] += vol[j] * exchange_r_sum;
            
            exchange_d_sums[j] += exchange_d_sum;
            exchange_r_sums[j] += exchange_r_sum;
            //Sum up etas to get J source terms
            
        }
    }
    
    if(debug >= 33) {
        
        cout<<"L ="<<endl;
        for(int i = 0; i < size_M; i++) {
            
            cout<<l.at(i)<<" ";
            
        }
        
        cout<<"D ="<<endl;
        for(int i = 0; i < size_M; i++) {
            cout<<d.at(i)<<" ";
            
        }
        
        cout<<"u ="<<endl;
        for(int i = 0; i < size_M; i++) {
  
            cout<<u.at(i)<<" ";
            
        }
        
        char stepstop;
        cin>>stepstop;
    }
    
    //
    // Solve!
    //
    tridiag.factor_matrix(&l[0], &d[0], &u[0]) ;
    tridiag.solve(&r[0], &r[0]) ; // Solve in place
    
    //
    // Check J for negative values and store the result
    //
    
    int Jswitch = 0;
    if( globalTime > 1e15)
        Jswitch = 1;
    
    for (int j=0; j <= num_cells+numcells_offset-1; j++) {
        for(int b=0; b<num_bands_out; b++) {
            if(solve_for_j)
                Jrad_FLD(j, b) = r[j*num_vars + b] ;
                        
            //Check for negative J
            if(Jrad_FLD(j, b) < 0. && couple_J_into_T) {
                cout<<" -J in j/steps "<<j<<"/"<<steps<<" rhokr = "<< arr_rhokr[j]<<" R = "<< arr_R[j]<<" D ="<< arr_D[j]<<" dJ = "<<arr_JDJ[j]<<" J = "<<Jrad_FLD(j, b)<<" Ji+Ji+1 = "<<Jrad_FLD(j+1, b)+Jrad_FLD(j, b)<<" exchange sums = "<<exchange_d_sums[j]<<"/"<<exchange_r_sums[j]<<endl;
                
                Jswitch = 1;
            }
            
            if(radiation_matter_equilibrium_test == 1) {
                Jrad_FLD(j, b) = Jrad_init(j,b);
            }
        }
    }
    
    //Stop for a check with the user, if we find negative J?
     if(Jswitch == 1) {
        
        //char a;
        //cin>>a;
    }
    
    //
    // Check T for negative values and store
    //
    
    int Tswitch = 0;
    
    //
    // Compute Ti-Tj terms via a separate coupling step
    //
    
    if (use_collisional_heating && num_species > 1) {
        //compute_collisional_heat_exchange();
        
        for (int j=0; j < num_cells+numcells_offset; j++) {
            double init_tmean  = return_T_mean(j);
            double init_etotal = return_e_total(j);

            if(1==0) { //Old heat solver

                fill_alpha_basis_arrays(j);
                compute_alpha_matrix(j);
                do_highenergy_cooling(j, std::max(init_tmean, 3.) );
                
                coll_heat_matrix.setZero();
                coll_heat_b.setZero();
                coll_heat_output.setZero();

                double   tau = total_opacity(j,0) * (x_i12[j+1]-x_i12[j]);
                double security_multiplier = tau<1e3?1.:1.e-4;
                //compute_collisional_heat_exchange_matrix(j);  21 May 2023: This function has been disabled in radiation_simple, the coefficients are computed in place there now

                for(int si=0; si<num_species; si++) {
                    
                    int idx_s = j * (num_species) + si;
                    
                    //coll_heat_matrix(si,si) += 1./ddt ;
                    //coll_heat_b(si)         += species[si].prim[j].temperature / ddt;
                    //
                    // Rad equilibrium terms
                    //
                    coll_heat_matrix(si,si) += denoms[idx_s];
                    if(couple_J_into_T)
                        coll_heat_b(si)     += eta2[idx_s]*Jrad_FLD(j,0);
                    coll_heat_b(si)         += eta1[idx_s];
                    
                    double diag_sum = 0;
                    double temp = 0;
                    
                    for(int sj=0; sj<num_species; sj++) {
                        temp = ddt * 1. * friction_coefficients(si,sj) * 3 * kb / (species[si].cv * (mass_vector(si) + mass_vector(sj)) );
                        diag_sum += temp;
                        coll_heat_matrix(si,sj) -= temp;
                    }

                    coll_heat_matrix(si,si) += diag_sum;
                }

                //double e_orig = return_total_e(j);

                LU.compute(coll_heat_matrix) ;
                coll_heat_output.noalias() = LU.solve(coll_heat_b);

                double final_tmean  = return_T_mean(j, coll_heat_output);
                double final_etotal = return_e_total(j, coll_heat_output);

                if((j==10) && steps >= 457) {
                    cout<<" Tfin = ";
                    for(int ss=0; ss<num_species; ss++){ cout<<" "<<coll_heat_output(ss); }
                    cout<<endl; 
                    cout<<" b = ";
                    for(int ss=0; ss<num_species; ss++){ cout<<" "<<coll_heat_b(ss); }
                    cout<<endl; 
                    cout<<" mat = "<<endl;
                    cout<<coll_heat_matrix; 
                    cout<<endl;
                    cout<<"steps "<<steps<<" Reporting change in mean temperature and internal energy: "<<(1-init_tmean/final_tmean)<<" / "<<(1-init_etotal/final_etotal)<<" heat ="<<species[0].dS(j)<<" cool "<<species[e_idx].dG(j)<<" "<<(species[e_idx].dG(j) + species[e_idx].dGdT(j)*100.)<<endl;
                    cout<<endl;
                }
            }

            int dd = 0;
            if((steps>=10000) && (j==122) && (steps <= 10010)){ //Choose a time and cell to look into
                dd=0;
            }
            int found_solution=0;
            int num_cycles = num_subcycles;
            const int max_heat_cycles = 5;
            while(!subcycle_heat_exchange(j, num_cycles++, debug=dd, dt)) {
                
                if(num_cycles>max_heat_cycles) break;
            }
            

            double avgT_nom = 0;
            double avgT_denom = 0;

            //
            // Make checks and corrections
            //
            for(int si=0; si<num_species; si++) {
                
                double tt = coll_heat_output(si); //old solver
                //double tt = std::min(coll_heat_output(si), 1e4); //old solver
                //double tt = eta1[idx_s]/denoms[idx_s] + eta2[idx_s]/denoms[idx_s]*Jrad_FLD(j, 0);
                //double tt = species[si].prim[j].temperature; // this temperature now comes freshly out of the subcycle

                int idx_s = j * (num_species) + si;
                
                if(tt<temperature_floor)
                        tt=temperature_floor;
                    
                if(tt> ( max_temperature + globalTime/max_temperature_time * 1e4) )
                        tt=(max_temperature + globalTime/max_temperature_time * 1e4);         

                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
                avgT_denom += species[si].u[j].u1 * species[si].cv;

                if(globalTime > 1e-30)
                    species[si].prim[j].temperature = tt ;
                else
                    species[si].prim[j].temperature = species[si].const_T_space;
            }

            if(use_avg_temperature && globalTime < avg_temperature_t1){
                
                if(globalTime < avg_temperature_t0) {

                    double avgtemp    = avgT_nom/avgT_denom;
                    double relaxtemp;

                    double rt = shadow_relaxation_time; //shadowed regions temperature relaxation timescale
                    double totalheat = 0;
                    for(int si=0; si<num_species; si++) 
                        totalheat += species[si].dS(j);

                    
                    if(use_shadow_relaxation && (totalheat < shadow_relaxation_threshold) && (x_iVC[j]< shadow_relaxation_radius) ) { //1e-50 is the arbitrary limit we set on the heating function throughout the code
                        double fac = ddt/rt;
                        double f1  = (1 + 1e-40) / (1 + fac + 1e-40);
                        double f2  = (1 + 1e-40) / (1/fac + 1 + 1e-40);

                        relaxtemp = avgtemp * f1 + species[0].const_T_space * f2;
                    }
                        else
                            relaxtemp = avgtemp;

                    for(int si=0; si<num_species; si++) {
                            species[si].prim[j].temperature = relaxtemp;
                    } 
                            
                } else {
                        for(int si=0; si<num_species; si++) {
                                species[si].prim[j].temperature = avgT_nom/avgT_denom * (1. - globalTime/avg_temperature_t1)  + species[si].prim[j].temperature * globalTime/avg_temperature_t1;
                        }
                }
            }//No else case need: temperatures already assigned

        } //end j loop

    } else { //dont use coll heating or num_species !> 1
        
            
        for (int j=0; j < num_cells + numcells_offset; j++){
            for(int s=0; s<num_species; s++) {
                
                int idx_s = j * (num_species) + s;
                double tt = eta1[idx_s]/denoms[idx_s]; 
                if(couple_J_into_T)
                    tt += eta2[idx_s]/denoms[idx_s]*Jrad_FLD(j, 0);
                
                if( j==2 && steps>3364 && false)
                    cout<<"t = "<<steps<<" j ="<<j<<" T = "<<tt<<" num_cells = "<<num_cells<<endl;
                
                
                if(tt < 0. && (j<num_cells+numcells_offset)) {
                //if(steps == 221160) {
                    cout<<" negative T in s = "<<species[s].speciesname<<" j/s = "<<j<<"/"<<s<<" eta1/eta2/J = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<"/"<<Jrad_FLD(j, 0)<<" denom/eta2*J = "<<denoms[idx_s]<<"/"<<eta2[idx_s]*Jrad_FLD(j,0)<<" t/dt/steps = "<<globalTime<<"/"<<ddt<<"/"<<steps<<endl;
                    Tswitch = 1;
		            cout<<"tempers[s] = ";
		            for(int ss=0; ss<num_species; ss++) { 
                           int idx_ss = j*num_species + ss;
                           cout<<eta1[idx_ss]/denoms[idx_ss]<<" ";
                     }
                }
                
                if(tt<temperature_floor)
                    tt=temperature_floor;
                
                if(tt>max_temperature)
                    tt=max_temperature;
                
                if(Jswitch == 0)
                    species[s].prim[j].temperature = tt ;
                

		//if( (steps==0 || steps==1) || j==300)
         //            cout<<"t = "<<steps<<" j ="<<j<<" FINAL T = "<<tt<<" num_cells = "<<num_cells<<endl;
            }
        }
        
        if( steps>10e99) {
            char a;
            cin>>a;
        }
       
    }
    
    //
    // Conduction
    //
    if(use_conduction && ( globalTime < do_cond_until)) {
        std::vector<double> temp_temperatures        = std::vector<double>(num_cells+1);
        //std::vector<double> flux_temperatures        = std::vector<double>(num_cells+2);

	    temp_temperatures[1] = species[0].const_T_space;        
        for (int j=2; j < num_cells+numcells_offset-1; j++){
                temp_temperatures[j] = species[0].prim[j].temperature;
                
                double n_tot = 0;
		        double n_neutrals =0;
                double mumean = 1.;
                double mumean_nom = 0.;
                for(int s=0; s<num_species; s++) {
                    n_tot      += species[s].prim[j].number_density;
                    mumean_nom += species[s].prim[j].density;
                    if(species[s].static_charge == 0)
                        n_neutrals += species[s].prim[j].number_density;
                }
                mumean = mumean_nom / n_tot;
                
                double dT1     = (species[0].prim[j-1].temperature - species[0].prim[j].temperature) / (x_i12[j-1]-x_i12[j]);
                double  T1avg  = (species[0].prim[j-1].temperature + species[0].prim[j].temperature) * 0.5;
                double dT2     = (species[0].prim[j].temperature -   species[0].prim[j+1].temperature) / (x_i12[j]-x_i12[j+1]);
                double  T2avg  = (species[0].prim[j].temperature +   species[0].prim[j+1].temperature) * 0.5;
                
                double kappa_cond = n_neutrals/n_tot * conductivity + (1. - n_neutrals/n_tot) * conductivity2; //1e-5
                //double kappa_cond = conductivity; //1e-5
                double vfactor1    = std::max((1. - std::fabs(species[0].prim[j-1].speed/species[0].prim[j-1].sound_speed) ), 0.);
                double vfactor2    = std::max((1. - std::fabs(species[0].prim[j].speed/species[0].prim[j].sound_speed) ), 0.);
                double c1 = kappa_cond * std::pow(T1avg, 0.7); //*n_tot
                double c2 = kappa_cond * std::pow(T2avg, 0.7); // *n_tot
                
                temp_temperatures[j] +=  - dt * (c1 * dT1 * vfactor1 * surf[j-1] - c2 * dT2 * vfactor2 * surf[j]) / vol[j] ;

		if(j==-20) {
			cout<<"inside conduction vfactor1 = "<<vfactor1<<" vf2 "<<vfactor2<<" kappa_cond "<<kappa_cond<<endl;
		}

                
        }
	//temp_temperatures[2] = species[0].const_T_space;
        
        for (int j=1; j < num_cells+numcells_offset-1; j++){
                
            if(j==-20 && globalTime > 1.) {
                cout<<" temperatures before and after: "<<species[0].prim[j].temperature<<" "<<temp_temperatures[j]<<" delta = "<<(1.-species[0].prim[j].temperature/temp_temperatures[j])<<endl;
		        char aa;
		        cin>>aa;
             }
                
                for(int s=0; s<num_species; s++) {
                      species[s].prim[j].temperature  =  temp_temperatures[j];
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
        bool electrons = 0;
        
        for (int j=0; j < num_cells+numcells_offset; j++){
            for(int s=0; s<num_species; s++) {
                
                double dx = (x_i12[j+1]-x_i12[j]) ;
                double rhoavg = (species[s].prim[j].density + species[s].prim[j+1].density) * 0.5;
                double Pavg   = (species[s].prim[j].pres + species[s].prim[j+1].pres) * 0.5;
                double Tavg   = (species[s].prim[j].temperature + species[s].prim[j+1].temperature) * 0.5;
                double dP     = (species[s].prim[j].pres - species[s].prim[j+1].pres)/Pavg;
                double dT     = (species[s].prim[j].temperature - species[s].prim[j+1].temperature)/Tavg / dP;
                //double dTabs  = (species[s].prim[j].temperature - species[s].prim[j+1].temperature);
                double glocal = -get_phi_grav(x_i[j], enclosed_mass[j])/x_i[j];
                            
                double nabla_ad = 1.-1./species[s].gamma_adiabat;
                double lam = Pavg / (species[s].prim[j].pres - species[s].prim[j+1].pres);  // The mixing length
                double DT =  (dT > nabla_ad ? dT - nabla_ad : 0.); //smooth(dT, nabla_ad); //   // Gradient comparison and switch for Lconv
                       DT = (dx * total_opacity(j,0)) > 2./3. ? DT : 0.; //Guardian to not use convection in optically thin areas
                            
                double alphaconv = 0.5 * species[s].cv * lam * lam * DT * rhoavg * std::sqrt(glocal/Tavg); //Prefactor

                if(electrons) {
                    
                    double denom = 0.;
                    
                    for(int sj=0; sj<num_species; sj++) {
                        if(s!=sj){ //Also if sj == neutral
                            
                            double Q = 1.;
                            
                            denom += species[sj].prim[j].number_density * Q;
                        }
                    }
                    
                    denom *= 3.22e4 * Tavg * Tavg * species[s].prim[j].number_density;  //Eqn. 5.146 in Schunk
                    denom += 1.;
                
                    lam = 7.7e5 * Tavg * Tavg * std::pow(Tavg, 0.5) / denom;
                    
                }
                else {
                    lam = 0.*alphaconv;//Placeholder code for the compiler to stfu//For general expression need 5.167 with 4.130a and the collision integrals 4.114
                }
                //TODO: Couple with solution matrix once we have a good idea how to do the convection in the simple radiation solver
            }
            
        }
        
    }
        

    
    
    // Update energies. 
    // TODO: We should add cv * (Tf - Ti) to u to conserve energy properly.
    for(int si=0; si<num_species; si++) {
        for(int j =0; j<= num_cells+1; j++) {
            /* 
            double tt = tmp.u3 - 0.5 * tmp.u2*tmp.u2/tmp.u1;
            tt /= (tmp.u1*species[si].cv);
            */
            AOS_prim tmpp = species[si].prim[j];
            AOS      tmp  = species[si].u[j];

            if(std::isnan(species[si].prim[j].temperature))
                cout<<" @end of T_solver_simple "<<steps<<" "<<j<<" "<<si<<" "<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tmpp.temperature<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;
            species[si].prim[j].temperature = std::min(std::max(species[si].prim[j].temperature, temperature_floor), max_temperature );
        }

        species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);

        species[si].eos->compute_auxillary(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        

        //species[si].eos->compute_primitive(&(species[si].u[0]), &(species[si].prim[0]), num_cells+2) ;    
        
    }


    for(int si=0; si<num_species; si++) {
        for(int j =0; j<= num_cells+1; j++) {
            AOS      tmp  = species[si].u[j];
            AOS_prim tmpp = species[si].prim[j];
            double tt = tmp.u3 - 0.5 * tmp.u2*tmp.u2/tmp.u1;
            tt /= (tmp.u1*species[si].cv);
            if(std::isnan(tmpp.temperature) || std::isnan(tt) )
                cout<<" after T_solver "<<steps<<" "<<j<<" "<<si<<" "<<tt<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tmpp.temperature<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;                
        }
    }

}



/*
Function subcycle heat exchange
    Splits the heat exchange and heating operator for all species into smaller substeps
*/
int c_Sim::subcycle_heat_exchange(int j, int num_cycles, int debug, double dt) {

            //if(num_cycles>num_subcycles)
            //    cout<<" in subycles, being called at j ="<<j<<" with num_cycles = "<<num_cycles<<endl;

            int scale_back = 1;
            double scl_fac = 1 ;
            double dgdt_mul = 0;

            for(int si=0; si<num_species; si++) {
                AOS      tmp  = species[si].u[j];
                AOS_prim tmpp = species[si].prim[j];
                double   tt     = tmpp.temperature; //tmp.u3 - 0.5 * tmp.u2*tmp.u2/tmp.u1;
                
                if(std::isnan(tt))
                    cout<<" @start T subcycling found NaN! "<<steps<<" "<<j<<" "<<si<<" "<<tt<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tt<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;                
            }
            
            coll_heat_matrix.setZero();
            coll_heat_matrix_fixed.setZero();
            coll_heat_b.setZero();
            coll_heat_b_fixed.setZero();
            coll_heat_output.setZero();
            tmp_temperatures.setZero();
            double e_init = return_e_total(j);

            Matrix_t documentation         = Matrix_t::Zero(num_species+1, num_cycles+1); //last species row is the mean temperature

            if( (update_coll_frequently && j>=2) || (friction_solver==0) ) {
                fill_alpha_basis_arrays(j); 
                compute_alpha_matrix(j);
            }
            
            //********************************************************************** */
            //********************************************************************** */
            //Build fixed matrix which we will reuse for each cycle
            //********************************************************************** */
            //********************************************************************** */
            for(int si=0; si<num_species; si++) {

                if(scale_back)
                    scl_fac = dt/species[si].cv;
                
                int idx_s                     = j * (num_species) + si;
                coll_heat_b_fixed(si)         = 0.;
                if(couple_J_into_T)
                    coll_heat_b_fixed(si)         += rhs_sc(idx_s) ;//eta1[idx_s]; 4.* pi * kappa * Jrad_FLD(j, 0);

                double diag_sum = 0;
                double temp = 0;
                for(int sj=0; sj<num_species; sj++) {
                    temp      = scl_fac * 1. * friction_coefficients(si,sj) * 3 * kb / (mass_vector(si) + mass_vector(sj)) ;
                    diag_sum += temp;
                    coll_heat_matrix_fixed(si,sj) -= temp;
                }

                coll_heat_matrix_fixed(si,si) += diag_sum; 
                coll_heat_matrix_fixed(si,si) += scl_fac * dgdt_mul * species[si].dGdT(j) / species[si].u[j].u1; //cooling gradient term
                
                //Initialize temperatures and document evolution
                tmp_temperatures(si) = species[si].prim[j].temperature;
                documentation(si,0)  = species[si].prim[j].temperature;
                documentation(num_species,0)  = return_T_mean(j); //Document mean temperature before first step
            }

            //cout<<endl<<friction_coefficients<<endl;

            //********************************************************************** */
            //********************************************************************** */
            // Begin subcycling
            //********************************************************************** */
            //********************************************************************** */

            double ddt = dt/((double)num_cycles);

            for(int c=0; c<num_cycles; c++) {
                
                //Reset matrices for this cycle
                coll_heat_matrix  = coll_heat_matrix_fixed;
                coll_heat_b       = coll_heat_b_fixed;
                double cooling_temp = documentation(num_species,c); //last mean temperature //std::max(documentation(e_idx,c), 3.)
                //
                // Temperature averaging in last iteration to obtain higher order estimate of the r.h.s of the heat equation
                //
                if(c==num_cycles-1) {
                    tmp_temperatures.setZero();

                    if(debug >= 1)
                        cout<<" in heat subscycling, order = "<<num_cycles<<" Adding the following steps to average: "<<steps<<" time "<<this->globalTime;
                    
                    
                    if(num_cycles==1) { //If more than 
                            for(int si=0; si<num_species; si++) {
                                tmp_temperatures(si) += documentation(si,0);
                            }
                    } else { //If more than 1, exclude 1 as this can have extreme t-differences with solutions
                        for(int cc=1; cc<=num_cycles-1; cc++) {
                            if(debug >= 1)
                                cout<<cc<<" ";
                            for(int si=0; si<num_species; si++) {
                                tmp_temperatures(si) += documentation(si,cc);
                            }
                        }
                        tmp_temperatures /= (double)(num_cycles-1);
                    }

                    if(debug >= 1)
                        cout<<" avg over "<<num_cycles-1<<" cycles. ";

                     //average the temperatures over the course of the subcycle
                    if(e_idx>-1)
                        cooling_temp = tmp_temperatures(e_idx); 

                    if(debug>=1)
                        cout<<" Resulting electron cooling temperature = "<<cooling_temp<<endl;
                    
                    ddt = dt;
                }
                if(debug >= 1)
                    cout<<" cycle "<<c<<" dt fraction "<<ddt/dt<<" cooling temp "<<cooling_temp<<endl;
                
                //Recompute cooling with new temp temperatures
                
                if(e_idx>-1) {
                    if(cooling_temp < 0) {
                        if(debug>=1)
                            cout<<" Exiting subcycling! ";
                        return 0;
                    }
                        
                    do_highenergy_cooling(j, cooling_temp); 
                }
                //Update with mean temperature, as intermediate electron temperatures can be extremely high
                //do_highenergy_cooling(j, std::max(documentation(e_idx,c), 3.) ); //TODO: Add also update for kappa_planck, to allow bolometric cooling to converge
                
                //
                // Add remaining terms for this cycle 
                //
                for(int si=0; si<num_species; si++) {

                    if(scale_back)
                        scl_fac = ddt/species[si].cv;

                    int idx_s = j * (num_species) + si;

                    double Tsold = documentation(si,0); //(si,c) //This should make the method second order: Tsold is on the r.h.s of the equation, while the updated fractional-step temperatures are used to estimate new gradients
                    double Ts    = documentation(si,c);
                    double Ts3   = Ts*Ts*Ts;
                    double kappa = 1.*species[si].opacity_planck(j, 0);
                    //
                    // Rad equilibrium terms

                    coll_heat_matrix(si,si) += scl_fac * (species[si].cv/ddt     + 16 * sigma_rad * kappa * Ts3);
                    coll_heat_b(si)         += scl_fac * (species[si].cv/ddt     + 12 * sigma_rad * kappa * Ts3) * Tsold ;
                    coll_heat_b(si)         += scl_fac * (species[si].dS(j)      - species[si].dG(j)  + dgdt_mul * species[si].dGdT(j) * Tsold )/species[si].u[j].u1;

                    
                    /////////////////////////////////////////////////////////////////////
                    /////////////////////////////////////////////////////////////////////
                    // Coll heat matrix reintroduced, but should be put back into "fixed" part of the code
                    // Variant of the code which recomputes alpha friction coefficients as well during subcycles -> very expensive
                    /*
                    double diag_sum = 0;
                                     double temp = 0;
                    
                    for(int sj=0; sj<num_species; sj++) {
                        temp      = friction_coefficients(si,sj) * 3 * kb / (mass_vector(si) + mass_vector(sj)) ;
                        diag_sum += temp;
                        coll_heat_matrix(si,sj) -= temp;
                    }

                    coll_heat_matrix(si,si) += diag_sum;
                    coll_heat_matrix(si,si) += dGdT_mul * species[si].dGdT(j) / species[si].u[j].u1;
                    */
                    /////////////////////////////////////////////////////////////////////
                    /////////////////////////////////////////////////////////////////////
                }

                //Still code variant which recombutes collision alphas, might need to reimplement if instabilities occur
                /* 
                LU.compute(coll_heat_matrix) ;
                coll_heat_output.noalias() = LU.solve(coll_heat_b);

                tmp_temperatures = coll_heat_output;
 
                Vector_t tmp_results = return_preconditioned_LU_solution(coll_heat_matrix, coll_heat_b, LU, j);
                */
                if(debug >=1) {
                    cout<<" b = ";
                    for(int ss=0; ss<num_species; ss++){ cout<<" "<<coll_heat_b(ss); }
                    cout<<endl; 
                }

                if(steps < 10e99) {
                    LU.compute(coll_heat_matrix) ;
                    coll_heat_output.noalias() = LU.solve(coll_heat_b);

                    //tmp_temperatures = coll_heat_output;
                } else {
                    coll_heat_output = return_preconditioned_LU_solution(coll_heat_matrix, coll_heat_b, temperature_vector, LU, j); //note: temperature_vector contains the initial temps going into this routine and are filled in fill_alpha_basis_matrix_thingies
                }

                //
                for(int si=0; si<num_species; si++) {
                    documentation(si,c+1) = coll_heat_output(si);
                }
                documentation(num_species,c+1)  = return_T_mean(j, coll_heat_output);
                //cycle complete
            }


        //********************************************************************** */
        //********************************************************************** */
        // All cycles complete   
        //********************************************************************** */
        //********************************************************************** */
       
        //********************************************************************** */
        //********************************************************************** */
        // Show time evolution of temperatures    
        //********************************************************************** */
        //********************************************************************** */
        
        if( (debug >= 1)) {
            cout<<endl<<"step "<<steps<<" j "<<j<<" dt = "<<dt <<" Showing time evolution over subcycles "<<endl;
            for(int si=0; si<num_species; si++) {
                
                cout<<si<<"   ";
                for(int c=-1;c<num_cycles; c++)
                    cout<<documentation(si,c+1)<<" ";
                cout<<" dS/dG/dGdT = "<<species[si].dS(j)<<" "<<species[si].dG(j)<<" "<<species[si].dGdT(j)<<endl;
                
            }
            cout<<"avg   "<<documentation(num_species,0)<<" ";
                for(int c=0;c<num_cycles; c++)
                    cout<<1-documentation(num_species,0)/documentation(num_species,c+1)<<" ";
            
            cout<<" kappa ";
            for(int si=0; si<num_species; si++)
                cout<<species[si].opacity_planck(j, 0)<<" ";

            // Heating
            double de_total_expected;
            for(int si=0; si<num_species; si++) {
                de_total_expected += dt  * (species[si].dS(j) - species[si].dG(j) - species[si].dGdT(j) * (documentation(si,num_cycles)-documentation(si,0)  )   );
            }
            cout<<" rel. de expected from heating: "<<de_total_expected/e_init<<endl ;
        }
    
 
        //********************************************************************** */
        //********************************************************************** */
        // After checks are complete, write last solution into original temperature   
        //********************************************************************** */
        //********************************************************************** */
       
    int allgood = 1;
    for(int si=0; si<num_species; si++) {
        AOS      tmp  = species[si].u[j];
        AOS_prim tmpp = species[si].prim[j];
        double tt = documentation(si, num_cycles-1); 
        
        if(std::isnan(tt) || (tt<0) || ( (j==300) && (steps==10000))) {
            cout<<" @end T subcycling found NaN or <0! "<<steps<<" "<<j<<" "<<si<<" "<<tt<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tt<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;   
            cout<<" Temperature history: in subcycling: ";
            for(int c=0; c<=num_cycles; c++) {
                cout<<documentation(si, c)<<" ";
            }
            cout<<endl;
            
            allgood =0;
        }
    }

    if(!allgood)
        return 0;
    
    if(debug <= 1) {

        for(int si=0; si<num_species; si++) {
            species[si].prim[j].temperature = documentation(si, num_cycles-1); //ignore last computation
        }
        
        if( (debug >= 1) && (e_idx > -1))
            cout<<" Found temperature solution with num_cycles="<<num_cycles<<", returning. FInal electron temperature ="<<species[e_idx].prim[j].temperature<<endl;
    }
    return 1; //allgood

}




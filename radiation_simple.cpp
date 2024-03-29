
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
    
    if(debug > 1)
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
    
    int numcells_offset = 1; //Nominally 1
    // Step 1: setup transport terms (J)
    for(int b=0; b<num_bands_out; b++) {
        for (int j=0; j < num_cells + numcells_offset; j++) {
            int idx = j*stride + b*(num_vars + 1) ;
            int idx_r = j*num_vars + b ;

            // Time dependent terms:
            d[idx] +=  vol[j] / (c_light * ddt) ;
            r[idx_r] += (vol[j] / (c_light * ddt)) * Jrad_FLD(j, b) ;

            // Flux across right boundary
            if (j > 0 && j < num_cells + numcells_offset) {
                double dx      = (x_i12[j+1]-x_i12[j]) ;                
                double rhokr   = 0;
                
                rhokr = max(2.*(total_opacity(j,b)*total_opacity(j+1,b))/(total_opacity(j,b) + total_opacity(j+1,b)), 4./3./dx );
                rhokr   = min( 0.5*( total_opacity(j,b) + total_opacity(j+1,b)) , rhokr);
                
                double tau_inv = 1. / (dx * rhokr) ;
                double dJ = std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b))/(Jrad_FLD(j, b) + 1e-300);
                    
                double R       = 1.* xi_rad * tau_inv *  dJ ; // Put in 1.0 as prefactor to get correct rad shock
                double D       = 1.* tau_inv * surf[j] * flux_limiter(R);
                
                arr_rhokr[j] = rhokr;
                arr_R[j]        = R;
                arr_D[j]        = D;
                arr_JDJ[j]      =  std::abs(Jrad_FLD(j+1,b) - Jrad_FLD(j,b));///(4.*pi*(Jrad_FLD(j+1,b) + Jrad_FLD(j, b) + 1e-300));
                
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
                
                if(j==48000){
                    cout<<"j=48, fac1, fac2 = "<<Ts * ( 1. + fac * 12. * sigma_rad * Ts3)<<"/"<<ddt * (species[s].dS(j) + species[s].dG(j)) / species[s].u[j].u1 / species[s].cv;
                    cout<<" Ts, fac = "<<Ts<<"/"<<fac<<endl;
                }
                    
                double tempeta = 0;
                
                denoms[idx_s] = denom;
                eta1[idx_s] += Ts * ( 1. + 12. * fac);
                eta1[idx_s] += 1. * ddt * (species[s].dS(j) + species[s].dG(j) - photocooling_expansion * species[s].dGdT(j)*Ts     ) / species[s].u[j].u1 / species[s].cv;

                eta2[idx_s] += 4.* pi * ddt * kappa * no_rad_trans / species[s].cv;
                
                //if(false) {
                if(steps == 3363 && j==2) {
                    //cout<<"reporting cooling terms["<<s<<"]: dG / dGdT * Ts "<<species[s].dG(j)<<" / "<< - photocooling_multiplier * species[s].dGdT(j)*Ts <<endl;
		    cout<<steps<<" etas = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<" Ts = "<<Ts<<endl;
		}
                if(false) {
                //if(steps == 3115 && s==2) {
                    cout<<"reporting cooling terms["<<s<<"]: dG / dGdT * Ts "<<species[s].dG(j)<<" / "<< - photocooling_expansion * species[s].dGdT(j)*Ts <<endl;
                    //cout<<"reporting denom cooling terms : 16*fac / "<<16.*fac<<" / "<<moredenom<<endl;
                    //char a;
                   // cin>>a;
                    
                }
                
                //if(j == 90 && steps > 41800)
                if(false)
                    cout<<"j ="<<j<<" denoms = "<<denom<<" eta1 = "<<eta1[idx_s]<<" eta1 parts = Ts/dS "<<Ts * ( 1. + 12. * fac)<<"/"<<1. * ddt * (species[s].dS(j) + species[s].dG(j)) / species[s].u[j].u1 / species[s].cv<<" Ts parts = Ts/(1+12fac)/facp1/facp2 = "<<Ts<<"/"<<(1.+12.*fac)<<"/"<<ddt * kappa / species[s].cv * sigma_rad<<"/"<<Ts3<<endl; 
                
                
                if(j==48000) {
                    cout<<" eta1/eta2/tempeta1 = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<"/"<<tempeta<<endl;
                }
                
                exchange_d_sum += no_rad_trans * rhos * kappa * (1 - 4.*sigma_rad*Ts3/pi * eta2[idx_s]/denoms[idx_s]);
                exchange_r_sum += no_rad_trans * rhos * kappa * sigma_rad * Ts3/pi * (4 * eta1[idx_s]/denoms[idx_s] - 3 * Ts );
            }
            
           // cout<<" radiation part2, t = "<<steps<<" cell["<<j<<"] exchange_d/r = "<<exchange_d_sum<<"/"<<exchange_r_sum<<" etas = "<<eta1[j*(num_species)]<<"/"<<eta2[j*(num_species)];
            //cout<<" dS/rho/T = "<<species[0].dS(j)<<"/"<<species[0].u[j].u1<<"/"<<species[0].prim[j].temperature<<endl;
            
            int idx_b  = j*stride;
            int idx_rb = j*num_vars;
            
            d[idx_b]  += vol[j] * exchange_d_sum;
            r[idx_rb] += vol[j] * exchange_r_sum;
            
            exchange_d_sums[j] += exchange_d_sum;
            exchange_r_sums[j] += exchange_r_sum;
            //Sum up etas to get J source terms
            
        }
    }
    
    if(debug >= 3) {
        
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
    
    for (int j=0; j <= num_cells+numcells_offset; j++) {
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
    
    if(Tswitch == 1) {
        
        //char a;
        //cin>>a;
    }
    
    //
    // Compute Ti-Tj terms via a separate coupling step
    //
    
    if (use_collisional_heating && num_species > 1) {
        //compute_collisional_heat_exchange();
        
    
        for (int j=0; j < num_cells+numcells_offset; j++){
            
            Matrix_t coll_heat_matrix   = Matrix_t::Zero(num_species, num_species);
            Vector_t coll_heat_b        = Vector_t::Zero(num_species);
            Vector_t coll_heat_output   = Vector_t::Zero(num_species);
            double   tau = total_opacity(j,0) * (x_i12[j+1]-x_i12[j]);
            double security_multiplier = tau<1e3?1.:1.e-4;
            
            if(update_coll_frequently) {
                fill_alpha_basis_arrays(j);
                compute_alpha_matrix(j);
            }
            
            ////compute_collisional_heat_exchange_matrix(j);

            for(int si=0; si<num_species; si++) {
                
                int idx_s = j * (num_species) + si;
                
                //coll_heat_matrix(si,si) += 1./ddt ;
                //coll_heat_b(si)         += species[si].prim[j].temperature / ddt;
                
                //
                // Rad equilibrium terms
                //
                coll_heat_matrix(si,si) += denoms[idx_s];
                if(couple_J_into_T)
                    coll_heat_b(si)         += eta2[idx_s]*Jrad_FLD(j,0);
                coll_heat_b(si)         += eta1[idx_s];
                
                //for(int sj=0; sj<num_species; sj++) {
                //   coll_heat_matrix(si,sj) -= security_multiplier * friction_coefficients(si,sj);
                //}
                double diag_sum = 0;
                double temp = 0;
//                cout<<endl;
                for(int sj=0; sj<num_species; sj++) {
                    temp = ddt * friction_coefficients(si,sj) * 3 * kb / (species[si].cv * (mass_vector(si) + mass_vector(sj)) );
                    diag_sum += temp;
                    coll_heat_matrix(si,sj) -= temp;
                }

                coll_heat_matrix(si,si) += diag_sum;

            }
            
            LU.compute(coll_heat_matrix) ;
            coll_heat_output.noalias() = LU.solve(coll_heat_b);
            
            if(steps == 500e99 && j==175) {
                cout<<"j==5, steps=="<<steps<<", Coll Matrix / b was "<<endl<<coll_heat_matrix<<endl<<coll_heat_b<<endl;

		cout<<"invtotalmasses = "<<inv_totmasses<<endl;

                for(int si=0; si<num_species; si++) {
		    int idx_s = j * (num_species) + si;
                    cout<<" After coll, species = "<<si<<" Tbefore/Tafter = "<<species[si].prim[j].temperature<<"/"<<coll_heat_output(si)<<" rel. difference = "<<1-coll_heat_output(si)/species[si].prim[j].temperature<<endl;
                    cout<<" +heating/+cooling/-dGdT = "<<species[si].dS(j)<<"/"<<species[si].dG(j)<<"/"<<photocooling_expansion * species[si].dGdT(j)*species[si].prim[j].temperature<<" sum = "<<(species[si].dS(j)+species[si].dG(j)-photocooling_multiplier * species[si].dGdT(j)*species[si].prim[j].temperature)/ species[si].u[j].u1 / species[si].cv<<endl;
		    cout<<"eta1/eta2/denom = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<"/"<<denoms[idx_s]<<" alpha_ij = ";
			for(int sj=0; sj<num_species; sj++) {
				cout<<" "<<friction_coefficients(si,sj);
			}
                    cout<<endl;
                    
                    //char a;
                    //cin>>a;
                }
                
            }
            double avgT_nom = 0;
	    double avgT_denom = 0;

            for(int si=0; si<num_species; si++) {
                
                double tt = coll_heat_output(si);
                

                int idx_s = j * (num_species) + si;
                //double tt = eta1[idx_s]/denoms[idx_s] + eta2[idx_s]/denoms[idx_s]*Jrad_FLD(j, 0);
                
                if(tt < 0. && j > 0) {
                    cout<<" negative T after TiTj s = "<<species[si].speciesname<<" j/s = "<<j<<"/"<<si<<" eta1/eta2/J = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<"/"<<Jrad_FLD(j, 0)<<" denom/eta2*J = "<<denoms[idx_s]<<"/"<<eta2[idx_s]*Jrad_FLD(j,0)<<" t/dt/steps = "<<globalTime<<"/"<<ddt<<"/"<<steps<<endl;
                    
                    for(int ss=0; ss<num_species; ss++)
                        cout<<" +heating/+cooling/-dGdT = "<<species[ss].dS(j)<<"/"<<species[ss].dG(j)<<"/"<<photocooling_expansion * species[ss].dGdT(j)*species[si].prim[j].temperature<<" sum = "<<(species[ss].dS(j)+species[ss].dG(j)-photocooling_multiplier * species[ss].dGdT(j)*species[si].prim[j].temperature)/ species[ss].u[j].u1 / species[ss].cv<<endl;
                        //cout<<" negative T after TiTj in s = "<<species[si].speciesname<<" j/s = "<<j<<"/"<<si<<" t/dt/steps = "<<globalTime<<"/"<<ddt<<"/"<<steps<<endl;
                        //Tswitch = 1;
                    
                    if(si==1)
                        cout<<"Matrix / b was "<<endl<<coll_heat_matrix<<endl<<coll_heat_b<<endl;
                }
                    
                if(tt<temperature_floor)
                        tt=temperature_floor;
                
                //if(j<105 && globalTime<1e4)
                //        tt=temperature_floor;
                    
                if(tt>max_temperature)
                        tt=max_temperature;
		if(j<=4)
			tt=species[si].const_T_space;                

		avgT_nom +=   species[si].u[j].u1 * species[si].cv * tt;
		avgT_denom += species[si].u[j].u1 * species[si].cv;

		
                if(globalTime > 1e-10)
                    species[si].prim[j].temperature = tt ;
                else
                    species[si].prim[j].temperature = species[si].const_T_space;
            }
        
            
	     if(use_avg_temperature && globalTime < avg_temperature_t1){
		if(globalTime < avg_temperature_t0) {
                   for(int si=0; si<num_species; si++) {
		       species[si].prim[j].temperature = avgT_nom/avgT_denom;
	               } 
                     } 
                  } else {
                   for(int si=0; si<num_species; si++) {
                         species[si].prim[j].temperature = avgT_nom/avgT_denom * (1. - globalTime/avg_temperature_t1)  + species[si].prim[j].temperature * globalTime/avg_temperature_t1;
                 }
			
              }

         } //end j loop
        
    } else {
        
            
        for (int j=0; j < num_cells + numcells_offset; j++){
            for(int s=0; s<num_species; s++) {
                
                int idx_s = j * (num_species) + s;
                double tt = eta1[idx_s]/denoms[idx_s]; 
		if(couple_J_into_T)
			tt += eta2[idx_s]/denoms[idx_s]*Jrad_FLD(j, 0);
                
                if( j==2 && steps>3364 && false)
                    cout<<"t = "<<steps<<" j ="<<j<<" T = "<<tt<<" num_cells = "<<num_cells<<endl;
                
                //if(j==140)
                //    cout<<"t = "<<steps<<" T = "<<tt<<endl;
                
                if(tt < 0.) {
                //if(steps == 221160) {
                    cout<<" negative T in s = "<<species[s].speciesname<<" j/s = "<<j<<"/"<<s<<" eta1/eta2/J = "<<eta1[idx_s]<<"/"<<eta2[idx_s]<<"/"<<Jrad_FLD(j, 0)<<" denom/eta2*J = "<<denoms[idx_s]<<"/"<<eta2[idx_s]*Jrad_FLD(j,0)<<" t/dt/steps = "<<globalTime<<"/"<<ddt<<"/"<<steps<<endl;
                    Tswitch = 1;
                }
                
                if(tt<temperature_floor)
                    tt=temperature_floor;
                
                if(Jswitch == 0)
                    species[s].prim[j].temperature = tt ;
                
            }
        }
        
        
        if( steps>10e99) {
            char a;
            cin>>a;
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
        species[si].eos->update_eint_from_T(&(species[si].prim[0]), num_cells+2);
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }

}

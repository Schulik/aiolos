/**
 *  advection.cpp
 *
 * This file contains the base routines to start and run a simulation, 
 * as well as hydrodynamics, i.e. the HLLC Riemann solver.
 */

#include "aiolos.h"

/**
 * Main simulation loop
 *  
 * After successful simulation init, execute the main simulation from globalTime = [0.,t_max]
 * Every loop iteration starts with an estimator for all species' mid-cell pressures, 
 * to be able to determine the sound-speed and hence the CFL factor.
 * Following are output checks and then the main code modules are executed in order.
 */
void c_Sim::execute() { 
    
    steps = 0;
    double output_counter = 0;
    double monitor_counter= 0;
    const double dt_initial = dt_min_init;
    double next_print_time  = dt_initial;
    const signed long long int maxsteps = 1e12;
    
    int crashed_T = 0, crashed_J = 0;
    int crash_T_imin = num_cells+2, crash_T_imax = 0, crash_T_numcells = 0;
    int crash_J_imin = num_cells+2, crash_J_imax = 0, crash_J_numcells = 0;
    double crashtime, crashed_temperature, crashed_meanintensity;
    int crashtime_already_assigned = 0;
    
    cout<<" VERSION 0.2"<<endl;
    cout<<endl<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfactor="<<cflfactor<<" and num_species = "<<num_species<<endl;
    if(num_species == 0) 
        throw std::invalid_argument("WARNING: No species specified! I cannot work like that. Aborting program.") ;
    
    //
    // Compute real, physical scales for orientation
    //
    
    scale_cs = std::sqrt(species[0].gamma_adiabat * (species[0].gamma_adiabat - 1.) * species[0].cv * species[0].prim[num_cells].temperature); // cm/s
    scale_rb = G*planet_mass / scale_cs / scale_cs; 
    scale_rh = planet_semimajor * au * pow(planet_mass / (3.* star_mass ),0.333333333333333333);
    
    scale_vk = std::sqrt(G*star_mass/(au*planet_semimajor));
    scale_time = scale_rb/scale_cs;
    
    cout<<"    Reporting base physical scales for selected problem in cgs units or other units that make sense;"<<endl;
    cout<<"    bondi_radius: "<<scale_rb<<" cm = "<<scale_rb/au<<" au\n" ;
    cout<<"    hill radius: " <<scale_rh << "cm = "<< scale_rh/au <<" au = "<<scale_rh/scale_rb<< " rb" << endl;
    cout<<"    velocities: vk = "<<scale_vk<<" cm/s, vk/cs = "<<scale_vk/scale_cs<<", cs = "<<scale_cs<< "cm/s" << endl;
    cout<<"    times, rb/cs/yr = "<<scale_time/year<<" rh/cs/yr"<<scale_rh/scale_cs/year<<endl;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {

        
        
    
          if(start_hydro_time > 0. && globalTime > start_hydro_time) {   //Comment in if a radiative equilibrium phase is desired before starting hydro
              do_hydrodynamics = 1;     
          }

        if(steps==0) {
            for(int s = 0; s < num_species; s++)
                species[s].compute_pressure(species[s].u);
            compute_total_pressure();
        }
            
        dt = get_cfl_timestep();
        dt = std::min(dt, timestep_rad2) ;
        dt = std::min(dt, t_max - globalTime) ;
        if(steps == 0)
            dt = std::min(dt, dt_initial);
        
        if( globalTime > next_print_time) {
            cout<<" Beginning step "<<steps<<" @ globalTime "<<globalTime<<" dt "<<dt;
            cout<< ", CFL " << cfl_step << ", energy dt " << timestep_rad2 << "\n";
            next_print_time *= 10.;
        }
         
        //
        // Save internal energy before we update it
        //
        for(int s = 0; s < num_species; s++)
            for(int i=num_cells; i>=0; i--)  {
                species[s].primlast[i].internal_energy = species[s].prim[i].internal_energy;
            }
        
        //
        // Step 0: Update gravity. Important for self-gravity and time-dependent smoothing length
        //
        if(debug >= 2)
            cout<<"Beginning timestep "<<steps<<endl;
        
        update_mass_and_pot();
        
	if(steps<10)
	        for(int s=0; s<num_species; s++)
        	    species[s].update_kzz_and_gravpot(s);  //Recomputes the homopause boundary and adjusts species-specific potentials
        
        if(debug >= 2)
            cout<<"Before fluxes... ";
        
        //
        // Output data, when required. Keep the output here, between the flux and conserved variable update, so that the 
        // zeroth output has the initialized conserved values, but already the first fluxes.
        //
        
        if(steps==0 || steps==99e99) {
            for(int s=0; s<num_species; s++) {
                species[s].print_AOS_component_tofile((int) output_counter);
            }
            print_monitor((int)monitor_counter);
            print_diagnostic_file((int)output_counter);
            
            monitor_counter+=1.;
            output_counter +=1.;
         }
         if(globalTime > output_counter*output_time + output_time_offset){
             if(debug >= 1)
                 cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
             
             print_diagnostic_file((int)output_counter);
             for(int s=0; s<num_species; s++)
                species[s].print_AOS_component_tofile((int)output_counter); 
             
             
             output_counter+=1.; 
         }
         if(globalTime > monitor_counter*monitor_time) {
            print_monitor(steps);
            
             monitor_counter += 1.;
        }
         
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
        //
        // Proper start
        // Do all explicit and implicit operations for one timestep on the entire grid for all species
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
        
        //
        // Step 0: Hydrodynamics, if so desired
        //

	//if (do_hydrodynamics == 1) 
        //    compute_drag_update() ;

        if(steps > debug_steps) {
            cout<<"t="<<steps<<" Pos 0 T[423]_s = ";
            for(int s = 0; s < num_species; s++) {
                cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
            }
            cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
        }
        
        if (do_hydrodynamics == 1) {
        
            //int cnt_broken_cells = 0;
            
            for(int k=0; k<=1; k++) {
                for(int s = 0; s < num_species; s++) {
                    species[s].execute(species[s].u, species[s].dudt[0], u_mask 0);
                }
                    
                
                //
                //March 19th 2024: Added get_cfl_timestep2() to sit here, to obtain the updated timestep based on the extrapolated left and right values - they can produce inconsistent fluxes with the cell-centered values, which the call of dt = get_cfl_timestep(); at the beginning of the timestep is based on;
                //
                dt = std::min(get_cfl_timestep2(), dt);
                
                for(int s=0; s < num_species; s++) {
                    species[s].u0 = species[s].u ;
                    for(int j=0; j < num_cells+2; j++)
                        species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
                }
                
                if( count_broken_cells(species[s].u, std::vector<double>&u_mask) == 0)
                    k=2;
            }
        }
        
        if(steps > debug_steps) {
            cout<<"t="<<steps<<" Pos 1 T[423]_s = ";
            for(int s = 1; s < 2; s++) {
                cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
            }
            cout<<" s1 numbers, p ="<<species[debug_species].prim[debug_cell].pres<<" eint = "<<species[debug_species].prim[debug_cell].internal_energy<<" rho = "<<species[debug_species].u[debug_cell].u1<<" mom = "<<species[debug_species].u[debug_cell].u2<<" E = "<<species[debug_species].u[debug_cell].u3;
            cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
            cout<<" fluxes in 423 and neighbours: "<<endl;
            for(int ll=-1; ll<=1; ll++)
                cout<<" "<<species[debug_species].dudt[0][debug_cell+ll].u1<<" "<<species[debug_species].dudt[0][debug_cell+ll].u2<<" "<<species[debug_species].dudt[0][debug_cell+ll].u3<<" E="<<species[debug_species].u[debug_cell+ll].u3<<" p="<<(species[debug_species].u[debug_cell+ll].u3 - 0.5*species[debug_species].u[debug_cell+ll].u2*species[debug_species].u[debug_cell+ll].u2/species[debug_species].u[debug_cell+ll].u1)<<" ekin = "<<0.5*species[debug_species].u[debug_cell+ll].u2*species[debug_species].u[debug_cell+ll].u2/species[debug_species].u[debug_cell+ll].u1<<endl;
        }
        
        if (order == IntegrationType::first_order) {
            globalTime += dt;
        } else if (order == IntegrationType::second_order) {
            // 2nd step evaluated at t+dt
            globalTime += dt;

            if(use_self_gravity==1)
                update_mass_and_pot();

            if (do_hydrodynamics == 1) {

                if (use_drag_predictor_step) {
                    for(int s = 0; s < num_species; s++)
                        species[s].compute_pressure(species[s].u);
                    compute_total_pressure();
                
                    compute_drag_update() ;
                    if (use_collisional_heating)
                        compute_collisional_heat_exchange() ; //Disable if radiation is used?
                }
                
                if(steps > debug_steps) {
                    cout<<"t="<<steps<<" Pos 1.05 T[423]_s = ";
                    for(int s = 0; s < num_species; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                    }
                    cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
                }

                //debug = 2;
                for(int s = 0; s < num_species; s++)
                    species[s].execute(species[s].u, species[s].dudt[1], 1);

                if(steps > debug_steps) {
                    cout<<"t="<<steps<<" Pos 1.1 T[423]_s = ";
                    for(int s = 0; s < num_species; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                    }
                    cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
                }
                
                for(int s = 0; s < num_species; s++){
                    for(int j = 0; j < num_cells+2; j++) {
                        // Need to re-compute u[j] if it was updated in the drag-predictor
                        if (use_drag_predictor_step)
                            species[s].u[j] = species[s].u0[j] + species[s].dudt[0][j]*dt ;
                        
                        species[s].u[j] += (species[s].dudt[1][j] - species[s].dudt[0][j])*dt / 2 ;  
                    }
                }
                
                
                
                if(steps > debug_steps) {
                    cout<<"t="<<steps<<" Pos 1.2 T[423]_s = ";
                    for(int s = 1; s < 2; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                    }
                    cout<<" s1 numbers, p ="<<species[debug_species].prim[debug_cell].pres<<" eint = "<<species[debug_species].prim[debug_cell].internal_energy<<" rho = "<<species[debug_species].u[debug_cell].u1<<" mom = "<<species[debug_species].u[debug_cell].u2<<" E = "<<species[debug_species].u[debug_cell].u3;
                    cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
                    cout<<" fluxes in 423 and neighbours: "<<endl;
                    for(int ll=-1; ll<=1; ll++)
                        cout<<" "<<species[debug_species].dudt[1][debug_cell+ll].u1<<" "<<species[debug_species].dudt[1][debug_cell+ll].u2<<" "<<species[debug_species].dudt[1][debug_cell+ll].u3<<" E="<<species[debug_species].u[debug_cell+ll].u3<<" p="<<(species[debug_species].u[debug_cell+ll].u3 - 0.5*species[debug_species].u[debug_cell+ll].u2*species[debug_species].u[debug_cell+ll].u2/species[debug_species].u[debug_cell+ll].u1)<<" ekin = "<<0.5*species[debug_species].u[debug_cell+ll].u2*species[debug_species].u[debug_cell+ll].u2/species[debug_species].u[debug_cell+ll].u1<<endl;
                }
                
                
            } else {
                for(int s = 0; s < num_species; s++) {
                    species[s].apply_boundary_left(species[s].u) ;
                    species[s].apply_boundary_right(species[s].u) ;
                }
            }
        }
        
        if(steps > debug_steps) {
                cout<<"t="<<steps<<" Pos 1.3 T[423]_s = ";
                for(int s = 0; s < num_species; s++) {
                    cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                }
                cout<<endl;
        }
        
        for(int s = 0; s < num_species; s++) 
            species[s].compute_pressure(species[s].u);

        if(steps > debug_steps) {
                cout<<"t="<<steps<<" Pos 1.5 T[423]_s = ";
                for(int s = 0; s < num_species; s++) {
                    cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                }
                cout<<endl;
            }
        
        //Computes the velocity drag update after the new hydrodynamic state is known for each species
        if (do_hydrodynamics == 1) 
            compute_drag_update() ;

        if (do_hydrodynamics == 0 && friction_solver > 0) 
                compute_drag_update() ;
        
        if(steps > debug_steps) {
                cout<<"t="<<steps<<" Pos 2 T[423]_s = ";
                for(int s = 0; s < num_species; s++) {
                    cout<<" ["<<s<<"]= "<<species[s].prim[debug_cell].temperature;
                }
                cout<<endl;
            }
        
        // If either switch is set we need to think more carefully about what should be done
        if( (photochemistry_level + use_rad_fluxes ) > 0 ) {
            
            update_opacities();
            if(photochemistry_level == 0 && use_rad_fluxes > 0)
                reset_dS();

            // Compute high-energy dS and ionization
            if(photochemistry_level == 1) {   //C2Ray scheme
                reset_dS();
                
                if(false) {
                    cout<<"Pos 1 dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                do_photochemistry();
                
                if(false) {
                    cout<<"Pos 2 dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                
            }
            else if(photochemistry_level == 2) { //Linearized time-dependent general chemistry scheme
                if(steps > 2) {
                    
                    if(steps % dt_skip_ichem == 0) {
                        dt_skip_dchem += dt;
                        reset_dS();
                        
                        //cout<<" Running chemistry with dt/dt_chem/dt_skip/steps = "<<dt<<" / "<<dt_skip_dchem<<" / "<<dt_skip_ichem<<" / "<<steps<<endl;
                        do_chemistry(dt_skip_dchem);
                        
                        dt_skip_dchem = 0.;
                    } else {
                        dt_skip_dchem += dt;
                    }
                }
                    
            }             
            update_dS();               //Compute low-energy dS
        
            
            if(false) {
                cout<<"Pos 3 dS_UV = "<<dS_band(num_cells-10,0)<<endl;
            }
            
            if(use_rad_fluxes==1) {
                update_fluxes_FLD();   //FLD Radiation transport, updating Temperatures and photon band energies
            }
            else if (use_rad_fluxes == 2){
                   
                update_fluxes_FLD_simple(dt); //'Simple' FLD solver
                
            }
        }
        
        if(steps > debug_steps) {
            for(int k=debug_cell;k<debug_cell+1;k++) {
                cout<<"t="<<steps<<"Pos 3 T["<<k<<"]_s = ";
                    for(int s = 0; s < num_species; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].prim[k].temperature;
                    }
                cout<<endl;
            }
            
        }
        
        if(steps==0)
            cout<<"Initial sound crossing time = "<<max_snd_crs_time<<", debug = "<<debug<<endl;
            
        steps++;
        
        if(debug > 1)
            cout<<"timestep in execute()="<<dt<<" stepnum "<<steps<<" totaltime"<<globalTime<<endl;
        
        //
        // Detection of negative J and T and soft exit
        //
        for(int s = 0; s < num_species; s++) {
            for(int i=num_cells-1; i>=0; i--)  {
                    if(species[s].prim[i].temperature < 0 || std::isnan(species[s].prim[i].temperature) ) {
                        
                        if(crashtime_already_assigned == 0) {
                            crashtime = globalTime;
                            crashtime_already_assigned = 1;
                        }
                        
                        crashed_temperature = species[s].prim[i].temperature;
                        crashed_T = s+1;
                        crash_T_imin = (i<crash_T_imin)? i : crash_T_imin;
                        crash_T_imax = (i>crash_T_imax)? i : crash_T_imax;
                        crash_T_numcells++;
                        
                        globalTime = 1.1*t_max; //This secures that the program exits smoothly after the crash and produces a -1 output file of the crashed state
                    }
            }
            
            if(crashed_T > 0) {
                cout<<endl<<">>> CRASH <<< DUE TO NEGATIVE TEMPERATURES, crash_imin/imax = "<<crash_T_imin<<"/"<<crash_T_imax<<" sample T = "<<crashed_temperature<<" num of crashed cells/total cells = "<<crash_T_numcells<<"/"<<num_cells<<"  crashed species name = "<<species[crashed_T-1].speciesname<<endl; 
                cout<<" @t/dt = "<<crashtime<<"/"<<dt<<" stepnum "<<steps<<endl;
                cout<<"Writing crash dump into last output and exiting program."<<endl;
            } 
                    
        }
        if(couple_J_into_T) {
        for(int b = 0; b < num_bands_out; b++) {
            for(int i=num_cells-1; i>=0; i--)  {
                    if(Jrad_FLD(i,b) < 0) {
                        
                        if(Jrad_FLD(i,b) > -1e-10) { //Ignore extremely small flux overshoots
                            Jrad_FLD(i,b) *= -1.; 
                        }
                        else {
                            
                            // If we're computing the planet's surface temperature
                            // allow J < 0 in the boundary.
                            if (i < num_cells && use_planetary_temperature &&
                                Jrad_FLD(num_cells,b) > 0)
                                continue ;
                            
                            if(crashtime_already_assigned == 0) {
                                crashtime = globalTime;
                                crashtime_already_assigned = 1;
                            }
                            
                            crashed_meanintensity = Jrad_FLD(i,b);
                            crashed_J = b+1;
                            crash_J_imin = (i<crash_J_imin)? i : crash_J_imin;
                            crash_J_imax = (i>crash_J_imax)? i : crash_J_imax;
                            crash_J_numcells++;
                            
                            globalTime = 1.1*t_max;    
                        }
                    }
            }
            
            if(crashed_J > 0 && couple_J_into_T) {
                cout<<endl<<">>> CRASH <<< DUE TO NEGATIVE J, crash_imin/imax = "<<crash_J_imin<<"/"<<crash_J_imax<<" sample J = "<<crashed_meanintensity<<" num of crashed cells/total cells = "<<crash_J_numcells<<"/"<<num_cells<<"  crashed band number = "<<crashed_J-1<<endl; 
                cout<<" @t/dt = "<<crashtime<<"/"<<dt<<" stepnum "<<steps<<endl;
                cout<<"Writing crash dump into last output and exiting program."<<endl;
            }
        }
        }
        
    }
    cout<<endl;
    
    //Print successful end result
    double checksum = 0.;
    for(int s = 0; s < num_species; s++) {
        for(int i=num_cells; i>=0; i--)  {
                checksum = (species[s].de_e[i] > checksum )?species[s].de_e[i]:checksum ;
        }
    }
    
    cout<<"Finished at time="<<globalTime<<" after steps="<<steps<<" with num_cells="<<num_cells<<" and checksum = "<<checksum<<endl;
    
    for(int s=0; s<num_species; s++) 
        species[s].print_AOS_component_tofile(-1);
    print_diagnostic_file(-1);
    print_monitor(-1);
        
}

/**
 * Compute the total left and right cell extrapolated cell-pressure of all species. Currently unused function.
 */
void c_Sim::compute_total_pressure() {
    
    for(int i=num_cells+1; i>=0; i--)  {
            
        total_press_l[i] = 0.;
        total_press_r[i] = 0.;
            
        for(int s = 0; s < num_species; s++) {
                //total_pressure[i] += species[s].prim[i].pressure;
                
            if(species[s].prim_l[i].pres > 0.) {
                total_press_l[i] += species[s].prim_l[i].pres;
                total_press_r[i] += species[s].prim_r[i].pres;
                
            }
                
            //if(i<5) cout<<" i/s = "<<i<<"/"<<s<<" pl/pr = "<<species[s].prim_l[i].pres<<"/"<<species[s].prim_r[i].pres<<" dens = "<<species[s].prim_l[i].density<<endl;
        }
            
        //cout<<"total pressure, cell "<<i<<" pl/pr = "<<total_press_l[i]<<"/"<<total_press_r[i]<<endl;
            
    }
}
    


////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  CLASS SPECIES
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * Compute first or second order time derivatives of conservative variables, i.e. flux gradients using intermediate primitives. Includes well-balanced gravity. 
 * Called from c_Sim::execute()
 * 
 * @param[in] u_in Vector of conservative variables over the entire grid.
 * @param[in] orderstep assures that the first order step is a strictly flux conserving step, by switching off the slope predictor for rho and v in the first step
 * @param[out] dudt Vector of time derivatives of conservative variables over entire grid.
 */
void c_Species::execute(std::vector<AOS>& u_in, std::vector<AOS>& dudt, std::vector<AOS>& u_mask, int orderstep) {
    
    
        if(base->steps > base->debug_steps && this_species_index == base->debug_species) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 1 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<endl;
        }
        
        //
        // Step 1: Boundary values
        //
        //
  
        apply_boundary_left(u_in) ;
        apply_boundary_right(u_in) ;
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 2 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
        }
        
        if(USE_WAVE==1) 
            add_wave(u_in, base->globalTime);
        
        if(debug >= 2)
            cout<<"Done."<<endl<<" Before compute pressure in species "<<speciesname<<"... ";
        
        if(base->debug > 0) { //Can put in custom debug info here. We left some example code in just for referencing the code usage.
            
            if(prim[21].temperature < 0. ) {
                cout<<" T<0 in execute, steps ="<<base->steps<< " beforeP, T = "<<prim[20].temperature<<" "<<prim[21].temperature<<" "<<prim[22].temperature<<" p_nominal = "<<(u[21].u3 - 0.5*std::pow(u[21].u2,2.)/u[21].u1)<<endl;
            }
        
            if( (base->steps == 3205600) || (base->steps == 3206000) ) { //Narrow down problematic time
                
                for(int j=19; j<=22; j++) { //Narrow down problematic cells, in this case rad fluxes do something unexpected
                    
                    double dx      = (base->x_i12[j+1]-base->x_i12[j]) ;
                    double rhokr   = max(2.*(base->total_opacity(j,0)*base->total_opacity(j+1,0))/(base->total_opacity(j,0) + base->total_opacity(j+1,0)), 4./3./dx );
                            rhokr   = min( 0.5*( base->total_opacity(j,0) + base->total_opacity(j+1,0)) , rhokr);
                    double tau_inv = 0.5 / (dx * rhokr) ;
                    double R       = 2 * tau_inv * std::abs(base->Jrad_FLD(j+1,0) - base->Jrad_FLD(j,0)) / (base->Jrad_FLD(j+1,0) + base->Jrad_FLD(j, 0) + 1e-300) ;
                    double flux_limiter;
                        if (R <= 2)
                            flux_limiter = 2 / (3 + std::sqrt(9 + 10*R*R)) ;
                        else 
                            flux_limiter = 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
                    
                    double D       = base->surf[j] * flux_limiter * tau_inv;
                    double flux    = - 4. * pi * D * (base->Jrad_FLD(j+1,0) - base->Jrad_FLD(j,0));
                
                    cout<<" i/i+1= "<<j<<"/"<<j+1<<" F = "<<flux<<" D = "<<D<<" limiter = "<<flux_limiter<<" R = "<<R<<" rhokr = "<<rhokr<<" J-B = "<<prim[j].density*opacity_planck(j,0)*(base->Jrad_FLD(j,0)-sigma_rad*pow(prim[j].temperature,4.) / pi )<<endl;
                }
                
                //Narrow down why fluxes are problematic. Is it the density gradients?
                cout<<"At steps ="<<base->steps<< " beforeP, dens = "<<prim[20].density<<" "<<prim[21].density<<" "<<prim[22].density<<" p_nominal = "<<(u[21].u3 - 0.5*std::pow(u[21].u2,2.)/u[21].u1)<<endl;
            }
        }
        
        //
        // Start!
        //
        compute_pressure(u_in);
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 3 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
        }
        
        if(debug >= 2)
            cout<<"Done. Starting edge-states."<<endl;
        
        reconstruct_edge_states(u_mask, orderstep) ;
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 4 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
        }
        //
        // Step 2: Compute fluxes and sources
        //
        if (not is_dust_like) {
            for(int j=0; j <= num_cells; j++)
                flux[j] = hllc_flux(j); 
        }
        else {
            for(int j=0; j <= num_cells; j++)
                flux[j] = dust_flux(j);     
        }
        
        if(debug >= 2)
            cout<<"Done. Starting sources."<<endl;
        
        for(int j=1; j<=num_cells; j++) {
            source[j]          = source_grav(u_in[j], j);
            source_pressure[j] = AOS(0, -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
                                          base->source_pressure_prefactor_right[j] * prim_r[j].pres)  ,0); 
        }
        
        //
        // Step 3: Add it all up to update the conserved variables
        //
        
        for(int j=1; j<=num_cells; j++) {
            dudt[j] = (flux[j-1] * base->surf[j-1] - flux[j] * base->surf[j]) / base->vol[j] + (source[j] + source_pressure[j]) ;
            
            if( debug > 3) { //Or put in your own conditions
                char alpha;
                cout<<"Debuggin fluxes in cell i= "<<j<<" for species "<<speciesname<<" at time "<<base->steps<<endl; 
                cout<<"     fl.u1 = "<<flux[j-1].u1<<": fr.u1 = "<<flux[j].u1<<endl;
                cout<<"     fl.u2 = "<<flux[j-1].u2<<": fr.u2 = "<<flux[j].u2<<endl;
                cout<<"     fl.u3 = "<<flux[j-1].u3<<": fr.u3 = "<<flux[j].u3<<endl;
                cout<<"     Cartesian fluxes: Fl-Fr+s = "<<((flux[j-1].u2 - flux[j].u2)/base->dx[j] + source[j].u2)<<endl;
                cout<<"     D_surface/volume="<<(0.5*(base->surf[j]-base->surf[j-1])/base->vol[j])<<" vs. 1/dx="<<(1./base->dx[j])<<" AL = "<<base->surf[j]<<" AR = "<<base->surf[j-1]<<" VOL = "<<base->vol[j]<<endl;
                cout<<endl;
                cout<<"     Reminder, the following three lines must be +- equal for well-balancing:"<<endl;
                cout<<"     s = "<<source[j].u1<<"/"<<source[j].u2<<"/"<<source[j].u3<<endl;
                cout<<"     Al*Fl - Ar*Fr + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + source_pressure[j].u2)<<endl;
                cout<<"     dP/dr = "<<((prim_l[j].pres - prim_r[j].pres)/base->dx[j])<<endl;
                cout<<" "<<endl;
                cout<<"     sP = "<<"/"<<source_pressure[j].u1<<"/"<<source_pressure[j].u2<<"/"<<source_pressure[j].u3<<endl;
                cout<<"     Al*Fl - Ar*Fr = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) /base->vol[j])<<endl;
                cout<<"     dP/dr + S = "<<((prim_l[j].pres - prim_r[j].pres)/base->dx[j] + source[j].u2)<<endl;
                cout<<endl;
                cout<<"     Al*Fl - Ar*Fr + s = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2))<<endl;
                cout<<"     u1 : Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u1 * base->surf[j-1] - flux[j].u1 * base->surf[j]) / base->vol[j] + (source[j].u1 +source_pressure[j].u1))<<endl;
                cout<<"     u2 : Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2 +source_pressure[j].u2))<<endl;
                cout<<"     u3 : Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u3 * base->surf[j-1] - flux[j].u3 * base->surf[j]) / base->vol[j] + (source[j].u3 +source_pressure[j].u3))<<endl;
                cout<<" "<<endl;
                cout<<"     u1 : dudt[1] = "<<(dudt[j].u1)<<endl;
                cout<<"     u2 : dudt[2] = "<<(dudt[j].u2)<<endl;
                cout<<"     u3 : dudt[3] = "<<(dudt[j].u3)<<endl;
                
                cin>>alpha;
            }
        }
    
} 

/**
 * The HLLC Riemann solver with pressure-based wave-speed estimate, according to Toro(2007) and references therein.
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @return flux at cell interface j
 */
AOS c_Species::hllc_flux(int j) 
{
    int jleft = j, jright = j+1;
    AOS flux;
    int option = 0;
    
    //Speed of gas
    double ul = prim_r[jleft].speed;  
    double ur = prim_l[jright].speed; 
    
    double pl = prim_r[jleft].pres;  
    double pr = prim_l[jright].pres;
    double pl_e = pl;
    double pr_e = pr;
    if(base->use_total_pressure) {
            pl = base->total_press_r[jleft];
            pr = base->total_press_l[jright];
    }
    
    double dl = prim_r[jleft].density;  
    double dr = prim_l[jright].density;

    double mom_l = dl*ul ;
    double mom_r = dr*ur ;

    double El = dl*prim_r[jleft].internal_energy + 0.5*mom_l*ul ;
    double Er = dr*prim_l[jright].internal_energy + 0.5*mom_r*ur ;

    if( (debug > 2) && (j>0 && j<13) ) {
        cout<<" IN HLLC, j = "<<j<<" dl/dr = "<<dl<<"/"<<dr<<" ul/ul = "<<ul<<"/"<<ur<<" pl/pr = "<<pl<<"/"<<pr<<" El/Er = "<<El<<"/"<<Er;
    }
    
    //Speed of shocks
    double SL = ul - prim_r[jleft].sound_speed ;
    double SR = ur + prim_l[jright].sound_speed ;
    
    //Intermediate values in the star region, equations 10.30 -10.39 in Toro
    double SS     = ( pr-pl+mom_l*(SL - ul)-mom_r*(SR-ur) )/(dl*(SL - ul)-dr*(SR-ur) );
    
    if ((SL <= 0) &&  (SS >= 0)) {
        AOS FL         = AOS (mom_l, mom_l * ul + pl, ul * (El + pl_e) );
        double comp3_L = El/dl + (SS-ul)*(SS + pl_e/(dl*(SL-ul)));
        AOS US_L       = AOS(1,SS, comp3_L) * dl * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - AOS(dl, mom_l, El))  * SL;    
           
        flux = FS_L;
        option = 1;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (mom_r, mom_r * ur + pr, ur * (Er + pr_e) );
        double comp3_R = Er/dr + (SS-ur)*(SS + pr_e/(dr*(SR-ur)));
        AOS US_R       = AOS(1,SS, comp3_R) * dr * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - AOS(dr, mom_r, Er)) * SR;
        
        flux = FS_R;
        option = 2;
    }
    else if (SL >= 0) {
        AOS FL = AOS (mom_l, mom_l * ul + pl, ul * (El + pl_e) );
        flux = FL;
        option = 3;
    }
    else if (SR <= 0) {
        AOS FR = AOS(mom_r, mom_r * ur + pr, ur * (Er + pr_e) );
        flux = FR;
        option= 4 ;
    }
    if((debug > 2) && (j>0 && j<13))
        cout<<" opt = "<<option<<" SL/SS/SR = "<<SL<<"/"<<SS<<"/"<<SR<<" f = "<<flux.u1<<"/"<<flux.u2<<"/"<<flux.u3<<" ((f.u2-Pl)/Pl)-1 = "<<((flux.u2-pl)/pl)<<" ((f.u2-Pr)/Pr)-1 = "<<((flux.u2-pr)/pr)<<endl;
    return flux;
}

/**
 * Riemann problem for (almost) pressure-less dust following Leveque (2004), Pelanti & Leveque (2006)
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @return flux at cell interface j
 */
AOS c_Species::dust_flux(int j) 
{
    AOS_prim Wl = prim_r[j] ;
    AOS_prim Wr = prim_l[j+1] ;

    // Case 1: Vacuum central state, no flux
    if (Wl.speed < 0 && Wr.speed > 0)
        return AOS(0,0,0) ;

    // Case 2: delta-shock with speed SS
    double R = std::sqrt(Wr.density/Wl.density) ;
    double SS = (Wl.speed + R*Wr.speed) / (1 + R) ;

    auto flux = [](const AOS_prim& W) {
        double mom = W.density * W.speed ;
        double en = W.density*(W.internal_energy + 0.5*W.speed*W.speed) ;

        return AOS(mom, W.speed*mom + W.pres, W.speed*(en + W.pres)) ;
    } ;

    if (SS > 0) {
        return flux(Wl) ;
    } else if (SS < 0) {
        return flux(Wr) ;
    } else {
        return (flux(Wl) + flux(Wr)) * 0.5 ;
    }
}

/**
 * The higher-order hydro solver can return unphysical values, i.e. negative pressures from the reconstruction p=E-0.5*rho*u^2, when E is dominated by kinetic energy
 * or even E itself can become negative for e.g. electrons.
 * We try to fix this by identifying offending cells, and re-running their entire hydro update after cells have been found.
 * The slope correction in the hydro reconstruction in cell j is always done with (1-u_mask) as multiplier.
 * The original hydro run is performed with all mask values being 0. This function counts breaking cells using the data u, and setting their u_mask values to 1 in the offending cell j
 * as well as the surrounding cells j-1 and j+1 (as we cannot know which flux edge is broken without going even deeper into analysis).
 * 
 * @param[in] u Input hydrodynamic data to be checked for negative E and negative p
 * @param[out] u_mask Output mask values. Any offending cell is set from 0 to 1
 * @return the sum of all u_mask values. The main loop will react if non-zero values are returned.
 */
int c_Species::count_broken_cells(std::vector<AOS>&u, std::vector<double>&u_mask) {
    
    int cnt =0;
    std::vector<double> p_temps = np_somevalue(num_cells+2, 0.);
    
    //Predict pressures from this hydro data
    for(int j=0; j<=num_cells+1; j++) {
            p_temps[j] = u[j].u3 - 0.5*u[j].u2*u[j].u2/u[j].u1;
        }
    for(int j=0; j<=num_cells; j++) {
            //Cheap -E detection and -p detection
            u_mask[j] = std::signbit( u[j].u3 ) + std::signbit( p_temps[j] );
            cnt      += (double)u_mask[j];
        }

        return cnt;
}


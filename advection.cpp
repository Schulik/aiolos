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
void c_Sim::execute(int restartnumber, double restarttime_cmdline) { 
    
    steps = 0;
    //if(restartnumber != 0)
    //    steps = 11;
    
    double output_counter = 0;
    double monitor_counter= 0;
    const double dt_initial = dt_min_init;
    double next_print_time  = dt_initial;
    double next_output_time = log_time_start==-20 ? output_time + output_time_offset : std::pow(log_time_factor, log_time_start);
    const signed long long int maxsteps = 1e12;
    
    int crashed_T = 0, crashed_J = 0;
    int crash_T_imin = num_cells+2, crash_T_imax = 0, crash_T_numcells = 0;
    int crash_J_imin = num_cells+2, crash_J_imax = 0, crash_J_numcells = 0;
    double crashtime, crashed_temperature, crashed_meanintensity;
    int crashtime_already_assigned = 0;
    int logtime = log_time_start==-20 ? 0 : 1;
    
    const int printstuff_steps = 0;
    
    cout<<" VERSION 0.2"<<endl;
    cout<<endl<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfactor="<<cflfactor<<" and num_species = "<<num_species<<endl;
    if(num_species == 0) 
        throw std::invalid_argument("WARNING: No species specified! I cannot work like that. Aborting program.") ;
    
    //
    // Compute real, physical scales for orientation
    //
    
    scale_cs = std::sqrt(species[0].gamma_adiabat * (species[0].gamma_adiabat - 1.) * species[0].cv * species[0].prim[num_cells].temperature); // cm/s
    scale_rb = G*planet_mass / scale_cs / scale_cs; 
    scale_rh = planet_semimajor * au * pow(planet_mass / (3.* (star_mass+1e-10) ),0.333333333333333333);
    
    scale_vk = std::sqrt(G*star_mass/(au*planet_semimajor));
    scale_time = scale_rb/scale_cs;
    
    cout<<"    Reporting base physical scales for selected problem in cgs units or other units that make sense;"<<endl;
    cout<<"    bondi_radius: "<<scale_rb<<" cm = "<<scale_rb/au<<" au\n" ;
    cout<<"    hill radius: " <<scale_rh << "cm = "<< scale_rh/au <<" au = "<<scale_rh/(scale_rb+1e-10)<< " rb" << endl;
    cout<<"    velocities: vk = "<<scale_vk<<" cm/s, vk/cs = "<<scale_vk/scale_cs<<", cs = "<<scale_cs<< "cm/s" << endl;
    cout<<"    times, rb/cs/yr = "<<scale_time/year<<" rh/cs/yr"<<scale_rh/scale_cs/year<<endl;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    
    

    for (globalTime = restarttime; (globalTime < t_max) && (steps < maxsteps); ) {
     
        //if(steps>10)
        //feenableexcept(FE_INVALID);
     
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
        dt = std::min(dt, t_max - globalTime);
        if(steps == 0)
            dt = std::min(dt, dt_initial);
        //if(steps > 400) {
	//    cout<<" steps  = "<<steps<<" shrinking dt = "<<dt<<endl;
	//    dt = dt * 0.5;
        //}

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
            if(output_counter >= restartnumber) {
                for(int s=0; s<num_species; s++) {
                    species[s].print_AOS_component_tofile((int) output_counter);
                }
            }
            //print_monitor((int)monitor_counter);
            //print_diagnostic_file((int)output_counter);
            
            monitor_counter+=1.;
            output_counter +=1.;
         }
         if(cont_output_steps > -1) {
            if(steps%cont_output_steps == 0) {
                for(int s=0; s<num_species; s++) {
                        species[s].print_AOS_component_tofile(-999);
                    }  
            }
         }
         //if(globalTime > output_counter*output_time + output_time_offset){
         if(globalTime > next_output_time) {
             if(debug >= 1)
                 cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
             
             //print_diagnostic_file((int)output_counter);
             if(output_counter >= restartnumber) {
                for(int s=0; s<num_species; s++)
                    species[s].print_AOS_component_tofile((int)output_counter); 
             }

             output_counter+=1.;
		
	     if(logtime==1) { next_output_time *= log_time_factor; }
	     else           { next_output_time = output_counter*output_time + output_time_offset; }
             
             output_chemistry = 1; // Setting this switch will trigger a filling of the reaction rate table at the end of chemistry. The switch is unset afterwards
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
        //    compute_drag_update(1.0*dt) ;

        if(steps >printstuff_steps) {    
                print_velocity_numberdens_ratios(" Pos 0:: ", 210); 
        }
        
        if (do_hydrodynamics == 1) {
        
            //int cnt_broken_cells = 0;
            //
                //March 19th 2024: Added get_cfl_timestep2() to sit here, to obtain the updated timestep based on the extrapolated left and right values - they can produce inconsistent fluxes with the cell-centered values, which the call of dt = get_cfl_timestep(); at the beginning of the timestep is based on;
                //
            //cout<<"num_sepcies = "<<num_species<<endl;
            for(int s = 0; s < num_species; s++) {
                species[s].u_mask           = np_zeros(num_cells+2);
                species[s].u0    = species[s].u ;
                species[s].u_tmp = species[s].u ;
                
                //cout<<"    running species "<<species[s].speciesname<<" s = "<<s<<" steps ="<<steps<<endl;
                 //Apply implicit electron solver for electrons only if so desired. Otherwise continue as usual with all other solvers.
                if( (solver == HydroSolver::implicitelectrons && s==e_idx)  ||   (solver == HydroSolver::implicitelectrons && problem_number==1)  ||  ( (solver == HydroSolver::implicitelectrons) && force_solving_with_implicit)   ) {    
                    double ff = 1.0;
                    if(order == IntegrationType::second_order )
                        ff = 0.5;

                    for(int ii=0; ii<num_implicit_substeps; ii++) {
                        species[s].implicit_incompressible2(ff/((double)num_implicit_substeps)*dt);
                    }

                } else {
                    
                    int ex_order = 1; //(s==e_idx)?0:1;
                    species[s].execute(species[s].u, species[s].dudt[0], species[s].u_mask, ex_order);
                    
                    //species[s].u0 = species[s].u ;
                    for(int j=0; j < num_cells+2; j++)
                        species[s].u_tmp[j] += species[s].dudt[0][j]*dt ;
                    
                    double ekin = 0.5 * species[s].u_tmp[100].u2*species[s].u_tmp[100].u2/species[s].u_tmp[100].u1 ;
                    double expe = ekin + species[s].u_tmp[100].u1 * species[s].cv * 1000.;
                    double expT = (species[s].u_tmp[100].u3-ekin)/species[s].u_tmp[100].u1/species[s].cv;
                    //cout<<species[s].speciesname<<" after +dudt * dt "<<s<<" "<<species[s].dudt[0][100].u3<<" dt = "<<dt<<" u3 = "<<species[s].u_tmp[100].u3<<" expected u3 = "<<expe<<" expected T "<<expT<<endl;

                    species[s].u_mask           = np_zeros(num_cells+2);
                    int numbroken = species[s].count_broken_cells(species[s].u_tmp, species[s].u_mask);
                    //if( numbroken == 0)
                    //    break;
                    
                    //cout<<s<<" before  "<<species[s].u_tmp[100].u3<<endl;

                    species[s].fix_negative_pressures_sometimes(species[s].u_tmp, 1);
                    
                    //cout<<s<<" after "<<species[s].u_tmp[100].u3<<endl;
                    //Done, now all values should be ok
                    for(int j=0; j < num_cells+2; j++) {
                        species[s].u[j] = species[s].u_tmp[j];
                        //species[s].positivity_preserving_step(j);
                    }

                }
                //cout<<"END running species "<<species[s].speciesname<<" s = "<<s<<" steps ="<<steps<<" num_species "<<num_species<<endl;

            }
        }// End of first order hydrodynamic step

        if(steps >printstuff_steps) {    
            print_velocity_numberdens_ratios(" Pos 1:: ", 210);
        }

        //if(use_inflow_damping==1)
        //    apply_inflow_damping();
        
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure(species[s].u);

        
        if(steps >printstuff_steps) {    
                print_velocity_numberdens_ratios(" Pos 1.3:: ", 210);
        }
        
        if (order == IntegrationType::first_order) {
            globalTime += dt;
        } else if (order == IntegrationType::second_order && (1==1)) {
            // 2nd step evaluated at t+dt
            globalTime += dt;

            if(use_self_gravity==1)
                update_mass_and_pot();

            if (do_hydrodynamics == 1) {

                if (use_drag_predictor_step) {
                    for(int s = 0; s < num_species; s++) {
                        for(int j=0; j < num_cells+2; j++) {
                            species[s].u0[j]    = species[s].u[j];
                        }
                         //Those values contain the fixed cells, they might not be identical to u + dudt[0]*dt
                    }
                    
                    for(int s = 0; s < num_species; s++) {
                        species[s].compute_pressure(species[s].u);
			            species[s].fix_negative_pressures_sometimes(species[s].u_tmp, 1);
		            }
                
                    compute_drag_update(0.99*dt) ;
                    
                    if (use_collisional_heating && (use_rad_fluxes==0))
                        compute_collisional_heat_exchange() ; //Disable if radiation is used?
                
                } else
                        compute_drag_update(0.99*dt); //MARCH 28 ONLY FOR DEBUGGING
                
                //if(steps > debug_steps && debug_cell < num_cells+1) {
                if(steps >printstuff_steps) {    
                    print_velocity_numberdens_ratios(" Pos 1.05:: ", 210);
                }
                
                compute_total_pressure();

                //*********************************************Feb 1hth 2025: Anomalous electron temperatures at shocks occur here, after Pos 1.05 and are then perpetuated. Try fix with a bit of friction **//
                //compute_drag_update(0.01*dt);
                //compute_total_pressure();
		        //*******************************************//

                for(int s = 0; s < num_species; s++) {
                    species[s].u_mask           = np_zeros(num_cells+2);
                    species[s].u_tmp = species[s].u ;
                    
                    //Apply implicit electron solver if wanted: either for electrons in applications or for shock tubes
                    if( (solver == HydroSolver::implicitelectrons && s==e_idx)  ||   (solver == HydroSolver::implicitelectrons && problem_number==1)  ||  ( (solver == HydroSolver::implicitelectrons) && force_solving_with_implicit)   ) {    

                        for(int ii=0; ii<num_implicit_substeps; ii++) {
                            species[s].implicit_incompressible2(0.5/((double)num_implicit_substeps)*dt);
                        }
                        

                    } else {
                    //for(int k=0; k<=0; k++) { //The k=0 run is the nominal run. k=1 is only triggered if some cells are broken
                        
                        int ex_order = 1;//(s<=2)? 1.:0; //1;// (s==e_idx)?0:1;
                        species[s].execute(species[s].u, species[s].dudt[1], species[s].u_mask, ex_order);
                        
                        for(int j=0; j < num_cells+2; j++) {
                            double scale_f =  1;//species[s].prim[j].pres/total_press[j];
                            if (use_drag_predictor_step)
                                species[s].u_tmp[j] = species[s].u0[j];// + species[s].dudt[0][j]*dt;// March28th 2024 changed this line, as u0 now contains the first-order correct, non-crashed values
                            
                            species[s].u_tmp[j] +=  (species[s].dudt[1][j] - species[s].dudt[0][j])*dt / 2  * scale_f;  
                        }
                            
                        species[s].u_mask           = np_zeros(num_cells+2);
                        int numbroken = species[s].count_broken_cells(species[s].u_tmp, species[s].u_mask);
                        //if( numbroken == 0)
                        //    break;
                    }
                    
                    species[s].fix_negative_pressures_sometimes(species[s].u_tmp, 2);
                    
                    //Done, now all values should be ok
	                //if(s != e_idx) { //Feb18th: switch off second order update for electrons, as that seems to cause the shock problem
	                if(s > -1) { //Feb18th: switch off second order update for electrons, as that seems to cause the shock problem
	                    for(int j=0; j < num_cells+2; j++) {
                        	species[s].u[j] = species[s].u_tmp[j];
                            //species[s].positivity_preserving_step(j);
                    	}
                    }   
                }
                
                for(int s = 0; s < num_species; s++)
                        species[s].compute_pressure(species[s].u);
                
                for(int s = 0; s < num_species; s++) {
                    species[s].apply_boundary_left(species[s].u) ;
                    species[s].apply_boundary_right(species[s].u) ;
                }

            } else {
                for(int s = 0; s < num_species; s++) {
                    species[s].apply_boundary_left(species[s].u) ;
                    species[s].apply_boundary_right(species[s].u) ;
                }
            }
        }// End of second order hydrodynamic step
        else {
            globalTime += dt;
        }

        if(steps >printstuff_steps) {    
                print_velocity_numberdens_ratios(" Pos 1.4:: ", 210);
        }

        //Damp unphysical inflow at boundaries
        if(use_inflow_damping==1)
            apply_inflow_damping();
        
        //
        // Begin other operators
        //
        for(int s = 0; s < num_species; s++) {
            species[s].compute_pressure(species[s].u);
            species[s].fix_negative_pressures_sometimes(species[s].u, 3);
        }

        if(steps >printstuff_steps) {    
                print_velocity_numberdens_ratios(" Pos 1.5:: ", 210);
        }

        //Misc sources: diffusion:
        if(diffusivity_style >= 0) {

            if(do_hydrodynamics==0){

                vector<double> dummy = {0};
                for(int s=0;s<num_species; s++) 
                    species[s].reconstruct_edge_states(dummy, 1);
            }

            execute_separate_diffusion_step();
        }

        for(int s = 0; s < num_species; s++) {
            species[s].fix_negative_pressures_sometimes(species[s].u, 31);
        }

        //Computes the velocity drag update after the new hydrodynamic state is known for each species
        if (do_hydrodynamics == 1) 
            compute_drag_update(0.99*dt) ;

        if (do_hydrodynamics == 0 && friction_solver > 0) 
                compute_drag_update(0.99*dt) ;

        for(int s = 0; s < num_species; s++) {
            species[s].fix_negative_pressures_sometimes(species[s].u, 4);
        }

        // If either switch is set we need to think more carefully about what should be done
        if( (photochemistry_level + use_rad_fluxes ) > 0 ) {
            
            update_opacities();
            if(photochemistry_level == 0 && use_rad_fluxes > 0)
                reset_dS();

            // Compute high-energy dS and ionization
            if(photochemistry_level == 1) {   //C2Ray scheme
                reset_dS();
                
                if(debug >= 2) {
                    cout<<"Before Photochem dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                do_photochemistry();
                
                if(debug >= 2) {
                    cout<<"After photochem dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                
            } else if(photochemistry_level == 2) { //Linearized time-dependent general chemistry scheme
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
            
            for(int s = 0; s < num_species; s++) {
                species[s].fix_negative_pressures_sometimes(species[s].u, 41);
            }
            
            update_dS();               //Compute low-energy dS
        
            if (do_hydrodynamics == 1) {
                compute_drag_update(0.01*dt) ;
                compute_total_pressure();
            }

            for(int s = 0; s < num_species; s++) {
                species[s].fix_negative_pressures_sometimes(species[s].u, 42);
            }

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

        for(int s = 0; s < num_species; s++) {
            //species[s].compute_pressure(species[s].u);
            species[s].fix_negative_pressures_sometimes(species[s].u, 5);
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
        if(steps >printstuff_steps) {    
                print_velocity_numberdens_ratios(" Pos 3:: ", 210);
        }
        
        if(output_chemistry==1 &&  photochemistry_level==2) {
            write_reaction_table(output_counter-1);
            empty_reaction_table();
            output_chemistry = 0;
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
    if(use_rad_fluxes)
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
        total_press[i] = 0.;
        total_numdens[i] = 0;
            
        for(int s = 0; s < num_species; s++) {
                //total_pressure[i] += species[s].prim[i].pressure;
            //if(species[s].prim_l[i].pres > 0.) {
                total_press[i] += species[s].prim[i].pres;
                total_numdens[i] += species[s].prim[i].number_density ;
                total_press_l[i] += species[s].prim_l[i].pres;
                total_press_r[i] += species[s].prim_r[i].pres;
            //}
                
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
void c_Species::execute(std::vector<AOS>& u_in, std::vector<AOS>& dudt, std::vector<double>& u_mask, int orderstep) {
    
    
        if(base->steps > base->debug_steps && this_species_index == base->debug_species  && base->debug_cell < num_cells+1) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 1 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<endl;
        }
        
        //
        // Step 1: Boundary values
        //
        //

        apply_boundary_left(u_in) ;
        apply_boundary_right(u_in) ;
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species  && base->debug_cell < num_cells+1) {
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
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species && base->debug_cell < num_cells+1) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 3 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
        }
        
        if(debug >= 2)
            cout<<"Done. Starting edge-states."<<endl;
        

        reconstruct_edge_states(u_mask, orderstep) ;
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
        if(base->steps > base->debug_steps && this_species_index == base->debug_species && base->debug_cell < num_cells+1) {
            cout<<"t="<<base->steps<<"IN EXECUTE Pos 4 T[423]_s = [1]= "<<prim[base->debug_cell].temperature<<" p ="<<prim[base->debug_cell].pres<<" eint = "<<prim[base->debug_cell].internal_energy<<" rho = "<<u[base->debug_cell].u1<<" mom = "<<u[base->debug_cell].u2<<" E = "<<u[base->debug_cell].u3<<" manual pressure u_in = "<<(u[base->debug_cell].u3 - 0.5*u[base->debug_cell].u2*u[base->debug_cell].u2/u[base->debug_cell].u1)<<endl;
        }
        
        if(base->steps %1000==0  && this_species_index == -1 ) {    
                    cout<<endl<<"         IN EXECUTE Pos 4:: v["<<20<<"]_s/v_e = "<<endl<<"         ";
                    for(int s = 0; s < base->num_species; s++) {
                        cout<<" ["<<s<<"]= "<<base->species[s].u[20].u2/base->species[s].u[20].u1 / (base->species[2].u[20].u2/base->species[2].u[20].u1);
                    }
                    cout<<endl<<"                           1-n["<<20<<"]_s/n_e = "<<endl<<"         ";
                    for(int s = 0; s < base->num_species; s++) {
                        cout<<" ["<<s<<"]= "<<1-base->species[s].prim[20].number_density/base->species[2].prim[20].number_density;
                    }
                    //cout<<" manual pressure u_in = "<<(species[debug_species].u[debug_cell].u3 - 0.5*species[debug_species].u[debug_cell].u2*species[debug_species].u[debug_cell].u2/species[debug_species].u[debug_cell].u1)<<endl;
                }
        //
        // Step 2: Compute fluxes and sources
        //
	    std::vector<double> grav_prefactors = np_ones(num_cells+2);
        if (not is_dust_like) {
            const double params[3] = {base->mix_p1, base->mix_p2, base->mix_p3};
            
            switch(base->solver) {
                case HydroSolver::hllc:
                    for(int j=0; j <= num_cells; j++) {
                        flux[j] =  get_hlle_flux2(j); //hllc_flux(j);
                        
                    }
                    
                    break;
                case HydroSolver::implicitelectrons:     //Every species which makes it into this loop is not electrons, hence solved with default hllc
                    for(int j=0; j < 3; j++) {
                        flux[j] = hllc_flux(j);
                    }
                    for(int j=3; j <= num_cells; j++) {
                        flux[j] = roe_flux(j);
                    }
                    
                    break;
                case HydroSolver::roe:
                    for(int j=0; j < 3; j++) {
                        flux[j] =  hllc_flux(j);
                    }
                    for(int j=3; j <= num_cells; j++) {
                        flux[j] =  roe_flux(j);
                    }
                    //cout<<"ROE"<<endl;
                    break;
                case  HydroSolver::mix:
                    
                    //if(this_species_index > 0)
                    //    for(int j=0; j <= num_cells; j++) { flux[j]  = roe_flux(j); }
                    //else
                    
                    //cout<<" In mix solver ";
                    //if(this_species_index == base->e_idx) {
                    if(this_species_index >= 0) {
                        for(int j=0; j <= num_cells; j++) {
                            double flim = 1.; //previously 0.1

                    if(this_species_index == 0)
                    flim = params[0];
                    if(this_species_index == base->e_idx)
                                    flim = params[2];
                    else {
                    if(j>base->mix_reset_i)
                                        flim = params[1]; //1e-100; 
                    else 
                        flim = params[1];
                    }
                    //if(j > homopause_boundary_i)
                    //if(j > base->grid2_transition_i)
                    if(base->x_i12[j] > 9e99)
                        flim = 1e-10;

                            double totpress = 0.;
                            int negpresscontributions = 0;
                                        
                            for(int s=0; s<base->num_species; s++) {
                                totpress += base->species[s].prim[j].pres;
                                if(base->species[s].prim[j].pres < 0)
                                    negpresscontributions++;
                            }
                                    //double f           = prim[j].pres/base->total_press[j];
                            double f           = prim[j].pres/totpress;
                                        //f = (f/flim)*(f/flim);
                            //if(f > flim)
                                        //    f=1;

                            f = 1.-1./std::exp( f*f/flim/flim );
                                        //f = std::pow(f, params[0]); //Allow for continuous-linear scaling below threshold
                            f = std::max(f,1e-10);//Cut at very low values to keep sound speeds from rapidly fluctuating
                            //AOS flux1 = hllc_flux(j);
                                        //AOS flux2 = passivescalar_flux2(j);
                                        //flux[j]   = (flux1 * f) +  (flux2 * (1.-f)); 
                            grav_prefactors[j] = f;
                            flux[j] = hllc_flux2(j, f);
                            if(base->steps%1000==0) {
                            //if(base->steps%5000==0 && base->steps > 33000e99) {
                            //if(0==0) {
                                //cout<<" s = "<<this_species_index<<" f ="<<f<<" 1.-f "<<(1.-f)<<" species "<<this_species_index<<" flim "<<flim<<endl;
                                //<<" hllc.u1*f = "<<flux[j].u1<< " "<<flux[j].u1 * f<<" pflux.u1*(1-f) = "<<flux[j].u1<<" "<<flux[j].u1 * (1.-f)<<" 1-flux.u1/hllc.u1 = "<<1.-flux[j].u1/flux1.u1<<" negpress = "<<negpresscontributions<<endl; 
                            }

                            //flux[j]  = hllc_flux(j); // * f + passivescalar_flux(j) * (1.-f); 
                            //flux[j] =  passivescalar_flux(j);
                        }
                    } else {
                        for(int j=0; j <= num_cells; j++) {
                            flux[j] =  hllc_flux(j);
                        }
                    }
                    
                    break;
                case  HydroSolver::laxfriedrich:
                    for(int j=0; j <= num_cells; j++) {
                        flux[j] =  laxfriedrich_flux(j);
                    }
                    break;
                case  HydroSolver::laxwendroff:
                    for(int j=0; j <= num_cells; j++) {
                        flux[j] =  laxwendroff_flux(j);
                    }
                    
                    break;
                break;
            }
                    
 
            if(base->steps % 1000==1 && base->steps < -1002) {
                for(int j=0; j <= num_cells; j++) {
                        AOS roe = roe_flux(j);
                        AOS hllc = hllc_flux(j);
                        AOS p =  laxwendroff_flux(j); //passivescalar_flux(j);
                        if((j>20) && (j<80))
                        //cout<<this_species_index<<" "<<j<<" HLLC = "<<hllc.u1<<" "<<hllc.u2<<" "<<hllc.u3<<" passive ="<<p.u1<<" "<<p.u2<<" "<<p.u3<<" momflux = "<<u[j].u2<<endl;
                            cout<<this_species_index<<" "<<j<<" dP "<<prim_r[j].pres-prim_l[j+1].pres<<" HLLC.u1 "<<hllc.u1<<" HLLC.u1/u.u1 "<<hllc.u1/u[j].u1<<" HLLC.u1/u.u2 "<<hllc.u1/u[j].u2<< " rho "<<u[j].u1<<" HLLC.u1/cs "<<hllc.u1/prim[j].sound_speed<<endl;
                        
                }
            }
            
        }
        else {
            for(int j=0; j <= num_cells; j++)
                flux[j] = dust_flux(j);     
        }
        
        if(debug >= 2)
            cout<<"Done. Starting sources."<<endl;
        
        for(int j=1; j<=num_cells; j++) {
            source[j]          = source_grav(u_in[j], j) * grav_prefactors[j];
            source_pressure[j] = AOS(0, -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
                                          base->source_pressure_prefactor_right[j] * prim_r[j].pres)  ,0); 

            if (this->mass_amu < 0.5)
                source_diffusion[j]= AOS(0,0,0);
            else
                //source_diffusion[j]= ( source_diffusion_flux(j) * base->surf[j] * base->omegaplus[j] * (-1.) + source_diffusion_flux(j-1) * base->surf[j-1] * base->omegaminus[j]) / base->vol[j];
                source_diffusion[j]= ( source_diffusion_flux(j) * base->surf[j] * (-1.) + source_diffusion_flux(j-1) * base->surf[j-1] ) / base->vol[j];
        }
        
        //
        // Step 2a: Compute separate internal energy conserving fluxes for the kinetic energy dominated regime
        //
        //for(int j=base->num_ghosts; j<=num_cells; j++) {
        //     u[j].u4 = prim[j].internal_energy;
        //}
        for(int j=base->num_ghosts; j<=num_cells; j++) {
            compute_e_fluxes(j); //This computes e>0 and if in the kinetic regime, and E-0.5rhou^2 < 0 then it sets new_e = e and updates E, i.e. u[j].u3
        }
        //
        // Step 3: Add it all up to update the conserved variables
        //
        for(int j=base->num_ghosts; j<=num_cells; j++) {
            dudt[j] = (flux[j-1] * base->surf[j-1] - flux[j] * base->surf[j]) / base->vol[j] + (source[j] + source_pressure[j] + source_diffusion[j]) ;
            
            AOS newaos   = (flux[j-1] * base->surf[j-1] - flux[j] * base->surf[j]) / base->vol[j];
            AOS sourceaos = (source[j] + source_pressure[j] + source_diffusion[j]);

            if(j == -100 &&  this->speciesname=="H2" ) {
                cout<<" flux part "<<newaos.u1<<" "<<newaos.u2<<" "<<newaos.u3<<" source part "<<sourceaos.u3<<endl;
                char a;
                cin>>a; 
            }

            //if( debug > 3) { //Or put in your own conditions
            if( (base->steps >= 0) && (j==-100) && (this_species_index==18)) { //Or put in your own conditions    
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
                
                //cin>>alpha;
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

    if( (debug > 2) && (j==100) && (this->speciesname=="H2")) {
        cout<<" IN HLLC, j = "<<j<<" dl/dr = "<<dl<<"/"<<dr<<" ul/ul = "<<ul<<"/"<<ur<<" pl/pr = "<<pl<<"/"<<pr<<" El/Er = "<<El<<"/"<<Er;
    }
    
    //Speed of shocks
    double SL = ul - prim_r[jleft].sound_speed ;
    double SR = ur + prim_l[jright].sound_speed ;
    
    //Intermediate values in the star region, equations 10.30 -10.39 in Toro
    double SS     = ( pr-pl+ mom_l*(SL - ul)-mom_r*(SR-ur) )/(dl*(SL - ul)-dr*(SR-ur) );
    
    double advectionfactor = 1;
    double uavg   = (std::sqrt(u[j].u1) * ul + std::sqrt(u[j].u1) * ur )/(std::sqrt(u[j].u1) + std::sqrt(u[j].u1));
    double hl  = (pl + u[j].u3)/u[j].u1;
    double hr  = (pr + u[j].u3)/u[j].u1;
    double h   = (std::sqrt(u[j].u1) * hl + std::sqrt(u[j].u1) * hr)/(std::sqrt(u[j].u1) + std::sqrt(u[j].u1));
    double c   = std::sqrt((gamma_adiabat-1.) * (h-0.5*uavg*uavg));
    double mach = (std::fabs(uavg)/c);
    
    if ((SL <= 0) &&  (SS >= 0)) {
        AOS FL         = AOS (mom_l, mom_l * ul + pl, ul * (El + pl_e) );
        double comp3_L = El/dl + (SS-ul)*(SS + pl_e/(dl*(SL-ul)));
        AOS US_L       = AOS(advectionfactor ,SS, comp3_L) * dl * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - AOS(dl * advectionfactor, mom_l, El))  * SL;    
           
        flux = FS_L;
        option = 1;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (mom_r, mom_r * ur + pr, ur * (Er + pr_e) );
        double comp3_R = Er/dr + (SS-ur)*(SS + pr_e/(dr*(SR-ur)));
        AOS US_R       = AOS(advectionfactor, SS, comp3_R) * dr * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - AOS(dr * advectionfactor, mom_r, Er)) * SR;
        
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
    if((debug > 2) && (j==100) && (this->speciesname=="H2"))
        cout<<" opt = "<<option<<" SL/SS/SR = "<<SL<<"/"<<SS<<"/"<<SR<<" f = "<<flux.u1<<"/"<<flux.u2<<"/"<<flux.u3<<" ((f.u2-Pl)/Pl)-1 = "<<((flux.u2-pl)/pl)<<" ((f.u2-Pr)/Pr)-1 = "<<((flux.u2-pr)/pr)<<endl;
    return flux;
}

AOS c_Species::hllc_flux2(int j, double f) 
{
    int jleft = j, jright = j+1;
    AOS flux;
    int option = 0;
    
    //Speed of gas
    double ul = prim_r[jleft].speed;  
    double ur = prim_l[jright].speed; 
    
    double pl = f*prim_r[jleft].pres;  
    double pr = f*prim_l[jright].pres;
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
    double SS     = ( pr-pl+ mom_l*(SL - ul)-mom_r*(SR-ur) )/(dl*(SL - ul)-dr*(SR-ur) );
    
    if ((SL <= 0) &&  (SS >= 0)) {
        AOS FL         = AOS (mom_l, mom_l * ul + pl, ul * (El + pl_e) );
        double comp3_L = El/dl + (SS-ul)*(SS + pl_e/(dl*(SL-ul)));
        AOS US_L       = AOS(1 ,SS, comp3_L) * dl * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - AOS(dl, mom_l, El))  * SL;    
           
        flux = FS_L;
        option = 1;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (mom_r, mom_r * ur + pr, ur * (Er + pr_e) );
        double comp3_R = Er/dr + (SS-ur)*(SS + pr_e/(dr*(SR-ur)));
        AOS US_R       = AOS(1, SS, comp3_R) * dr * (SR - ur)/(SR-SS);
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
 * The Lax Friedrichs Flux
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @param[in] dt timestep, as LF needs dt/dx
 * @return flux at cell interface j
 */
AOS c_Species::laxfriedrich_flux(int j) 
{
    int jleft = j, jright = j+1;
    AOS state_l = u[j];
    AOS state_r = u[j+1];
    AOS flux_l  = exact_flux(u[j]);
    AOS flux_r  = exact_flux(u[j+1]);    
    
    double dt  = base->dt;
    double dx1 = 0.5*(base->x_i[j]-base->x_i[j-1]);
    double dx2 = 0.5*(base->x_i[j+1]-base->x_i[j]);
    //double a0  = (0.5*(dx1+dx2)/dt);
    double am1 = std::fabs(prim[j-1].speed) + prim[j-1].sound_speed;
    double am0 = std::fabs(prim[j].speed)   +   prim[j].sound_speed;
    double ap1 = std::fabs(prim[j+1].speed) + prim[j+1].sound_speed;
    double a0 = std::max(am1,std::max(am0,ap1));
    //cout<<" cell j = "<<j<<" fluxes_l ="<<flux_l.u1<<" "<<flux_l.u2<<" "<<flux_l.u3<<" "<<" fluxes_r ="<<flux_r.u1<<" "<<flux_r.u2<<" "<<flux_r.u3<<" "<<endl;
    //cout<<"             "<<" state_l ="<<state_l.u1<<" "<<state_l.u2<<" "<<state_l.u3<<" "<<" fluxes_r ="<<state_r.u1<<" "<<state_r.u2<<" "<<state_r.u3<<" "<<endl;
    


    AOS result = ((flux_l + flux_r) * 0.5) +  ((state_l - state_r) * a0); 
    
    if(j<= 1 || j>=num_cells-1) result = (0,0,0);
    //cout<<"          result  j  "<<j<<" "<<result.u1<<" "<<result.u2<<" "<<result.u3<<endl;
    
    return result;
}

/**
 * The Lax Wendroff Flux
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @param[in] dt timestep, as LF needs dt/dx
 * @return flux at cell interface j
 */
AOS c_Species::laxwendroff_flux(int j) 
{
    int jleft = j, jright = j+1;
    AOS state_l = u[j];
    AOS state_r = u[j+1];
    AOS flux_l  = exact_flux(u[j]);
    AOS flux_r  = exact_flux(u[j+1]);
 
    double dt  = base->dt;
    double dx1 = 0.5*(base->x_i[j]-base->x_i[j-1]);
    double dx2 = 0.5*(base->x_i[j+1]-base->x_i[j]);
    
    AOS state_half = (state_l + state_r) * 0.5 + (flux_l - flux_r) * (0.5*dt/(dx1+dx2));
    AOS result = exact_flux(state_half);
    if(j<= 1 || j>=num_cells-1) result = (0,0,0);
    //cout<<" cell j = "<<j<<" fluxes_l ="<<flux_l.u1<<" "<<flux_l.u2<<" "<<flux_l.u3<<" "<<" fluxes_r ="<<flux_r.u1<<" "<<flux_r.u2<<" "<<flux_r.u3<<" "<<endl;
    //cout<<"             "<<" state_l ="<<state_l.u1<<" "<<state_l.u2<<" "<<state_l.u3<<" "<<" fluxes_r ="<<state_r.u1<<" "<<state_r.u2<<" "<<state_r.u3<<" "<<endl;
    //cout<<"         LW result  j  "<<j<<" "<<result.u1<<" "<<result.u2<<" "<<result.u3<<endl;
    
    return result;
}

/**
 * The Roe Flux
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @param[in] dt timestep, as LF needs dt/dx
 * @return flux at cell interface j
 */
AOS c_Species::roe_flux(int j) 
{
    int jleft = j, jright = j+1;
    AOS_prim prim_l  = this->prim_r[jleft];
    AOS_prim prim_r  = this->prim_l[jright];
    
    AOS state_l      = AOS(prim_l.density, prim_l.speed*prim_l.density, prim_l.density*prim_l.internal_energy + 0.5*prim_l.density*prim_l.speed*prim_l.speed); //u[jleft];
    AOS state_r      = AOS(prim_r.density, prim_r.speed*prim_r.density, prim_r.density*prim_r.internal_energy + 0.5*prim_r.density*prim_r.speed*prim_r.speed); // = u[jright];
    
    AOS flux_l       = exact_flux(state_l);
    AOS flux_r       = exact_flux(state_r);
    
    double drho = state_r.u1 - state_l.u1;
    double dp   = prim_r.pres - prim_l.pres;
    double du   = prim_r.speed - prim_l.speed;
    
    //Roe averages
    double rho = std::sqrt(state_l.u1*state_r.u1);
    double u   = (std::sqrt(state_l.u1) * prim_l.speed + std::sqrt(state_r.u1) * prim_r.speed )/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double hl  = (prim_l.pres + state_l.u3)/state_l.u1;
    double hr  = (prim_r.pres + state_r.u3)/state_r.u1;
    double h   = (std::sqrt(state_l.u1) * hl + std::sqrt(state_r.u1) * hr)/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double c   = std::sqrt((gamma_adiabat-1.) * (h-0.5*u*u));
 
    //double lambdas[3] = {std::abs(u-c), u, std::abs(u+c) };
    double lambdas[3] = {std::fabs(u-c), std::fabs(u), std::fabs(u+c) };
    double eps = 5e-1;
    for(int k=0; k<=2; k+=2)
        if(lambdas[k]/c<eps)
            //lambdas[k] = 0.5*(lambdas[k]*lambdas[k]/eps+eps);
            lambdas[k] = 0.5*(lambdas[k]*lambdas[k]/(c*eps)+c*eps);
    
    double d[3] = {drho, state_r.u2 - state_l.u2, state_r.u3 - state_l.u3};
    double alphas[3];
    
    alphas[0] = (dp - rho * c * du)/(2*c*c);
    alphas[1] = drho - dp/(c*c);
    alphas[2] = (dp + rho * c * du)/(2*c*c);
    
    if(j==-80 && base->steps<10)
        cout<<" alphas "<<alphas[0]<<" "<<alphas[1]<<" "<<alphas[2]<<endl;
    
    AOS ev0 = AOS(1, u-c, h - u*c);
    AOS ev1 = AOS(1, u, 0.5*u*u);
    AOS ev2 = AOS(1, u+c, h + u*c);
    
    AOS result = (flux_l + flux_r) * 0.5 - ( (ev0*(alphas[0]*lambdas[0]))  + (ev1*(alphas[1]*lambdas[1])) + (ev2*(alphas[2]*lambdas[2]))) * 0.5;
    
    return result;
}

Vector3d c_Species::roe_flux_vec(int j) {
    AOS tmp    = roe_flux(j);

    return Vector3d(tmp.u1, tmp.u2, tmp.u3);
}

Vector3d c_Species::u_to_vec(int j) {
    Vector3d v = {u[j].u1, u[j].u2, u[j].u3}; 
    return v;
}

//
// returns Roe averages at interface j - sits between cell j and j+1
//
AOS c_Species::get_roe_averages(int j) {
 int jleft = j, jright = j+1;
    AOS_prim prim_l  = this->prim[jleft];
    AOS_prim prim_r  = this->prim[jright];
    
    AOS state_l      = AOS(prim_l.density, prim_l.speed*prim_l.density, prim_l.density*prim_l.internal_energy + 0.5*prim_l.density*prim_l.speed*prim_l.speed); //u[jleft];
    AOS state_r      = AOS(prim_r.density, prim_r.speed*prim_r.density, prim_r.density*prim_r.internal_energy + 0.5*prim_r.density*prim_r.speed*prim_r.speed); // = u[jright];
    
    AOS flux_l       = exact_flux(state_l);
    AOS flux_r       = exact_flux(state_r);
    
    double drho = state_r.u1 - state_l.u1;
    double dp   = prim_r.pres - prim_l.pres;
    double du   = prim_r.speed - prim_l.speed;
    
    //Roe averages
    double rho = std::sqrt(state_l.u1*state_r.u1);
    double u   = (std::sqrt(state_l.u1) * prim_l.speed + std::sqrt(state_r.u1) * prim_r.speed )/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double hl  = (prim_l.pres + state_l.u3)/state_l.u1;
    double hr  = (prim_r.pres + state_r.u3)/state_r.u1;
    double h   = (std::sqrt(state_l.u1) * hl + std::sqrt(state_r.u1) * hr)/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double c   = std::sqrt((gamma_adiabat-1.) * (h-0.5*u*u));

    return AOS(rho, u, h);
}

/**
 * A Flux for passive scalars - derived from Roe
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @param[in] dt timestep, as LF needs dt/dx
 * @return flux at cell interface j
 */
AOS c_Species::passivescalar_flux(int j) 
{
    
    int jleft = j, jright = j+1;
    AOS state_l = u[j];
    AOS state_r = u[j+1];
    AOS flux_l  = exact_advection_flux(u[j]);
    AOS flux_r  = exact_advection_flux(u[j+1]);
 
    double dt  = base->dt;
    double dx1 = 0.5*(base->x_i[j]-base->x_i[j-1]);
    double dx2 = 0.5*(base->x_i[j+1]-base->x_i[j]);
    
    AOS state_half = (state_l + state_r) * 0.5 + (flux_l - flux_r) * (0.5*dt/(dx1+dx2));
    AOS result     = exact_advection_flux(state_half);
    //double speed = state_half.u2/state_half.u1;
    //AOS result     = AOS(state_half.u2, state_half.u2*speed, state_half.u3*speed); //Assuming the exact, passive scalar flux
    
    if(j<= 1 || j>=num_cells-1) result = (0,0,0);
    //cout<<" cell j = "<<j<<" fluxes_l ="<<flux_l.u1<<" "<<flux_l.u2<<" "<<flux_l.u3<<" "<<" fluxes_r ="<<flux_r.u1<<" "<<flux_r.u2<<" "<<flux_r.u3<<" "<<endl;
    //cout<<"             "<<" state_l ="<<state_l.u1<<" "<<state_l.u2<<" "<<state_l.u3<<" "<<" fluxes_r ="<<state_r.u1<<" "<<state_r.u2<<" "<<state_r.u3<<" "<<endl;
    //cout<<"         LW result  j  "<<j<<" "<<result.u1<<" "<<result.u2<<" "<<result.u3<<endl;
    
    return result;
}



/**
 * Passive scalar flux 2 - derived from HLLC
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @return flux at cell interface j
 */
AOS c_Species::passivescalar_flux2(int j) 
{
 int jleft = j, jright = j+1;
    AOS state_l = u[j];
    AOS state_r = u[j+1];
    AOS flux_l  = exact_advection_flux(u[j]);
    AOS flux_r  = exact_advection_flux(u[j+1]);
 
    double uhalf = 0.5*(state_l.u2/state_l.u1 + state_r.u2/state_r.u1);
    AOS result = AOS(0,0,0);
    if(uhalf > 0)
        result = AOS(state_l.u1, state_l.u2, state_l.u3) * uhalf;
    else if(uhalf < 0)
        result = AOS(state_r.u1, state_r.u2, state_r.u3) * uhalf;
    
    if(j<= 1 || j>=num_cells-1) result = AOS(0,0,0);
    
    return result;
}

/**
 * The Exact Flux, prominent ingredient in Riemann-free flux-conservative methods
 * 
 * @param[in] u the hydrodynamic state u=(rho, momentu, total energy), of which to compute the flux from
 * @param[in] prim the primitive corresponding to u. Needed to compute the exact flux
 * @return flux at cell interface j
 */
AOS c_Species::exact_flux(AOS u) 
{
    double speed = u.u2/u.u1;
    double press = (gamma_adiabat - 1.) *(u.u3 - 0.5 * u.u2 * speed);
    double flux2 = (press + u.u2*speed);
    double flux3 = (speed*(u.u3 + press));
    AOS result = AOS(u.u2, flux2, flux3);
    
    return result;
}

Vector3d c_Species::exact_flux_as_vector(int j) 
{
    AOS uu = this->u[j];
    double speed = uu.u2/uu.u1;
    double press = (gamma_adiabat - 1.) *(uu.u3 - 0.5 * uu.u2 * speed);
    double flux2 = (press + uu.u2*speed);
    double flux3 = (speed*(uu.u3 + press));
    Vector3d result(uu.u2, flux2, flux3);
    
    return result;
}

AOS c_Species::exact_advection_flux(AOS u) 
{
    double speed = u.u2/u.u1;
    AOS result = AOS(u.u2, u.u2*speed, u.u3*speed);
    return result;
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
 * @param[in] debuginfotime timestep number from which detailed info on the internal workings is printed
 * @param[out] u_mask Output mask values. Any offending cell is set from 0 to 1
 * @return the sum of all u_mask values. The main loop will react if non-zero values are returned.
 */
int c_Species::count_broken_cells(std::vector<AOS>&u, std::vector<double>&u_mask) {
    
    int cnt =0;
    std::vector<double> p_temps = np_somevalue(num_cells+2, 0.);
    
    //Predict pressures from this hydro data
    for(int j=0; j<=num_cells+1; j++) {
            p_temps[j] = (u[j].u3 - 0.5*u[j].u2*u[j].u2/u[j].u1)  * (gamma_adiabat-1);
        }
    for(int j=0; j<=num_cells+1; j++) {
            //Cheap -E detection and -p detection
            u_mask[j] += std::signbit( u[j].u3 ) + std::signbit( p_temps[j] );
            if(u_mask[j] > 0.5)
                u_mask[j] = 1.;
            if( std::signbit( u[j].u3 ) + std::signbit( p_temps[j] ) > 0)
                cnt++;
            if(u_mask[j] > 0.5) {
                /*
                if(j+1<=num_cells+1)
                    u_mask[j+1] += 1.*u_mask[j];
                if(j-1>=0)
                    u_mask[j-1] += 1.*u_mask[j]; 

                */
            }
        }
        
        //if(base->steps >= base->debug_steps && this_species_index==1) {
        if(base->steps >= 585e99 && this_species_index==1 && cnt>0) {
            cout<<" IN COUNT BROKEN CELLS, steps"<<base->steps<<" manual pressures "<<endl;
            for(int j=0; j<=num_cells+2; j++) {
               if(u_mask[j] > 0.5)
                   cout<<"           p["<<j<<"] = "<< u[j].u3 - 0.5*u[j].u2*u[j].u2/u[j].u1<<" E = "<<u[j].u3<<" Ekin/E"<<0.5*u[j].u2*u[j].u2/u[j].u1/u[j].u3<<endl;
               
            }
            cout<<" and cntresult = "<<cnt<<endl;
            
        }


        return cnt;
}

/**
 * Fixes negative pressures in the cases where Ekin > E, by taking the last valid, positive value of the pressure and assigning it as new pressure. The total Energy is then updated.
 * The name of this function "sometimes" refers to that we should only fix negative pressures, if the causes are well understood, as we might violate energy conservation.
 * 
 */
int c_Species::fix_negative_pressures_sometimes(std::vector<AOS>&u_temp, int flag) {
    
    //if(this_species_index==base->e_idx)
    //    cout<<" electron temperature, steps "<<base->steps<<" fl "<<flag<<" T ="<<prim[100].temperature<<endl; 

    double ptemp = 0;
    double Ttemp = 0;
    double plast = 0;
    double eratio = 0;
    int fixed_cells = 0;
    double ekin = 0;
    double ekinold=0;
    double temper=0;
    double Tfix = 0;
    for(int j=0; j<=num_cells+1; j++) {

        int fixed = 0;
        ekin   = 0.5*u_temp[j].u2*u_temp[j].u2/u_temp[j].u1;
        ekinold= 0.5*u[j].u2*u[j].u2/u[j].u1;
        eratio = ekin/u_temp[j].u3;
        ptemp  = (u_temp[j].u3 - ekin) * (gamma_adiabat-1);
        Ttemp  = (u_temp[j].u3 - ekin)/u_temp[j].u1/this->cv;
        plast  = prim[j].pres;
        temper = prim[j].temperature;
        

        if( (u_temp[j].u1 < 0)) {    
            double rtmp = prim[j].density;
            if( (j>1) && (j<num_cells))
                rtmp = 0.5*std::log10(u[j-1].u1) + 0.5*std::log10(u[j+1].u1) ;
            
            u_temp[j].u1 = std::pow(10.,rtmp);
            u_temp[j].u2 = prim[j].speed * u_temp[j].u1;
            ekin = 0.5*u_temp[j].u2*u_temp[j].u2/u_temp[j].u1;
            u_temp[j].u3 = ekin + u_temp[j].u1 * prim[j].internal_energy;

            
            eos->compute_primitive(&(u_temp[j]), &(prim[j]), 1) ;    
            eos->compute_auxillary(&(prim[j]), 1);

            fixed = 1;
        }
        //First check: Negative E - big problem - use last value and fix pressure in case its also negative 
        if( 1==0) {    
            //u_temp[j].u3 = u[j].u3;// + (ekinold-ekin);

            //if(this->speciesname=="H2" && j==10)
            //    cout<<"Repaired E in cell/species = "<<j<<" "<<this->speciesname<<" steps "<<base->steps<<" flag "<<flag<<"  Es "<<u_temp[j-1].u3<<" "<<u_temp[j].u3<<" "<<u_temp[j+1].u3<<"  rho "<<u_temp[j-1].u1<<" "<<u_temp[j].u1<<" "<<u_temp[j+1].u1<<" Ttemp "<<Ttemp<<"  v "<<u_temp[j-1].u2/u_temp[j-1].u1<<" "<<u_temp[j].u2/u_temp[j].u1<<" "<<u_temp[j+1].u2/u_temp[j+1].u1<<endl;
            
            double etmp = ekin;
            if( (j>1) && (j<num_cells))
                etmp = 0.5*std::log10(u[j-1].u3) + 0.5*std::log10(u[j+1].u3) ;
            
            u_temp[j].u3 = std::pow(10.,etmp);
            
            eos->compute_primitive(&(u_temp[j]), &(prim[j]), 1) ;    
            eos->compute_auxillary(&(prim[j]), 1);

            fixed = 1;
        }
        //Second check: negative pressure
        if((Ttemp < 0) || (temper < 0) || (temper > base->max_temperature) || std::isnan(temper)  ||  (u_temp[j].u3 < 0) )  { //
            base->update_T_mean(j, flag);
            double T_mean = base->T_mean[j];
            if(std::isnan(T_mean) || (T_mean<0))
                cout<<" in repair T, T_mean is "<<T_mean<<endl;

            //if((this->speciesname=="H2") && j==10) {
            //if(temper > base->max_temperature)
            //cout<<"Repaired T in cell/species = "<<j<<" "<<this->speciesname<<" steps "<<base->steps<<" implied temp "<<Ttemp<<" current temper "<<temper<<" T_mean "<<T_mean<<" dt "<<base->dt<<" maxT = "<<base->max_temperature<<endl;
                /* cout<<"conditions "<<(Ttemp < 0)<<(temper < 0)<<std::isnan(temper) << (temper > base->max_temperature)<<endl;
                cout<<((Ttemp < 0) || (temper < 0) || std::isnan(temper) || (temper > base->max_temperature) )<<endl;
                if((Ttemp < 0) || (temper < 0) || std::isnan(temper) || (temper > base->max_temperature) )
                    cout<<" total expression true "<<endl;
                else
                    cout<<" total expression false "<<endl;

                if( 0 )
                    cout<<" 0 expression true "<<endl;
                else
                    cout<<" 0 expression false "<<endl;
            //} */
            //cout<<"Repaired T in cell/species = "<<j<<" "<<this->speciesname<<" steps "<<base->steps<<" flag "<<flag<<"  Es "<<u_temp[j-1].u3<<" "<<u_temp[j].u3<<" "<<u_temp[j+1].u3<<"  rho "<<u_temp[j-1].u1<<" "<<u_temp[j].u1<<" "<<u_temp[j+1].u1<<" Ttemp "<<Ttemp<<"  v "<<u_temp[j-1].u2/u_temp[j-1].u1<<" "<<u_temp[j].u2/u_temp[j].u1<<" "<<u_temp[j+1].u2/u_temp[j+1].u1<<endl;

            double T_tmp = std::min(std::max(T_mean, base->temperature_floor), base->max_temperature);
            volatile double enew  = u_temp[j].u1*this->cv*T_tmp;
            u_temp[j].u3 = enew + ekin;

            eos->compute_primitive(&(u_temp[j]), &(prim[j]), 1) ;    
            eos->compute_auxillary(&(prim[j]), 1);
            fixed = 1;
        }
        else{
            double M;
            double ret = return_entropy_with_jump(M);
        }
        
    }
     
    return 0;//fixed_cells;
}

//
// A routine that prints a lot of stuff
//
void c_Sim::print_velocity_numberdens_ratios(string position, int dcell)
{
    if(0==1) {
        cout<<endl<<position<<" E["<<dcell<<"]_s = "<<endl<<"           ";
                    for(int s = 0; s < num_species; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].u[dcell].u3;
                    }
        cout<<endl<<position<<" v["<<dcell<<"]_s/v_e = "<<endl<<"           ";
                        for(int s = 0; s < num_species; s++) {
                            cout<<std::setprecision (15)<<" ["<<s<<"]= "<<(species[s].u[dcell].u2/species[s].u[dcell].u1) / (species[2].u[dcell].u2/species[2].u[dcell].u1);
                        }
        cout<<endl<<position<<" 1-n["<<dcell<<"]_s/n_e = "<<endl<<"           ";
                        for(int s = 0; s < num_species; s++) {
                            cout<<std::setprecision (15)<<" ["<<s<<"]= "<<1.-(species[s].u[dcell].u1/species[s].mass_amu) / (species[2].u[dcell].u1/species[2].mass_amu);
                        }
                        
        cout<<endl<<position<<" T["<<dcell<<"]_s = "<<endl<<"           ";
                    for(int s = 0; s < num_species; s++) {
                        cout<<" ["<<s<<"]= "<<species[s].prim[dcell].temperature;
                    }
                    cout<<endl<<"_____________________________________________________"<<endl;
    }
}

//
// Damp negative velocities under certain conditions
//
void c_Sim::apply_inflow_damping() {

    if(globalTime > 1e-10) {
        double damping_time = 1e-2;

        for (int s=0; s<num_species; s++) {
            for(int j=1; j<=3; j++) {

                //if(species[s].u[j].u2 < 0) {
                    
                    species[s].u[j].u2 /= (1+dt/damping_time);
                    species[s].u[j].u3 = species[s].prim[j].pres/(species[s].gamma_adiabat-1.) + 0.5 * species[s].u[j].u2 * species[s].u[j].u2/species[s].u[j].u1;
                //}
            }
            //Recompute primitive variables after damping
            species[s].compute_pressure(species[s].u);
        }

    }


}

//Compute explicit fluxes for small e, to hopefully preserve positivity
void c_Species::compute_e_fluxes(int j) { 

}

//
// Dual energy hydro step for the internal energy
//
void c_Species::positivity_preserving_step(int j) { //From Zhang & Shu 2011

    double e_pred = u[j].u3 - 0.5 * u[j].u2 * u[j].u2/u[j].u1;
    if(e_pred < 0) {
        double eold = prim[j].internal_energy;
        double enew = 1;

        double ul = 0.5*(prim[j-1].speed + prim[j].speed);
        double ur = 0.5*(prim[j].speed + prim[j+1].speed);
        double fac = base->dt * (gamma_adiabat-1.5);

        enew = eold / (1 +fac * ( ul* base->surf[j-1] - ur * base->surf[j]) / base->vol[j] ) ; //Implicit, kinetic energy dominated step
        
        u[j].u3 = enew * u[j].u1 + 0.5 * u[j].u2 * u[j].u2/u[j].u1;
        
        if(enew < 0)
            cout<<speciesname<<" j = "<<j<<" in positivity preserving step, but e<0 !"<<endl;
    }
}


//
// Returns the dominant entropy jump across from left or right cell
//
double c_Species::return_entropy_with_jump(double k) {
/* 
    double g = gamma_adiabat;
    double dentropy[2] = {0,0};
    double entropy_l  = prim[k-1].pres/std::pow(prim[k-1].density, g);
    double entropy_r  = prim[k+1].pres/std::pow(prim[k+1].density, g);
    
    //Do both interfaces for cell k
    for(int j=0; j<=1; j++) {

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
        double SS     = ( pr-pl+ mom_l*(SL - ul)-mom_r*(SR-ur) )/(dl*(SL - ul)-dr*(SR-ur) );
        
        double advectionfactor = 1;
        double uavg   = (std::sqrt(u[j].u1) * ul + std::sqrt(u[j].u1) * ur )/(std::sqrt(u[j].u1) + std::sqrt(u[j].u1));
        double hl  = (pl + u[j].u3)/u[j].u1;
        double hr  = (pr + u[j].u3)/u[j].u1;
        double h   = (std::sqrt(u[j].u1) * hl + std::sqrt(u[j].u1) * hr)/(std::sqrt(u[j].u1) + std::sqrt(u[j].u1));
        double c   = std::sqrt((gamma_adiabat-1.) * (h-0.5*uavg*uavg));
        double mach = (std::fabs(uavg)/c);
        
        if ((SL <= 0) &&  (SS >= 0)) {

        }
        double SL = ul - prim_r[jleft].sound_speed ;
        double SR = ur + prim_l[jright].sound_speed ;

        double u_up, c_up;
        if( SS > 0 )

        double M = std::fabs(u_up-SL)/c_up;

        double nom   = 1+(2*g/(1+g))*(M*M-1);
        double denom = std:pow((g+1)*M*M/((g-1)*M*M + 2.), g)
        double dentropy[j] = cv*std::log(nom/denom);
    }
 */
    return 0 ;
}


/**
 * Core function for the HLLE approach - feeds wave speed estimates and left and right states into both explicit and implicit routines.
 */
void c_Species::write_hlle_wavespeeds(AOS_prim prim_l_in, AOS_prim prim_r_in, AOS &state_l_out, AOS &state_r_out, double &SL_out, double &SR_out) {

    AOS_prim prim_l  = prim_l_in;
    AOS_prim prim_r  = prim_r_in;
    
    //Speed of gas
    double ul = prim_l.speed;  
    double ur = prim_r.speed; 
    double cl = prim_l.sound_speed;  
    double cr = prim_r.sound_speed;
    double pl = prim_l.pres;  
    double pr = prim_r.pres;
    double pl_e = pl;
    double pr_e = pr;
    double dl = prim_l.density;  
    double dr = prim_r.density;
    double mom_l = dl*ul ;
    double mom_r = dr*ur ;
    double El = dl*prim_l.internal_energy + 0.5*mom_l*ul ;
    double Er = dr*prim_r.internal_energy + 0.5*mom_r*ur ;

    AOS state_l      = AOS(prim_l.density, mom_l, El); //u[jleft];
    AOS state_r      = AOS(prim_r.density, mom_r, Er); // = u[jright];

    //Roe averages
    double denom = (std::sqrt(dl) + std::sqrt(dr));
    double rho   = std::sqrt(dl*dr);
    double v     = (std::sqrt(dl) * ul + std::sqrt(dr) * ur )/denom;
    double hl    = (pl + state_l.u3)/dl;
    double hr    = (pr + state_r.u3)/dr;
    double h     = (std::sqrt(dl) * hl + std::sqrt(dr) * hr)/denom;
    double c     = std::sqrt((gamma_adiabat-1.) * (h-0.5*v*v));
    
    //Davis estimates
    /* SL_out = std::min(ul - cl, ur-cr) ;
    SR_out = std::max(ul + cl, ur+cr) ; */
    /*
    SL_out = v - c ;
    SR_out = v + c ;
    */
    /* 
    //Einfeld estimates
    double eta2 = 0.5*(std::sqrt(dl*dr)/(denom*denom));
    double dbar = (std::sqrt(dl) * cl*cl + std::sqrt(dr) * cr*cr )/denom + eta2 * (ur-ul)*(ur-ul);
    SL_out = v - dbar;
    SR_out = v + dbar;
 */
    SL_out = v - c ;
    SR_out = v + c ; //Those occasionally violate the consistency condition F_num(U,U) = F_exact(U)
    SL_out = std::min(ul - cl, ur-cr) ;
    SR_out = std::max(ul + cl, ur+cr) ; 

    double eta2 = 0.5*(std::sqrt(dl*dr)/(denom*denom));
    double dbar = (std::sqrt(dl) * cl*cl + std::sqrt(dr) * cr*cr )/denom + eta2 * (ur-ul)*(ur-ul);
    SL_out = v - dbar;
    SR_out = v + dbar;

    state_l_out = state_l;
    state_r_out = state_r;
}

/**
 * The HLLE Riemann solver with pressure-based wave-speed estimate, according to Toro(2007) and references therein.
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @return flux at cell interface j
 */
AOS c_Species::get_hlle_flux2(int j) 
{
    int jleft = j, jright = j+1;

    //if(base->steps>100)
    //    jright = 1;

    AOS flux, state_l, state_r;
    int option = 0;
    double SL, SR, denom;
    AOS_prim prim_l  = this->prim_r[jleft];
    AOS_prim prim_r  = this->prim_l[jright];

    write_hlle_wavespeeds(prim_l, prim_r, state_l, state_r, SL, SR); //in, in, out,out,out,out
    AOS FL = exact_flux(state_l);
    AOS FR = exact_flux(state_r);
    denom = (SR-SL);

    if ((SL <= 0) &&  (SR >= 0)) {
        flux = (FL * SR - FR * SL + (state_r - state_l)*SL*SR)/denom;
        option = 1;
    } else if (SL >= 0) {
        flux = FL;
        option = 2;
    } else {
        flux = FR;
        option= 3 ;
    }
    
    //if(base->steps>100)
    //    cout<<" checking consistency, f_numeric(U,U) =  "<<flux.u1<<" "<<flux.u2<<" "<<flux.u3<<" "<<" f_exact(U) = "<< FL.u1<<" "<< FL.u2<<" "<< FL.u3<<" "<<endl;

    return flux;
}

/*
* Write Jacobian matrices consistent with our HLLE implementation
*/
void c_Species::write_hlle_jacobians(Matrix3d &left_m, Matrix3d &right_m, int j) {

    int jleft = j, jright = j+1;
    AOS flux, state_l, state_r;
    int option = 0;
    double SL, SR, denom;
    AOS_prim prim_l  = this->prim_r[jleft];
    AOS_prim prim_r  = this->prim_l[jright];

    write_hlle_wavespeeds(prim_l, prim_r, state_l, state_r, SL, SR); //in, in, out,out,out,out
    AOS FL = exact_flux(state_l);
    AOS FR = exact_flux(state_r);
    denom = (SR-SL);

    Matrix3d A_left  = get_exact_Jacobian(state_l);
    Matrix3d A_right = get_exact_Jacobian(state_r);

    Matrix3d uni;
    uni<<1,0,0,  0,1,0, 0,0,1;
    Matrix3d A_tilda   = uni*(SL*SR/denom);
    //cout<<" returnting Atilda and Atilda2: "<<endl<<A_tilda<<endl<<endl<<A_tilda2<<endl;
 
    left_m  = A_left * (SR/denom) - A_tilda;
    right_m = A_right*(-SL/denom) + A_tilda;
 
    //left_m  = A_left;
    //right_m = A_right;
}


/**
 * A wrapper for he HLLE Riemann solver with pressure-based wave-speed estimate, according to Toro(2007) and references therein.
 * 
 * @param[in] j cell interface number at which to compute the flux. 
 * @return flux at cell interface j
 */
Vector3d c_Species::get_hlle_flux(int j) 
{
    AOS flux = get_hlle_flux2(j) ;
    Vector3d vflux;
    vflux << flux.u1, flux.u2, flux.u3;
    return vflux;
}


void c_Species::write_exact_jacobians(Matrix3d &left_m, Matrix3d &right_m, int j) {
    left_m  = get_exact_Jacobian(u[j]);
    right_m = get_exact_Jacobian(u[j+1]);
}
///////////////////////////////////////////////////////////
//
//
//  advection.cpp
//
// This file contains the base routines to start and run a simulation, 
// as well as hydrodynamics, i.e. the HLLC Riemann solver.
//
//
//
///////////////////////////////////////////////////////////

#include "aiolos.h"


////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  CLASS SIMULATION
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void c_Sim::execute() { 
     
    steps = 0;
    double output_counter = 0;
    double monitor_counter= 0;
    const double dt_initial = dt_min_init;
    double next_print_time  = dt_initial;
    const int maxsteps = 1e9;
    
    int crashed_T = 0, crashed_J = 0;
    int crash_T_imin = num_cells+2, crash_T_imax = 0, crash_T_numcells = 0;
    int crash_J_imin = num_cells+2, crash_J_imax = 0, crash_J_numcells = 0;
    double crashtime, crashed_temperature, crashed_meanintensity;
    int crashtime_already_assigned = 0;
        
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
    //cout<<" FIrst species mass: "<<species[0].mass_amu<<endl;
//     cout<<" "<<endl;
//     cout<<" Initiating malygin opacities"<<endl;
//     init_malygin_opacities();
//     cout<<" pressure_manual = "<<species[0].prim[3].density * species[0].prim[3].temperature * kb / (40.*amu)<<" pressure_true = "<<species[0].prim[3].pres <<endl;
//     cout<<endl;
//    kappa_landscape();
//     cout<<" test planck opacity = "<<opacity_semenov_malygin(0,    species[0].prim[3].temperature, species[0].prim[3].density, species[0].prim[3].pres);
//     cout<<" test rosseland opacity = "<<opacity_semenov_malygin(1, species[0].prim[3].temperature, species[0].prim[3].density, species[0].prim[3].pres);
//     
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {
        
        //if(steps > 17964) 
        //    cout<<" steps "<<steps<<endl;
        
        //if(steps > 18024) 
        //    debug = 2;
        
        if(steps==0) {
            for(int s = 0; s < num_species; s++)
                species[s].compute_pressure(species[s].u);
            compute_total_pressure();
        }
            
            
        dt = get_cfl_timestep();
        dt = std::min(dt, t_max - globalTime) ;
        if(steps == 0)
            dt = std::min(dt, dt_initial);
        
//         if(steps == 1e6) 
//             cout<<"Radiative equilibrium phase over."<<endl;
        
        if( globalTime > next_print_time) {
            cout<<" Beginning step "<<steps<<" @ globalTime "<<globalTime<<" dt "<<dt;
            cout<< ", CFL " << cfl_step << ", radiative dt " << timestep_rad2 << "\n";
            next_print_time *= 10.;
        }
            
//         
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
        //init_grav_pot();
        if(debug >= 2)
            cout<<"Beginning timestep "<<steps<<endl;
        
        //if(use_self_gravity==1)
        update_mass_and_pot();
        
        if(debug >= 2)
            cout<<"Before fluxes... ";
        
        //
        // Output data, when required. Keep the output here, between the flux and conserved variable update, so that the 
        // zeroth output has the initialized conserved values, but already the first fluxes.
        //
        
        if(steps==0 || ((steps==1 || steps==2) && debug > 0)) {
            for(int s=0; s<num_species; s++) {
                species[s].print_AOS_component_tofile((int) output_counter);
            }
            print_monitor((int)monitor_counter);
            //if(steps==1 || steps==2)
            print_diagnostic_file((int)output_counter);
            
            user_output_function(output_counter) ;

            monitor_counter+=1.;
            output_counter +=1.;
         }
         if(globalTime > output_counter*output_time){
             if(debug >= 1)
                 cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
             
             print_diagnostic_file((int)output_counter);
             for(int s=0; s<num_species; s++)
                species[s].print_AOS_component_tofile((int)output_counter); 
             
             user_output_function(output_counter) ;
             
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
        
        
        //if(steps >= 4345)
        //    debug = 2;
        
        
        //
        // Step 0: If no gas movement is desired, set all velocity changes to 0
        //
        
        if (do_hydrodynamics == 1) {
        
            for(int s = 0; s < num_species; s++)
                species[s].execute(species[s].u, species[s].dudt[0]);
            
            for(int s=0; s < num_species; s++) {
                species[s].u0 = species[s].u ;
                for(int j=0; j < num_cells+2; j++)
                    species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
            }
        }
        
        if (order == IntegrationType::first_order) {
            globalTime += dt;
        }
        else if (order == IntegrationType::second_order) {
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
                        compute_collisional_heat_exchange() ;
                }

                for(int s = 0; s < num_species; s++)
                    species[s].execute(species[s].u, species[s].dudt[1]);

                for(int s = 0; s < num_species; s++){
                    for(int j = 0; j < num_cells+2; j++) {
                        // Need to re-compute u[j] if it was updated in the drag-predictor
                        if (use_drag_predictor_step)
                            species[s].u[j] = species[s].u0[j] + species[s].dudt[0][j]*dt ;
                        
                        species[s].u[j] += (species[s].dudt[1][j] - species[s].dudt[0][j])*dt / 2 ;  
                    }
                }
            } else {
                for(int s = 0; s < num_species; s++) {
                    species[s].apply_boundary_left(species[s].u) ;
                    species[s].apply_boundary_right(species[s].u) ;
                }
            }
        }
        
        for(int s = 0; s < num_species; s++) 
            species[s].compute_pressure(species[s].u);
            
        //compute_total_pressure();
        
        if (do_hydrodynamics == 1)
            compute_drag_update() ;
        
        if( (photochemistry_level + use_rad_fluxes ) > 0 ) {
            
            update_opacities();
            reset_dS();
            
            // Compute high-energy dS and ionization
            if(photochemistry_level > 0)
                do_photochemistry();
            
            update_dS();               //Compute low-energy dS

            // Apply any additional heating / cooling
            user_heating_function() ;
            
            if(use_rad_fluxes==1) {
                update_fluxes_FLD();   //FLD Radiation transport, updating Temperatures and photon band energies
            }
            else {
                //update_temperatures_simple();//Fast eulerian temperature update, no stability guaranteed!
                //if (use_collisional_heating)
                //    compute_collisional_heat_exchange();
            }
        }

        user_loop_function() ;
        
        steps++;
        
        if(steps==1)
            cout<<"Initial sound crossing time = "<<max_snd_crs_time<<endl;
        
        if(debug > 1)
            cout<<"timestep in execute()="<<dt<<" stepnum "<<steps<<" totaltime"<<globalTime<<endl;
        
        //
        // Detection of negative J and T and soft exit
        //
        for(int s = 0; s < num_species; s++) {
            for(int i=num_cells; i>=0; i--)  {
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
                        
                        //cout<<" NEGATIVE TEMPERATURE at s/i = "<<s<<"/"<<i<<" in timestep="<<dt<<" stepnum "<<steps<<" totaltime "<<globalTime<<endl;
                        globalTime = 1.1*t_max; //This secures that the program exits smoothly after the crash
                    }
            }
            
            if(crashed_T > 0) {
                cout<<endl<<">>> CRASH <<< DUE TO NEGATIVE TEMPERATURES, crash_imin/imax = "<<crash_T_imin<<"/"<<crash_T_imax<<" sample T = "<<crashed_temperature<<" num of crashed cells/total cells = "<<crash_T_numcells<<"/"<<num_cells<<"  crashed species name = "<<species[crashed_T-1].speciesname<<endl; 
                cout<<" @dt="<<dt<<" stepnum "<<steps<<" Crashtime "<<crashtime<<endl;
                cout<<"Writing crash dump into last output and exiting program."<<endl;
            } 
                    
        }
        for(int b = 0; b < num_bands_out; b++) {
            for(int i=num_cells+2-num_ghosts; i>=num_ghosts; i--)  {
                    if(Jrad_FLD(i,b) < 0) {
                        
                        if(Jrad_FLD(i,b) > -1e-10) {
                           Jrad_FLD(i,b) *= -1.; 
                        }
                        else {
                            
                            if(crashtime_already_assigned == 0) {
                                crashtime = globalTime;
                                crashtime_already_assigned = 1;
                            }
                            
                            crashed_meanintensity = Jrad_FLD(i,b);
                            crashed_J = b+1;
                            crash_J_imin = (i<crash_J_imin)? i : crash_J_imin;
                            crash_J_imax = (i>crash_J_imax)? i : crash_J_imax;
                            crash_J_numcells++;
                            
                            //cout<<" NEGATIVE RAD DENSITY at b/i = "<<b<<"/"<<i<<"  at timestep "<<dt<<" stepnum "<<steps<<" totaltime "<<globalTime<<endl;
                            //cout<<" Writing crash dump into last output and exiting program."<<endl;
                            globalTime = 1.1*t_max;    
                            
                        }
                    }
            }
            
            if(crashed_J > 0) {
                cout<<endl<<">>> CRASH <<< DUE TO NEGATIVE J, crash_imin/imax = "<<crash_J_imin<<"/"<<crash_J_imax<<" sample J = "<<crashed_meanintensity<<" num of crashed cells/total cells = "<<crash_J_numcells<<"/"<<num_cells<<"  crashed band number = "<<crashed_J-1<<endl; 
                cout<<" @dt="<<dt<<" stepnum "<<steps<<" Crashtime "<<crashtime<<endl;
                cout<<"Writing crash dump into last output and exiting program."<<endl;
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
            
    //checksum /= num_species;
    //checksum /= (num_cells+2);
    
    cout<<"Finished at time="<<globalTime<<" after steps="<<steps<<" with num_cells="<<num_cells<<" and checksum = "<<checksum<<endl;
    
    for(int s=0; s<num_species; s++) 
        species[s].print_AOS_component_tofile(-1);
    print_diagnostic_file(-1);
    print_monitor(-1);
    user_output_function(-1) ;

        
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


//
// Species::execute: This species very own riemann solution, before its modified globally by applying frictional forces
//

void c_Species::execute(std::vector<AOS>& u_in, std::vector<AOS>& dudt) {
    
        
        //
        // Step 1: Boundary values
        //
        //
  
        apply_boundary_left(u_in) ;
        apply_boundary_right(u_in) ;
        
        if(USE_WAVE==1) 
            add_wave(u_in, base->globalTime);
        
        //if(steps==0)
        //    cout<<" initial ghost momentum after start: "<<left_ghost.u2<<endl;
        
        if(debug >= 2)
            cout<<"Done."<<endl<<" Before compute pressure in species "<<speciesname<<"... ";
        
        compute_pressure(u_in);
        
        
        if(debug >= 2)
            cout<<"Done. Starting edge-states."<<endl;
        
        reconstruct_edge_states() ;
        
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
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
            
            //if( (debug > 1) && ( j==1 || j==num_cells || j==(num_cells/2) )) {
            if( (debug > 3) && ( j==5 || j==1 || j==2 || j==3 || j==4 || j==6 || j==7 || j==8 || j==9 ) && (base->steps==1 || base->steps==0)) {
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
                cin>>alpha;
            }
        }
        
        
        if(debug >= 2)
            cout<<" Before updating rad fluxes... ";
    
} 

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
    //double pS     = (pr*dl*(SL-ul) - pl*dr*(SR-ur) - dl*dr*(SL-ul)*(SR-ur)*(ur-ul)) / (dl*(SL - ul)-dr*(SR-ur)) ;
    
    if ((SL <= 0) &&  (SS >= 0)) {
        AOS FL         = AOS (mom_l, mom_l * ul + pl, ul * (El + pl) );
        double comp3_L = El/dl + (SS-ul)*(SS + pl/(dl*(SL-ul)));
        AOS US_L       = AOS(1,SS, comp3_L) * dl * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - AOS(dl, mom_l, El))  * SL;    
           
        flux = FS_L;
        option = 1;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (mom_r, mom_r * ur + pr, ur * (Er + pr) );
        double comp3_R = Er/dr + (SS-ur)*(SS + pr/(dr*(SR-ur)));
        AOS US_R       = AOS(1,SS, comp3_R) * dr * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - AOS(dr, mom_r, Er)) * SR;
        
        flux = FS_R;
        option = 2;
    }
    else if (SL >= 0) {
        AOS FL = AOS (mom_l, mom_l * ul + pl, ul * (El + pl) );
        flux = FL;
        option = 3;
    }
    else if (SR <= 0) {
        AOS FR = AOS(mom_r, mom_r * ur + pr, ur * (Er + pr) );
        flux = FR;
        option= 4 ;
    }
    if((debug > 2) && (j>0 && j<13))
        cout<<" opt = "<<option<<" SL/SS/SR = "<<SL<<"/"<<SS<<"/"<<SR<<" f = "<<flux.u1<<"/"<<flux.u2<<"/"<<flux.u3<<" ((f.u2-Pl)/Pl)-1 = "<<((flux.u2-pl)/pl)<<" ((f.u2-Pr)/Pr)-1 = "<<((flux.u2-pr)/pr)<<endl;
    return flux;
}


// Riemann problem for (almost) pressure-less dust following Leveque (2004), Pelanti & Leveque (2006)
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

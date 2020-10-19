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
    const int maxsteps = 1e9;
    //double pressure_temp;
        
    cout<<endl<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfacotr="<<cflfactor<<" and num_species = "<<num_species<<endl;
    if(num_species == 0) 
        cout<<"WARNING: No species specified! I cannot work like that."<<endl;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure();
        
        dt = get_cfl_timestep();
        
        //
        // Step 0: Update gravity. Important for self-gravity and time-dependent smoothing length
        //
        //init_grav_pot();
        if(debug >= 2)
            cout<<"Beginning timestep "<<steps<<endl;
        
        if(use_self_gravity==1)
            update_mass_and_pot();
        
        if(debug >= 2)
            cout<<"Before fluxes... ";
        
        //
        // Output data, when required. Keep the output here, between the flux and conserved variable update, so that the 
        // zeroth output has the initialized conserved values, but already the first fluxes.
        //
        
        //TODO: DO for all species...
         if(steps==0) {
             for(int s=0; s<num_species; s++)
                species[s].print_AOS_component_tofile((int) output_counter);
             
             output_counter+=1.;
         }
         if(globalTime > output_counter*output_time){
             if(debug >= 1)
                 cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
             
             for(int s=0; s<num_species; s++)
                species[s].print_AOS_component_tofile((int)output_counter); 
             
             output_counter+=1.; 
         }
         
        //
        // Do all explicit and implicit operations for one timestep on the entire grid for all species
        //
        for(int s = 0; s < num_species; s++)
            species[s].execute();
        
        if(num_species > 2)
            compute_friction_step(); 
        
        globalTime += dt;
        steps++;
        
        if(steps==1)
            cout<<"Initial sound crossing time = "<<max_snd_crs_time<<endl;
        
        if(debug > 1)
            cout<<"timestep in execute()="<<dt<<" stepnum "<<steps<<" totaltime"<<globalTime<<endl;
        
    
    }
    cout<<endl;
    
    //Print successful end result
    cout<<"Finished at time="<<globalTime<<" after steps="<<steps<<" with num_cells="<<num_cells<<endl;
    
    for(int s=0; s<num_species; s++)
        species[s].print_AOS_component_tofile(-666);
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

void c_Species::execute() {
    
        if(use_rad_fluxes==1)
            update_radiation();
        
        //
        // Step 1: Boundary values
        //
        //
  
        apply_boundary_left() ;
        apply_boundary_right() ;
        
        if(USE_WAVE==1) 
                add_wave(base->globalTime, WAVE_PERIOD, WAVE_AMPLITUDE);
        
        //if(steps==0)
        //    cout<<" initial ghost momentum after start: "<<left_ghost.u2<<endl;
        
        if(debug >= 2)
            cout<<"Done."<<endl<<" Before compute pressure... ";
        
        compute_pressure();
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
        //
        // Step 2: Compute fluxes and sources
        //
        for(int j=1; j <= num_cells; j++) {
            
            if(j==1) 
                flux[0] = hllc_flux(u[j-1], u[j], j-1, j); //Flux with ghost cell on left
            
            flux[j] = hllc_flux(u[j], u[j+1], j, j+1); //Everything else is taken into account by the right flux

            /*if(j==2) {
                char alp;ms are which form of objects which have a velocity standard deviation comparable to their m
                cout<<"flux[0]="<<flux[0].u1<<" / "<<flux[0].u2<<" / "<<flux[0].u3<<endl;
                cout<<"flux[1]="<<flux[1].u1<<" / "<<flux[1].u2<<" / "<<flux[1].u3<<endl;
                //cout<<"flux[2]="<<flux[2].u1<<" / "<<flux[2].u2<<" / "<<flux[2].u3<<endl;
                
                cout<<"u[0]="<<u[0].u1<<" / "<<u[0].u2<<" / "<<u[0].u3<<endl;
                cout<<"u[1]="<<u[1].u1<<" / "<<u[1].u2<<" / "<<u[1].u3<<endl;
                cout<<"u[2]="<<u[2].u1<<" / "<<u[2].u2<<" / "<<u[2].u3<<endl;
                cout<<"phi[0,1,2]="<<phi[0]<<" / "<<phi[1]<<" / "<<phi[2]<<endl;
                cin>>alp;
            }*/
        }
        
        if(debug >= 2)
            cout<<"Done. Starting sources."<<endl;
        
        for(int j=1; j<=num_cells; j++) {
            source[j]          = source_grav(u[j], j);
            source_pressure[j] = AOS(0, -(base->source_pressure_prefactor_left[j] * pressure_l[j] - base->source_pressure_prefactor_right[j] * pressure_r[j])  ,0); 
        }
        
        
        //
        // Step 3: Add it all up to update the conserved variables
        //
        
        //#pragma omp simd
        for(int j=1; j<=num_cells; j++) {
            u[j] = u[j] + (flux[j-1] * base->surf[j-1] - flux[j] * base->surf[j]) * base->dt / base->vol[j] + (source[j] + source_pressure[j]) * base->dt;
            
            if( (debug > 0) && ( j==1 || j==num_cells || j==(num_cells/2) )) {
                char alpha;
                cout<<"Debuggin fluxes in cell i= "<<j<<endl; 
                cout<<"fl.u1 = "<<flux[j-1].u1<<": fr.u1 = "<<flux[j].u1<<endl;
                cout<<"fl.u2 = "<<flux[j-1].u2<<": fr.u2 = "<<flux[j].u2<<endl;
                cout<<"fl.u3 = "<<flux[j-1].u3<<": fr.u3 = "<<flux[j].u3<<endl;
                cout<<"Cartesian fluxes: Fl-Fr+s = "<<((flux[j-1].u2 - flux[j].u2)/base->dx[j] + source[j].u2)<<endl;
                cout<<"D_surface/volume="<<(0.5*(base->surf[j]-base->surf[j-1])/base->vol[j])<<" vs. 1/dx="<<(1./base->dx[j])<<endl;
                cout<<endl;
                cout<<"s = "<<source[j].u2<<endl;
                cout<<"sP = "<<source_pressure[j].u2<<endl;
                
                cout<<"Al*Fl - Ar*Fr = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) /base->vol[j])<<endl;
                cout<<"Al*Fl - Ar*Fr + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + source_pressure[j].u2)<<endl;
                cout<<"dP/dr = "<<((pressure_l[j] - pressure_r[j])/base->dx[j])<<endl;
                cout<<"dP/dr + S = "<<((pressure_l[j] - pressure_r[j])/base->dx[j] + source[j].u2)<<endl;
                cout<<endl;
                cout<<"Al*Fl - Ar*Fr + s = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2))<<endl;
                cout<<"Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2 +source_pressure[j].u2))<<endl;
                cin>>alpha;
            }
        }
        
        
        if(debug >= 2)
            cout<<" Before updating rad fluxes... ";
        
        if(use_rad_fluxes) {
            for(int j=1; j<=num_cells; j++) {
                u[j].u3 = u[j].u3 + (radiative_flux[j-1] - radiative_flux[j]); //Dummy physics for now
            }
        }
        
    
    
} 

AOS c_Species::hllc_flux(AOS &leftval, AOS &rightval, const int &jleft, const int& jright) 
{
    AOS flux;
    
    //Speed of gas
    double ul = speed[jleft]; 
    double ur = speed[jright];
    
    double pl = pressure_r[jleft];  
    double pr = pressure_l[jright]; 
    
    //Speed of shocks
    double SL = ul - std::sqrt(gamma_adiabat[jleft]*pl/leftval.u1);
    double SR = ur + std::sqrt(gamma_adiabat[jright]*pr/rightval.u1);
    
    //Intermediate values in the star region, equations 10.30 -10.39 in Toro
    double SS     = ( pr-pl+leftval.u2*(SL - ul)-rightval.u2*(SR-ur) )/(leftval.u1*(SL - ul)-rightval.u1*(SR-ur) );
    
    if ((SL <= 0) &&  (SS >= 0)) {
        
        AOS FL         = AOS (leftval.u2, leftval.u2 * ul  + pl, ul * (leftval.u3 + pl) );
        double comp3_L = leftval.u3/leftval.u1 + (SS-ul)*(SS + pl/(leftval.u1*(SL-ul)));
        AOS US_L       = AOS(1,SS, comp3_L) * leftval.u1  * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - leftval)  * SL;    
           
        flux = FS_L;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (rightval.u2,rightval.u2 * ur + pr, ur * (rightval.u3+ pr) );
        double comp3_R = rightval.u3/rightval.u1 + (SS-ur)*(SS + pr/(rightval.u1*(SR-ur)));
        AOS US_R       = AOS(1,SS, comp3_R) * rightval.u1 * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - rightval) * SR;
        
        flux = FS_R;
    }
    else if (SL >= 0) {
        AOS FL = AOS (leftval.u2, leftval.u2 * ul  + pl, ul * (leftval.u3 + pl) );
        flux = FL;
    }
    else if (SR <= 0) {
        AOS FR = AOS (rightval.u2,rightval.u2 * ur + pr, ur * (rightval.u3 + pr) );
        flux = FR;
    }
    return flux;
}




AOS c_Species::analytic_flux(AOS &input_vec, const int &j) {
     
    return AOS (input_vec.u2, 
                input_vec.u2*input_vec.u2/input_vec.u1 + pressure[j], 
                input_vec.u2/input_vec.u1 * (input_vec.u3 + pressure[j]) );
}

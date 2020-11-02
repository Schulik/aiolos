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
    const int maxsteps = 1e9;
    //double pressure_temp;
        
    cout<<endl<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfactor="<<cflfactor<<" and num_species = "<<num_species<<endl;
    if(num_species == 0) 
        throw std::invalid_argument("WARNING: No species specified! I cannot work like that.") ;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {
        
        if(steps==0)
            for(int s = 0; s < num_species; s++)
                species[s].compute_pressure(species[s].u);
        
        dt = get_cfl_timestep();
        dt = std::min(dt, t_max - globalTime) ;
        
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
        
         if(steps==0 || steps==1) {
             for(int s=0; s<num_species; s++) {
                species[s].print_AOS_component_tofile((int) output_counter);
            }
            print_monitor((int)monitor_counter);
                
             monitor_counter+=1.;
             output_counter +=1.;
         }
         if(globalTime > output_counter*output_time){
             if(debug >= 1)
                 cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
             
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
        for(int s = 0; s < num_species; s++)
            species[s].execute(species[s].u, species[s].dudt[0]);
        
        // Just do the basic update here
        for(int s=0; s < num_species; s++) {
            for(int j=0; j < num_cells+1; j++)
                species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
        }
        
        

        if (order == IntegrationType::first_order) {
            globalTime += dt ;
        }
        else if (order == IntegrationType::second_order) {
            // 2nd step evaluated at t+dt
            globalTime += dt ;

            if(use_self_gravity==1)
                update_mass_and_pot();

            for(int s = 0; s < num_species; s++)
                species[s].execute(species[s].u, species[s].dudt[1]);

            for(int s = 0; s < num_species; s++)
                for(int j = 0; j < num_cells+2; j++) {
                    species[s].u[j] += (species[s].dudt[1][j] -  species[s].dudt[0][j])*dt / 2 ; 
                    
                }
        }
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure(species[s].u);
        
        if(alpha_collision > 0 && num_species > 1) {
            if(friction_solver == 0)
                compute_friction_analytical(); 
            else if(friction_solver == 1) 
                compute_friction_numerical_static();
            else
                compute_friction_numerical_dynamic();
                //compute_friction_numerical_sparse();
        }
        
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
        species[s].print_AOS_component_tofile(-1);
    print_monitor(-1);
        
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
    
        if(base->use_rad_fluxes==1)
            update_radiation(u_in);
        
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
            cout<<"Done."<<endl<<" Before compute pressure... ";
        
        compute_pressure(u_in);
        reconstruct_edge_states() ;
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
        //
        // Step 2: Compute fluxes and sources
        //
        for(int j=0; j <= num_cells; j++) {
            
            flux[j] = hllc_flux(j); 

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
            source[j]          = source_grav(u_in[j], j);
            source_pressure[j] = AOS(0, -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
                                          base->source_pressure_prefactor_right[j] * prim_r[j].pres)  ,0); 
        }
        
        
        //
        // Step 3: Add it all up to update the conserved variables
        //
        
        //#pragma omp simd
        for(int j=1; j<=num_cells; j++) {
            dudt[j] = (flux[j-1] * base->surf[j-1] - flux[j] * base->surf[j]) / base->vol[j] + (source[j] + source_pressure[j]) ;
            
            if( (debug > 1) && ( j==1 || j==num_cells || j==(num_cells/2) )) {
                //char alpha;
                cout<<"Debuggin fluxes in cell i= "<<j<<endl; 
                cout<<"     fl.u1 = "<<flux[j-1].u1<<": fr.u1 = "<<flux[j].u1<<endl;
                cout<<"     fl.u2 = "<<flux[j-1].u2<<": fr.u2 = "<<flux[j].u2<<endl;
                cout<<"     fl.u3 = "<<flux[j-1].u3<<": fr.u3 = "<<flux[j].u3<<endl;
                cout<<"     Cartesian fluxes: Fl-Fr+s = "<<((flux[j-1].u2 - flux[j].u2)/base->dx[j] + source[j].u2)<<endl;
                cout<<"     D_surface/volume="<<(0.5*(base->surf[j]-base->surf[j-1])/base->vol[j])<<" vs. 1/dx="<<(1./base->dx[j])<<endl;
                cout<<endl;
                cout<<"     s = "<<source[j].u2<<endl;
                cout<<"     sP = "<<source_pressure[j].u2<<endl;
                
                cout<<"     Al*Fl - Ar*Fr = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) /base->vol[j])<<endl;
                cout<<"     Al*Fl - Ar*Fr + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + source_pressure[j].u2)<<endl;
                cout<<"     dP/dr = "<<((prim_l[j].pres - prim_r[j].pres)/base->dx[j])<<endl;
                cout<<"     dP/dr + S = "<<((prim_l[j].pres - prim_r[j].pres)/base->dx[j] + source[j].u2)<<endl;
                cout<<endl;
                cout<<"     Al*Fl - Ar*Fr + s = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2))<<endl;
                cout<<"     Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u2 * base->surf[j-1] - flux[j].u2 * base->surf[j]) / base->vol[j] + (source[j].u2 +source_pressure[j].u2))<<endl;
                //cin>>alpha;
            }
        }
        
        
        if(debug >= 2)
            cout<<" Before updating rad fluxes... ";
        
        if(base->use_rad_fluxes) {
            for(int j=1; j<=num_cells; j++) {
                dudt[j].u3 = (radiative_flux[j-1] - radiative_flux[j]); //Dummy physics for now
            }
        }
        
    
    
} 

AOS c_Species::hllc_flux(int j) 
{
    int jleft = j, jright = j+1;
    AOS flux;
    
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

    //Speed of shocks
    double SL = ul - prim_r[jleft].sound_speed ;
    double SR = ur + prim_l[jright].sound_speed ;
    
    //Intermediate values in the star region, equations 10.30 -10.39 in Toro
    double SS     = ( pr-pl+mom_l*(SL - ul)-mom_r*(SR-ur) )/(dl*(SL - ul)-dr*(SR-ur) );
    
    if ((SL <= 0) &&  (SS >= 0)) {
        AOS FL         = AOS (mom_l, mom_l * ul + pl, ul * (El + pl) );
        double comp3_L = El/dl + (SS-ul)*(SS + pl/(dl*(SL-ul)));
        AOS US_L       = AOS(1,SS, comp3_L) * dl * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - AOS(dl, mom_l, El))  * SL;    
           
        flux = FS_L;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (mom_r, mom_r * ur + pr, ur * (Er + pr) );
        double comp3_R = Er/dr + (SS-ur)*(SS + pr/(dr*(SR-ur)));
        AOS US_R       = AOS(1,SS, comp3_R) * dr * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - AOS(dr, mom_r, Er)) * SR;
        
        flux = FS_R;
    }
    else if (SL >= 0) {
        AOS FL = AOS (mom_l, mom_l * ul + pl, ul * (El + pl) );
        flux = FL;
    }
    else if (SR <= 0) {
        AOS FR = AOS(mom_r, mom_r * ur + pr, ur * (Er + pr) );
        flux = FR;
    }
    return flux;
}


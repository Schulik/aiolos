#include "main.h"

//extern AOS* init_AOS(int num);

void hydro_run::execute() { 
     
    int steps = 0; //, i,j;
    double output_counter = 0;
    const int maxsteps = 1e9;
    double dt = 1e-3; //get_cfl_timestep();
    
    cout<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfacotr="<<cflfactor<<endl;
    
    compute_pressure();
    //print_AOS_component_tofile(x_i12, u, flux, 0);
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop, let'se goooooooooooooooo! *Insert Super Mario*    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {
        
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
        
        if(use_rad_fluxes==1)
            update_radiation();
        
        //
        // Step 1: Boundary values
        //
        // For noflux-boundaries, boundaries 1, we symmetrize gravity in order to assure zero flux with the source function
        //
        ghost_xi_left   = x_i12[0] - (x_i12[1] - x_i12[0]);
        ghost_xi_right  = x_i12[num_cells-1] + (x_i12[num_cells-1] - x_i12[num_cells-2]);
        if(boundaries_number == 4) {
            phi_left_ghost  = get_phi_grav(x_i12[0],           planet_mass);
            phi_right_ghost = get_phi_grav(x_i12[num_cells-1], planet_mass);
        } else {
            phi_left_ghost  = get_phi_grav(ghost_xi_left, planet_mass);
            phi_right_ghost = get_phi_grav(ghost_xi_right, enclosed_mass[num_cells-1]);
        }
        
        
        //Compute the boundary values, WARNING: ALWAYS DO THAT FIRST
        if(boundaries_number == 1) {
            if(steps==0)
                cout<<"Const boundaries on both sides, given by inital conditions"<<endl;
            boundaries_const_both(left_ghost, u[0], u[num_cells-1], right_ghost );
        }
        else if(boundaries_number == 2) {
            if(steps==0)
                cout<<"Open boundaries on both sides / Wall+open if using gravity."<<endl;
            boundaries_open_both(left_ghost, u[0], u[1], u[num_cells-2], u[num_cells-1], right_ghost );
            
        }
        else if(boundaries_number == 3) {
            if(steps==0)
                cout<<"Wall and inflow boundaries."<<endl;
            boundaries_planet_mdot(left_ghost, u[0], u[num_cells-1], right_ghost );
        }
        else if(boundaries_number == 4) {
            if(steps==0)
                cout<<"Wall boundaries on both sides."<<endl;
            boundaries_wall_both(left_ghost, u[0], u[num_cells-1], right_ghost );
        }
        else {
            if(steps==0)
                cout<<"WARNING: default boundaries!"<<endl;
            boundaries_const_both(left_ghost, AOS(1,1,1), AOS(1,1,1), right_ghost );
        }
            
        
        //if(steps==0)
        //    cout<<" initial ghost momentum after start: "<<left_ghost.u2<<endl;
        
        //boundval_left  = conserved2[0];          //In general, this should be a function of the conserved[0] value
        //boundval_right = conserved2[num_cells-1];
        if(debug >= 2)
            cout<<"Done."<<endl<<" Before compute pressure... ";
        compute_pressure();
        
        if(debug >= 2)
            cout<<"Done. Starting fluxes."<<endl;
        
        //
        // Step 2: Compute fluxes and sources
        //
        for(int j=0; j<num_cells; j++) {
            
            if(j==0) {
                flux[0] = hllc_flux(left_ghost, u[j],  phi_left_ghost, phi[j],   0, j); //Flux with ghost cell on left
                flux[1] = hllc_flux(u[j],       u[j+1],phi[j],         phi[j+1], j, j+1); //First regular flux, jumpstarting the chain of fluxes
            }
            else if(j==num_cells-1)
                flux[j+1] = hllc_flux(u[j],     right_ghost, phi[j], phi_right_ghost, j, j); //Last flux with ghost cell on right
            else {
                flux[j+1] = hllc_flux(u[j], u[j+1], phi[j], phi[j+1], j, j+1); //Regular flux for most of the domain
            }
            
        }
        
        if(debug >= 2)
            cout<<"Done. Starting sources."<<endl;
        
        for(int j=0; j<num_cells; j++) {
        
            //cout<<" "<<j;
            
            if(j == 0)
                source[j]    = source_grav(u[j], phi_left_ghost, phi[j+1]);
            else if (j == num_cells-1) {
                source[j]    = source_grav(u[j], phi[j-1], phi_right_ghost);
            }
            else
                source[j]    = source_grav(u[j], phi[j-1], phi[j+1]);
        }
        
        //
        // Output data, when required. Keep the output here, between the flux and conserved variable update, so that the 
        // zeroth output has the initialized conserved values, but already the first fluxes.
        //
        if(steps==0) {
            print_AOS_component_tofile(x_i12, u, flux, (int) output_counter);
            output_counter+=1.;
        }
        if(globalTime > output_counter*output_time){
            if(debug >= 1)
                cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
            print_AOS_component_tofile(x_i12, u, flux, (int) output_counter); 
            output_counter+=1.; 
        }
        
        //
        // Step 3: Add it all up to update the conserved variables
        //
        
        //#pragma omp simd
        for(int i=0; i<num_cells; i++) {
            
            u[i] = u[i] + (flux[i] - flux[i+1] + source[i]) * dt/dx[i];
            
            //u[j].u1 += dt/dx[j] * ( flux[j].u1 - flux[j+1].u1);
            //u[j].u2 += dt/dx[j] * ( flux[j].u2 - flux[j+1].u2);
            //u[j].u3 += dt/dx[j] * ( flux[j].u3 - flux[j+1].u3);
        }
        
        if(debug >= 2)
            cout<<" Before updating rad fluxes... ";
        
        if(use_rad_fluxes) {
            for(int i=0; i<num_cells-1; i++) {
                u[i].u3 = u[i].u3 + (radiative_flux[i] - radiative_flux[i+1]) * 0.1; //Just tryout numbers for now
            }
        }
        
        
        globalTime += dt;
        steps++;
        
        //This is already the next timestep. Note that the initial timestep is properly initialized!
        dt = get_cfl_timestep();
        
        if(steps==1)
            cout<<"Initial sound crossing time = "<<snd_crs_time<<endl;
        
        if(debug >= 1)
            cout<<"timestep in execute()="<<dt<<" stepnum "<<steps<<" totaltime"<<globalTime<<endl;
        
        
        
    }
    cout<<endl;
    
    //Print successful end result
    cout<<"Finished at time="<<globalTime<<" after steps="<<steps<<" with num_cells="<<num_cells<<endl;
    print_AOS_component_tofile(x_i12, u, flux, -666);
}


AOS hydro_run::hllc_flux(AOS &leftval, AOS &rightval, double &phi_l, double &phi_r, const int &jleft, const int& jright) 
{
    //int decision = -1;
    AOS flux;
 
    /*double SL = wave_speeds2(leftval,  d_minusone, jleft); 
    double SR = wave_speeds2(rightval, d_plusone,  (jleft+1));
    AOS    FL = analytic_flux2(leftval, jleft);
    AOS    FR = analytic_flux2(rightval,(jleft+1)); */
    
    //Speed of gas
    double ul = speed[jleft]; //leftval.u2  / leftval.u1; 
    double ur = speed[jright];      //rightval.u2 / rightval.u1;
    
    //Pressures, TODO: Replace with grav_pressure(pressure2[jleft])
    //double pl = pressure[jleft]; //(gamma_adiabat-1.)*(leftval.u3 - 0.5* leftval.u2*leftval.u2 / leftval.u1 ); //pressure[jleft]; 
    //double pr = pressure[jright]; //(gamma_adiabat-1.)*(rightval.u3 - 0.5* rightval.u2*rightval.u2 / rightval.u1 ); //pressure[jright]; 
    double pl = get_p_hydrostatic(leftval, phi_r, phi_l, jleft);
    double pr = get_p_hydrostatic(rightval, phi_l, phi_r, jright);
    
    //Speed of shocks
    double SL = ul - std::sqrt(gamma_adiabat*pl/leftval.u1);
    double SR = ur + std::sqrt(gamma_adiabat*pr/rightval.u1);
    
    //double SL = ul - std::sqrt(gamma_adiabat*pressure2[jleft]/leftval.u1);
    //double SR = ur + std::sqrt(gamma_adiabat*pressure2[jright]/rightval.u1);
    
    //Analytic fluxes
    //AOS FL = AOS (leftval.u2, leftval.u2 * ul  + pl, ul * (leftval.u3 + pl) );
    //AOS FR = AOS (rightval.u2,rightval.u2 * ur + pr, ur * (rightval.u3+ pr) );
     
    //Intermediate values in the star region, equations 10.30 -10.39 in Toro
    double SS     = ( pr-pl+leftval.u2*(SL - ul)-rightval.u2*(SR-ur) )/(leftval.u1*(SL - ul)-rightval.u1*(SR-ur) );
    
    if ((SL <= 0) &&  (SS >= 0)) {
        
        AOS FL         = AOS (leftval.u2, leftval.u2 * ul  + pl, ul * (leftval.u3 + pl) );
        double comp3_L = leftval.u3/leftval.u1 + (SS-ul)*(SS + pl/(leftval.u1*(SL-ul)));
        AOS US_L       = AOS(1,SS, comp3_L) * leftval.u1  * (SL - ul)/(SL-SS);
        AOS FS_L       = FL + (US_L - leftval)  * SL;    
           
        flux = FS_L;
        //decision = 0;
    }
    else if ((SS <= 0) && (SR >= 0)) {
        
        AOS FR         = AOS (rightval.u2,rightval.u2 * ur + pr, ur * (rightval.u3+ pr) );
        double comp3_R = rightval.u3/rightval.u1 + (SS-ur)*(SS + pr/(rightval.u1*(SR-ur)));
        AOS US_R       = AOS(1,SS, comp3_R) * rightval.u1 * (SR - ur)/(SR-SS);
        AOS FS_R       = FR + (US_R - rightval) * SR;
        
        flux = FS_R;
        //decision = 1;
    }
    else if (SL >= 0) {
        AOS FL = AOS (leftval.u2, leftval.u2 * ul  + pl, ul * (leftval.u3 + pl) );
        flux = FL;
        //decision = 2;
    }
    else if (SR <= 0) {
        AOS FR = AOS (rightval.u2,rightval.u2 * ur + pr, ur * (rightval.u3 + pr) );
        flux = FR;
        //decision = 3;
    }
    return flux;
}

#include "main.h"

void hydro_run::execute() { 
     
    steps = 0;
    double output_counter = 0;
    const int maxsteps = 1e9;
    dt = 1e-3; //get_cfl_timestep();
    
    cout<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfacotr="<<cflfactor<<endl;
    
    compute_pressure();
    //print_AOS_component_tofile(x_i12, u, flux, 0);
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
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
        
        //ghost_xi_left   = x_i12[0] - (x_i12[1] - x_i12[0]);
        //ghost_xi_right  = x_i12[num_cells-1] + (x_i12[num_cells-1] - x_i12[num_cells-2]);
        if(boundaries_number == 4) {
            phi[0]           = get_phi_grav(x_i12[1],         planet_mass);
            phi[num_cells+1] = get_phi_grav(x_i12[num_cells], planet_mass);
        } else {
            phi[0]           = get_phi_grav(x_i12[1],  planet_mass);
            phi[num_cells+1] = get_phi_grav(x_i12[num_cells+1], planet_mass);
        }
        
        
        //Compute the boundary values
        if(boundaries_number == 1) {
            if(steps==0)
                cout<<"Const boundaries on both sides, given by inital conditions"<<endl;
            boundaries_const_both(u[0], u[1], u[num_cells], u[num_cells+1] );
        }
        else if(boundaries_number == 2) {
            if(steps==0)
                cout<<"Open boundaries on both sides / Wall+open if using gravity."<<endl;
            boundaries_open_both(u[0], u[1], u[2], u[num_cells-1], u[num_cells], u[num_cells+1] );
            
        }
        else if(boundaries_number == 3) {
            if(steps==0)
                cout<<"Wall and inflow boundaries."<<endl;
            boundaries_planet_mdot(u[0], u[1], u[num_cells], u[num_cells+1] );
        }
        else if(boundaries_number == 4) {
            if(steps==0)
                cout<<"Wall boundaries on both sides."<<endl;
            boundaries_wall_both(u[0], u[1], u[num_cells], u[num_cells+1] );
        }
        else {
            if(steps==0)
                cout<<"WARNING: default boundaries!"<<endl;
            boundaries_const_both(u[0], AOS(1,1,1), AOS(1,1,1), u[num_cells+1] );
        }
            
        
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
                flux[0] = hllc_flux(u[j-1], u[j],  phi[j-1], phi[j],   j-1, j); //Flux with ghost cell on left
            
            flux[j] = hllc_flux(u[j], u[j+1], phi[j], phi[j+1], j, j+1); //Everything else is taken into account by the right flux

            /*if(j==2) {
                char alp;
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
            source_pressure[j] = AOS(0, pressure[j] * (surf[j]-surf[j-1])/vol[j] ,0); //Metric term 2P/r in a discretization that well-balances constant states
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
        for(int j=1; j<=num_cells; j++) {
            
            //Cartesian
            //u[j] = u[j] + (flux[j-1] - flux[j]) * dt/dx[j] + source[j] * dt;
            
            //Spherical
            u[j] = u[j] + (flux[j-1] * surf[j-1] - flux[j] * surf[j]) * dt/vol[j] + (source[j] + source_pressure[j]) * dt;
            
            /*if(i==1 || i==num_cells || i==20) {
                char alpha;
                cout<<"Debuggin fluxes in cell i= "<<i<<": fl-fr = "<<((flux[i-1].u2 - flux[i].u2)/dx[i])<<endl;
                cout<<"Debuggin fluxes: s = "<<source[i].u2<<endl;
                cout<<"Debuggin fluxes: fl-fr+s = "<<((flux[i-1].u2 - flux[i].u2)/dx[i] + source[i].u2)<<endl;
                cin>>alpha;
            }*/
        }
            
        

        
        
        if(debug >= 2)
            cout<<" Before updating rad fluxes... ";
        
        if(use_rad_fluxes) {
            for(int j=1; j<=num_cells; j++) {
                u[j].u3 = u[j].u3 + (radiative_flux[j-1] - radiative_flux[j]); //Dummy physics for now
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
    
    //Pressures with hydrostatic reconstruction
    //double pl = get_p_hydrostatic(leftval, phi_r, phi_l, jleft);
    //double pr = get_p_hydrostatic(rightval, phi_l, phi_r, jright);
    
    double pl = get_p_hydrostatic_nonuniform(jleft,  +1);
    double pr = get_p_hydrostatic_nonuniform(jright, -1);
    
    //Speed of shocks
    double SL = ul - std::sqrt(gamma_adiabat*pl/leftval.u1);
    double SR = ur + std::sqrt(gamma_adiabat*pr/rightval.u1);
    
     
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

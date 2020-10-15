 #include "aiolos.h"

void hydro_run::execute() { 
     
    steps = 0;
    int output_counter = 0;
    const int maxsteps = 1e9;
    dt = get_cfl_timestep();
        
    cout<<"Beginning main loop with num_cells="<<num_cells<<" and timestep="<<dt<<" cflfacotr="<<cflfactor<<endl;
    
    compute_pressure(u);
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    //                                                                         //
    // Simulation main loop                                                    //
    //                                                                         //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~////
    for (globalTime = 0; (globalTime < t_max) && (steps < maxsteps); ) {
        

        if(debug >= 2)
            cout<<"Beginning timestep "<<steps<<endl;
        

        //
        // Output data, when required. Keep the output here, between the flux and conserved variable update, so that the 
        // zeroth output has the initialized conserved values, but already the first fluxes.
        //
        if(globalTime >= output_counter*output_time){
            if(debug >= 1)
                cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
            print_AOS_component_tofile(x_i12, u, flux, output_counter); 
            output_counter++; 
        }
        

        // Do the step
        runge_kutta_stage(u, u, dt) ;

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

    if(globalTime >= output_counter*output_time){
        if(debug >= 1)
            cout<<" Globaltime is "<<globalTime<<" and comparevalue is "<<output_counter<<" "<<output_time<<endl;
        print_AOS_component_tofile(x_i12, u, flux, output_counter); 
        output_counter++; 
    }

    //Print successful end result
    cout<<"Finished at time="<<globalTime<<" after steps="<<steps<<" with num_cells="<<num_cells<<endl;
    print_AOS_component_tofile(x_i12, u, flux, -666);
}

void hydro_run::runge_kutta_stage(std::vector<AOS>& u_in, std::vector<AOS>& u_out, double dt_stage){
    
    //
    // Step 0: Potential
    //
    //
    if(use_self_gravity==1)
        update_mass_and_pot(u_in);

    //
    // Step 1: Boundary values
    //
    //
    
    // Handled in apply boundary_right
    /*
    if(boundaries_number == 4) {
        phi[num_cells+1] = get_phi_grav(x_i12[num_cells], planet_mass);
    } else {
        phi[num_cells+1] = get_phi_grav(x_i12[num_cells+1], planet_mass);
    }
    */
    
    apply_boundary_left(u_in) ;
    apply_boundary_right(u_in) ;
    
    if(USE_WAVE==1) 
        add_wave(globalTime, u_in);
    
    //if(steps==0)
    //    cout<<" initial ghost momentum after start: "<<left_ghost.u2<<endl;
    
    if(debug >= 2)
        cout<<"Done."<<endl<<" Before compute pressure... ";
    compute_pressure(u_in);
    
    if(debug >= 2)
        cout<<"Done. Starting fluxes."<<endl;
    
    //
    // Step 2: Compute fluxes and sources
    //
    for(int j=1; j <= num_cells; j++) {
        
        if(j==1) 
            flux[0] = hllc_flux(u_in[j-1], u_in[j], j-1, j); //Flux with ghost cell on left
        
        flux[j] = hllc_flux(u_in[j], u_in[j+1], j, j+1); //Everything else is taken into account by the right flux
    }
    
    if(debug >= 2)
        cout<<"Done. Starting sources."<<endl;
    
    for(int j=1; j<=num_cells; j++) {
        source[j] = source_grav(u_in[j], j);
        double pressure_temp = 
            - (source_pressure_prefactor_left[j] * pressure_l[j] - source_pressure_prefactor_right[j] * pressure_r[j]);
        source_pressure[j] = AOS(0, pressure_temp  ,0); 
    }
        
    //
    // Step 3: Add it all up to update the conserved variables
    //
    
    //#pragma omp simd
    for(int j=1; j<=num_cells; j++) {
        u_out[j] = u_in[j] 
            + (flux[j-1] * surf[j-1] - flux[j] * surf[j]) * dt_stage/vol[j] 
            + (source[j] + source_pressure[j]) * dt_stage;
        
        
        if( (debug > 0) && ( j==1 || j==num_cells || j==(num_cells/2) )) {
            char alpha;
            cout<<"Debuggin fluxes in cell i= "<<j<<endl; 
            cout<<"fl.u1 = "<<flux[j-1].u1<<": fr.u1 = "<<flux[j].u1<<endl;
            cout<<"fl.u2 = "<<flux[j-1].u2<<": fr.u2 = "<<flux[j].u2<<endl;
            cout<<"fl.u3 = "<<flux[j-1].u3<<": fr.u3 = "<<flux[j].u3<<endl;
            cout<<"Cartesian fluxes: Fl-Fr+s = "<<((flux[j-1].u2 - flux[j].u2)/dx[j] + source[j].u2)<<endl;
            cout<<"D_surface/volume="<<(0.5*(surf[j]-surf[j-1])/vol[j])<<" vs. 1/dx="<<(1./dx[j])<<endl;
            cout<<endl;
            cout<<"s = "<<source[j].u2<<endl;
            cout<<"sP = "<<source_pressure[j].u2<<endl;
            
            cout<<"Al*Fl - Ar*Fr = "<<((flux[j-1].u2 * surf[j-1] - flux[j].u2 * surf[j]) /vol[j])<<endl;
            cout<<"Al*Fl - Ar*Fr + sP = "<<((flux[j-1].u2 * surf[j-1] - flux[j].u2 * surf[j]) /vol[j] + source_pressure[j].u2)<<endl;
            cout<<"dP/dr = "<<((pressure_l[j] - pressure_r[j])/dx[j])<<endl;
            cout<<"dP/dr + S = "<<((pressure_l[j] - pressure_r[j])/dx[j] + source[j].u2)<<endl;
            cout<<endl;
            cout<<"Al*Fl - Ar*Fr + s = "<<((flux[j-1].u2 * surf[j-1] - flux[j].u2 * surf[j]) /vol[j] + (source[j].u2))<<endl;
            cout<<"Al*Fl - Ar*Fr + s + sP = "<<((flux[j-1].u2 * surf[j-1] - flux[j].u2 * surf[j]) /vol[j] + (source[j].u2 +source_pressure[j].u2))<<endl;
            cin>>alpha;
        }
    }

    if(use_rad_fluxes) {
        for(int j=1; j<=num_cells; j++) {
            u_out[j].u3 = u_out[j].u3 + (radiative_flux[j-1] - radiative_flux[j]); //Dummy physics for now
        }
    }
}


AOS hydro_run::hllc_flux(AOS &leftval, AOS &rightval, const int &jleft, const int& jright) 
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
    
    double pl = pressure_r[jleft];  //get_p_hydrostatic_nonuniform(jleft,  +1);
    double pr = pressure_l[jright]; //get_p_hydrostatic_nonuniform(jright, -1);
    
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

#include "aiolos.h"

//extern AOS* init_AOS(int num);

/*
Initializer for a new simulation. Stores all parameters in the object, given to by the wrapper.
    
        par_control: List of control parameters
    
        par_numerics: List of numerical parameters
    
        par_physics: List of physics parameters
*/
hydro_run::hydro_run(string filename, double debug_) {

        debug = debug_ ;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //  Numerical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double dx0;
        dx0              = read_parameter_from_file<double>(filename,"PARI_DOMAIN_DX", debug).value;
        cells_per_decade = read_parameter_from_file<double>(filename,"PARI_CELLS_PER_DECADE", debug).value;
        type_of_grid     = read_parameter_from_file<int>(filename,"PARI_GRID_TYPE", debug).value;
        domain_min       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MIN", debug).value;
        domain_max       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MAX", debug).value;
        
        cout<<"pos1"<<endl;
        cout<<"log10(domain_min) = "<<log10f(domain_min)<<endl;\
        cout<<"log10(domain_max) = "<<log10f(domain_max)<<endl;
        cout<<"cells per decade = "<<cells_per_decade<<endl;
        
        if(type_of_grid == 0)
            num_cells        = (int)((domain_max - domain_min)/dx0);
        else {
            num_cells        = (int) ( (log10f(domain_max) - log10f(domain_min)) * cells_per_decade );
            dx0              = domain_min;
        }
            
        //num_cells must always be read in pretty early, as all other memory allocations in this function depend on it.
        //Note-to-self: Don't change the order of read_parameter_from_file commands.
        
        dt          = read_parameter_from_file<double>(filename,"PARI_TIME_DT", debug).value;
        cflfactor   = read_parameter_from_file<double>(filename,"PARI_CFLFACTOR", debug).value;
        t_max       = read_parameter_from_file<double>(filename,"PARI_TIME_TMAX", debug).value;
        output_time = read_parameter_from_file<double>(filename,"PARI_TIME_OUTPUT", debug).value; 
        globalTime = 0.0;    
    
        cout<<"pos2, cell number = "<<num_cells<<endl;
    
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //    Control parameters for users
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //In main, execution is about to start.

        simname = filename;
        
        boundary_left  = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_LEFT", debug).value;
        boundary_right = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_RIGHT", debug).value;
        //cout<<"Pos2"<<endl;
        problem_number    = read_parameter_from_file<int>(filename,"PARI_PROBLEM_NUMBER", debug).value;
        //cout<<"Pos3"<<endl;
 
        cout<<"boundaries used: "<<boundary_left<<" / "<<boundary_right<<endl;
        
        use_self_gravity  = read_parameter_from_file<int>(filename,"PARI_SELF_GRAV_SWITCH", debug).value;
        use_linear_gravity= read_parameter_from_file<int>(filename,"PARI_LINEAR_GRAV", debug).value;
        use_rad_fluxes    = read_parameter_from_file<int>(filename,"PARI_USE_RADIATION", debug).value;
        
        if(use_self_gravity)
            cout<<"WARNING: Self-gravity switched ON !!!"<<endl;
        else
            cout<<"Self-gravity switched OFF."<<endl;
        
        //  Boundaries explained:
        //	1 = Inflow left, fixed wall right boundaries
        //	2 = Outflow on both ends
        //	3 = Inflow left, outflow right
        //finalplotnumber = 1
        solver = 0;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Simulation data: Grid and variables
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        cout<<"pos3"<<endl;
        
        x_i        = np_zeros(num_cells+1);		//The cell boundaries
        x_i12      = np_zeros(num_cells+2);		//The cell mid positions
        surf       = np_zeros(num_cells+1);     // intercell surfaces
        vol        = np_zeros(num_cells+2);     // cell volumes
        timesteps  = np_zeros(num_cells+2);		    
        timesteps_cs= np_zeros(num_cells+2);	
        finalstep   = np_zeros(num_cells+2);
        dx         = np_zeros(num_cells+2);
        omegaplus  = np_zeros(num_cells+2);
        omegaminus = np_zeros(num_cells+2);
        source_pressure_prefactor_left    = np_zeros(num_cells+2);
        source_pressure_prefactor_right   = np_zeros(num_cells+2);
        
        cout<<"pos4"<<endl;
        
        //
        // Compute cell wall boundaries
        //         First and last two cells near boundaries are uniform
        //
        x_i[0] = domain_min;
        if(type_of_grid==1) {

            double dlogx = pow(10.,1./cells_per_decade);
            
            x_i[1] = x_i[0] * dlogx;
            x_i[2] = x_i[1] + (x_i[1] - x_i[0]);
            
            for(int i=3; i< num_cells-1; i++) {
                x_i[i]   = x_i[i-1] * dlogx;
            }
            double dxlast = x_i[num_cells-2] - x_i[num_cells-3];
            x_i[num_cells-1] = x_i[num_cells-2] + dxlast;
            x_i[num_cells]   = x_i[num_cells-1] + dxlast;
            
            //Assign the last boundary as domain maximum as long as nonuniform grid is in the test-phase
            domain_max = x_i[num_cells];
            cout<<"We have a DOMAIN MAX = "<<domain_max<<endl;
                
        }
        //Uniform grid
        else {
            num_cells = (int)((domain_max-domain_min)/dx0);
            for(int i=1; i<= num_cells; i++) {
                x_i[i] = x_i[i-1] + dx0;
            }
        }
       
        //Compute inter-sphere surfaces
        // Read geometry from file:
        geometry = read_parameter_from_file<Geometry>(filename, "PARI_GEOMETRY", debug, Geometry::cartesian).value;
        switch (geometry) {
            case Geometry::cartesian:
                cout<<"Initializing cartesian geometry."<<endl;
                for(int i=0; i<num_cells+1; i++) {
                    surf[i] = 1 ;
                }
                break;
            case Geometry::cylindrical:
                cout<<"Initializing cylindrical geometry."<<endl;
                for(int i=0; i<num_cells+1; i++) {
                    surf[i] = 2*M_PI * x_i[i] ;
                }
                break;
            case Geometry::spherical:
                cout<<"Initializing spherical geometry."<<endl;
                for(int i=0; i<num_cells+1; i++) {
                    surf[i] = 4*M_PI * x_i[i]*x_i[i] ;
                }
                break;
        }
        
        //Compute shell volumes
        for(int i=1; i<num_cells+1; i++) {
            switch (geometry) {          
            case Geometry::cartesian:
                vol[i] = x_i[i] - x_i[i-1] ;
                break;
            case Geometry::cylindrical:
                vol[i] = M_PI * (x_i[i]*x_i[i] - x_i[i-1]*x_i[i-1]);
                break;
            case Geometry::spherical:
                vol[i] = (4*M_PI/3) * 
                    (x_i[i]*x_i[i]*x_i[i] - x_i[i-1]*x_i[i-1]*x_i[i-1]);
                break;
            }
        }
        
        //Compute cell mid positions. Ghost cells also have mid positions in order to balance their pressure gradients
        // but need to be calculated after this loop
        for(int i=1; i<num_cells+1; i++) {
            x_i12[i] = (0.5 * (x_i[i] + x_i[i-1]) );
        }
        
        //Differences
        for(int i=1; i<num_cells+1; i++) {
            dx[i]  = x_i[i]-x_i[i-1];
        }
        //Ghost dxs, those have only to be generated such that hydrostatic eq is preserved, they don't have to be physical
        dx[0] = dx[1];
        dx[num_cells+1] = dx[num_cells] + (dx[num_cells] - dx[num_cells-1]);
        
        cout<<" FIRST AND LAST DX = "<<dx[0]<<" / "<<dx[num_cells+1]<<endl;
        
        //Ghost cells
        x_i12[0]           = domain_min - (x_i12[1] - x_i[0]);
        x_i12[num_cells+1] = domain_max + (x_i[num_cells] - x_i12[num_cells]); 
        
        //
        //Computing the metric factors for 2nd order non-uniform differentials
        // 
        // Keep in mind that the 0th and num_cells+1 elements are and must stay undefined!
        //
        for(int i=1; i<num_cells+1; i++) {
            omegaminus[i] = 2. * (dx[i-1] + dx[i]) / (dx[i-1] + 2. * dx[i] + dx[i+1]);
            omegaplus[i]  = 2. * (dx[i+1] + dx[i]) / (dx[i-1] + 2. * dx[i] + dx[i+1]);
        }
        //Ghost metric factors, those have only to be generated such that hydrostatic eq is preserved, they don't have to be physical
        omegaminus[0] = omegaminus[1];
        omegaplus[0]  = omegaplus[1];
        
        omegaminus[num_cells+1] = omegaminus[num_cells] + (omegaminus[num_cells] - omegaminus[num_cells-1]);
        omegaplus[num_cells+1]  = omegaplus[num_cells]  + (omegaplus[num_cells] - omegaplus[num_cells-1]);
        
        for(int i=1; i<num_cells+1; i++) {
            source_pressure_prefactor_left[i]  = (surf[i-1]/vol[i] - 1./dx[i]); 
            source_pressure_prefactor_right[i] = (surf[i]  /vol[i] - 1./dx[i]); 
        }
        
                
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Physical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        const_gamma_adiabat = read_parameter_from_file<double>(filename,"PARI_GAMMA", debug).value; //ratio of specific heats
        const_cv            = read_parameter_from_file<double>(filename,"PARI_CONST_CV", debug).value; //ratio of specific heats
        const_T             = read_parameter_from_file<double>(filename,"PARI_CONST_TEMP",  debug).value;
        T_increment         = read_parameter_from_file<double>(filename,"PARI_TEMP_INCR", debug).value;
        
        
        if(use_rad_fluxes==1)
            const_opacity   = read_parameter_from_file<double>(filename,"PARI_CONST_OPAC", debug).value;
        else
            const_opacity   = 1.;
        
        u               = init_AOS(num_cells+2); //{ np_ones(num_cells*3);   //Conserved hyperbolic variables: density, mass flux, energy density
        phi             = np_zeros(num_cells+2);  //Parabolic Variables: gravitational potential
        enclosed_mass   = np_zeros(num_cells+2);
        source          = init_AOS(num_cells+2);  //Parabolic Variables: gravitational potential
        source_pressure = init_AOS(num_cells+2);  //Parabolic Variables: gravitational potential
        flux            = init_AOS(num_cells+1);
        
        gamma_adiabat   = np_somevalue(num_cells+2, const_gamma_adiabat);
        cv              = np_somevalue(num_cells+2, const_cv);
        opacity        = np_somevalue(num_cells+2, const_opacity);
        opticaldepth   = np_zeros(num_cells+2);
        radiative_flux = np_zeros(num_cells+1);
        temperature    = np_somevalue(num_cells+2, const_T);
        
        //
        // Determine which equation of states we are using
        //
        eos_pressure_type       = EOS_pressure_type::adiabatic;
        eos_internal_energy_type = EOS_internal_energy_type::thermal;
        
        //Read EOS specifications from file
        if(eos_pressure_type == EOS_pressure_type::tabulated) {
            
            //Read in data into tab_p;
            read_tabulated_eos_data_pressure("eos_pressure.input");
            
        }
        if(eos_internal_energy_type == EOS_internal_energy_type::tabulated) {
            
            //Read in data into tab_p;
            read_tabulated_eos_data_eint("eos_internal_energy.input");
        }
        
        
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Last thing to do: Initialize derived quantities
        //
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        pressure        = np_zeros(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        pressure_l      = np_zeros(num_cells+2); 
        pressure_r      = np_zeros(num_cells+2); 
        internal_energy = np_zeros(num_cells+2);
        speed           = np_zeros(num_cells+2);
        cs              = np_zeros(num_cells+2);
        
        timecount = 0;
        
        if(debug > 0) {
                cout<<"DEBUGLEVEL 1: Samples of generated values"<<endl;
                //cout<<"rho[0], v[0], e[0] = "<<u[0]<<" "<<u[num_cells]<<" "<<u[2*num_cells]<<endl;
                cout<<"First cell coordinates: |<--"<<x_i[0]<<" // "<<x_i12[0]<<" // "<<x_i[1]<<"-->|"<<endl;
                cout<<"Last cell coordinates:  |<--"<<x_i[num_cells-1]<<" // "<<x_i12[num_cells-1]<<" // "<<x_i[num_cells]<<"-->|"<<endl;
                cout<<"Value of gamma: "<<gamma_adiabat[num_cells]<<endl;
        } 
        
        //
        // Problem 1: two states, left and right, separated at a mid-position
        //
        if(problem_number == 1) {
            
            double u1l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1L", debug).value;
            double u2l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2L", debug).value;
            double u3l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3L", debug).value;
            
            double u1r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1R", debug).value;
            double u2r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2R", debug).value;
            double u3r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3R", debug).value;
            
            SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug).value;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            SHOCK_TUBE_UL = AOS(u1l, u1l*u2l, 0.5*u1l*u2l*u2l + u3l/(gamma_adiabat[0]-1.) );
            SHOCK_TUBE_UR = AOS(u1r, u1r*u2r, 0.5*u1r*u2r*u2r + u3r/(gamma_adiabat[num_cells]-1.));
            
            initialize_shock_tube_test(SHOCK_TUBE_UL, SHOCK_TUBE_UR);
            
            if(debug > 0) 
                cout<<"Successfully initialized problem 1."<<endl;
        }
        //
        // Problem 2: A constant background state, with a planet embedded into it, free to evolve as physics dictates
        //
        else if(problem_number == 2) {
            
            cout<<"Problem 2 pos0."<<endl;
            
            double u1 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1", debug).value;
            double u2 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2", debug).value;
            double u3 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3", debug).value;
            
            cout<<"Problem 2 pos1."<<endl;
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug).value; //in Earth masses
            planet_position = read_parameter_from_file<double>(filename,"PARI_PLANET_POS", debug).value;  //inside the simulation domain
            rs              = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_LENGTH", debug).value; //Gravitational smoothing length in hill 
            rs_time         = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_TIME", debug).value; //Time until we reach rs starting at rs_at_moment
            rs_at_moment    = 0.2;
            init_static_atmosphere = read_parameter_from_file<int>(filename,"PARI_INIT_STATIC", debug).value; //Yesno, isothermal at the moment
            
            cout<<"Problem 2 pos2."<<endl;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            BACKGROUND_U = AOS(u1, u1*u2, 0.5*u1*u2*u2 + u3/(gamma_adiabat[0]-1.) );

     
            /* mdot boundaries
             * needs fixing... Should just be specified as a BoundaryType::fixed boundary
            if(boundaries_number == 3) {
                mdot = read_parameter_from_file<double>(filename,"ACCRETION_RATE", debug).value; //Accretion rate in numerical units
    
                double rho0 = BACKGROUND_U.u1;
                double e0   = BACKGROUND_U.u3 - 0.5 * BACKGROUND_U.u2 * BACKGROUND_U.u2 / rho0;
    
                SHOCK_TUBE_UR = AOS(rho0, -mdot, e0 + 0.5 * mdot * mdot / rho0 );

                boundary_left = BoundaryType::reflecting ;
                boundary_right = BoundaryType::fixed ;
            }
            */
            
            
            cout<<"Problem 2 pos3."<<endl;
            
            initialize_background(BACKGROUND_U);
            init_grav_pot();
            
            //
            // If we want to initialize a hydrostatic atmosphere, then we overwrite the so far given density/pressure data and set every velocity to 0
            //
            if(init_static_atmosphere == 1) {
                initialize_hydrostatic_atmosphere_nonuniform();
            }

            if(debug > 0) 
                cout<<"Successfully initialized problem 2."<<endl;
        }
        else 
            initialize_default_test();
        
        
        USE_WAVE = read_parameter_from_file<int>(filename,"WAVE_USE", debug, 0).value; 
        if(USE_WAVE==1) {
            WAVE_AMPLITUDE = read_parameter_from_file<double>(filename,"WAVE_AMPLITUDE", debug).value; 
            WAVE_PERIOD    = read_parameter_from_file<double>(filename,"WAVE_PERIOD", debug).value; 
        }
        
        // TEMPERATURE_BUMP_STRENGTH = read_parameter_from_file<double>(filename,"TEMPERATURE_BUMP_STRENGTH", debug).value; 
        
        if(debug > 0)
            cout<<"Ended init."<<endl;
}

//
//
//
void hydro_run::initialize_hydrostatic_atmosphere_nonuniform() {
    
    cout<<"ATTENTION: Initializing nonuniform hydrostatic atmosphere and overwriting prior initial values."<<endl;
    
    long double temp_rhofinal;
    long double factor_inner, factor_outer;
    //long double factor_grav;
    //long double delta_phi;
    long double T_inner;
    long double T_outer;
    long double metric_inner;
    long double metric_outer;
    //
    // Start with the outermost cell and build up a hydrostatic atmosphere
    // Fulfilling KKM16, Eqn. 17
    //
    //long double residual  = 0.;
    
    //
    // First, initialize (adiabatic) temperature
    //
    for(int i=num_cells+1; i>=0; i--) {
            //temperature[i] = planet_mass / x_i12[i] / (cv * gamma_adiabat) + 1.;
            temperature[i] = - 1.0 * phi[i] / (cv[i] * gamma_adiabat[i]) + temperature[num_cells+1];
            //temperature[i] = 100.;
        
            //Add temperature bumps and troughs
            temperature[i] += TEMPERATURE_BUMP_STRENGTH * 4. * exp( - pow(x_i12[i] - 1.e-1 ,2.) / (0.1) );
            
            temperature[i] -= TEMPERATURE_BUMP_STRENGTH * 15. * exp( - pow(x_i12[i] - 7.e-3 ,2.) / (1.5e-3) );

    }
    
    //Last normal cell has to be awkwardly initialized
    
    //At this point, the right ghost cell is already initialized, so we can just build up a hydrostatic state from there
    for(int i=num_cells; i>=0; i--)  {
        
        //
        // Construct next density as to fulfil the hydrostatic condition
        //
        T_outer = temperature[i+1];
        T_inner = temperature[i];  

        factor_outer = (gamma_adiabat[i+1]-1.) * cv[i+1] * T_outer; //TODO: Replace with T_outer for non-isothermal EOS
        factor_inner = (gamma_adiabat[i]  -1.) * cv[i]   * T_inner; //TODO: Replace with T_inner for non-isothermal EOS
        metric_outer = (phi[i+1] - phi[i]) * omegaplus[i+1] * dx[i+1] / (dx[i+1] + dx[i]);
        metric_inner = (phi[i+1] - phi[i]) * omegaminus[i]  * dx[i]   / (dx[i+1] + dx[i]);
        
        temp_rhofinal = u[i+1].u1 * (factor_outer + metric_outer) / (factor_inner - metric_inner);
        
        u[i] = AOS(temp_rhofinal, 0., cv[i] * temp_rhofinal * T_inner) ;
        
        //
        // Debug info
        // 
        //if( (i==20 || i== num_cells-20) && debug >= 0) {
        //if( i>0 ) {
        if(1==0) {
            char a;
            cout.precision(16);
            cout<<"In INIT STATIC i="<<i<<endl;
            cout<<"In hydostatic init: factor_outer/metric_outer = "<< factor_outer / metric_outer << " factor_inner/metric_inner = "<< factor_inner / metric_inner <<endl;
            cout<<"factor_outer+metric_outer = "<< (factor_outer + metric_outer) << " factor_inner-metric_inner = "<<( factor_inner - metric_inner) <<endl;
            cout<<"metric_outer = "<< metric_outer << " metric_inner = "<<metric_inner <<endl;
            cout<<"factor_outer = "<< factor_outer << " factor_inner = "<<factor_inner <<endl;
            //cout<<"In hydostatic init: factor_dens = "<< (2.* factor_outer / delta_phi + 1.) / (2. * factor_inner / delta_phi - 1.) <<endl;
            cout<<"Ratio of densities inner/outer = "<< temp_rhofinal/u[i+1].u1 <<endl;
            cout<<"Ratio of temperatures inner/outer = "<<T_inner/T_outer<<" t_inner ="<<T_inner<<" t_outer ="<<T_outer<<endl;
            cout<<"Ratio of pressures inner/outer = "<<cv[i] * temp_rhofinal * T_inner /u[i+1].u3<<endl;
            //tempgrad = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            //tempgrad2 = ((gamma_adiabat-1.)*cv*T_outer + 0.5 * delta_phi) * u[i+1].u1 - ((gamma_adiabat-1.)*cv*T_inner  + 0.5 * delta_phi ) * u[i].u1 ;
            //residual = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            //cout<<"pressure diff "<<(gamma_adiabat-1.)*u[i+1].u3 - (gamma_adiabat-1.)*(u[i].u3)<<endl;
            cout<<"pressure diff "<<( ((gamma_adiabat[i+1]-1.)*cv[i+1]*u[i+1].u1*T_outer - (gamma_adiabat[i]-1.)*cv[i]*u[i].u1*T_inner)/dx[i])<<endl;
            //cout<<"density sum with potential "<<(0.5 * (u[i].u1 + u[i+1].u1) * delta_phi/dx[i])<<endl;
            cout<<"density diff "<<(u[i].u1 - u[i+1].u1)<<endl;
            cout<<"density sum " <<(u[i].u1 + u[i+1].u1)<<endl;
            //cout<<"residual = "<<residual<<endl;
            
            //cout<<"residual2 = "<<residual<<endl;
            //cout<<"sum of hydrostatic gradients = "<<tempgrad<<endl;
            //cout<<"sum2 of hydrostatic gradients = "<<tempgrad2<<endl;
            cout<<"Resulting density == "<<temp_rhofinal<<endl;
            cout<<" density before "<<u[i+1].u1<<endl;
            cin>>a;
        }
    }
    
    cout<<"Ended hydrostatic construction."<<endl;

}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void hydro_run::initialize_default_test() {
    
    cout<<"Initialized default simulation."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = AOS(1,1,1);
}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void hydro_run::initialize_shock_tube_test(const AOS &leftval,const AOS &rightval) {
    
    cout<<"Initialized shock tube test."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) {
            
        if(x_i12[i] < SHOCK_TUBE_MID) 
            u[i] = leftval;
        else
            u[i] = rightval;
    }
}

//
// Initial conditions for one constant background state
//
void hydro_run::initialize_background(const AOS &background) {
    
    cout<<"Initialized constant background state. state = "<<background.u1<<" "<<background.u2<<" "<<background.u3<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = background;
}


void hydro_run::apply_boundary_left() {
    double E_kinetic, pressure_active, pressure_bound ;
    switch(boundary_left) {
        case BoundaryType::user:
            user_boundary_left();
            break;
        case BoundaryType::open:
            u[0]            = u[1]; 
            E_kinetic       = 0.5 * u[1].u2 * u[1].u2 / u[1].u1 ;
            pressure_active = (gamma_adiabat[1]-1.) * (u[1].u3 - E_kinetic);
            // Hydrostatic pressure extrapolation 
            pressure_bound = pressure_active - u[1].u1 * (phi[0] - phi[1]) ;
            pressure_bound = std::max(pressure_bound, 0.0) ;
            u[0].u3        = pressure_bound/(gamma_adiabat[0]-1.) + E_kinetic ;
            break ;
        case BoundaryType::reflecting:
            u[0]     = u[1]; 
            u[0].u2 *= -1;
            phi[0]   = phi[1] ;
            break;
        case BoundaryType::fixed:
            u[0]     = SHOCK_TUBE_UL;
            break;
        case BoundaryType::periodic:
            u[0]     = u[num_cells];
            phi[0]   = phi[num_cells] ;
            break;
    }
}

void hydro_run::apply_boundary_right() {
    double E_kinetic, pressure_active, pressure_bound ;
    switch(boundary_right) {
        case BoundaryType::user:
            user_boundary_right();
            break;
        case BoundaryType::open:
            u[num_cells+1]    = u[num_cells]; 
            E_kinetic         = 0.5 * u[num_cells].u2 * u[num_cells].u2 / u[num_cells].u1 ;
            pressure_active   = (gamma_adiabat[num_cells]-1.) * (u[num_cells].u3 - E_kinetic);
            // Hydrostatic pressure extrapolation 
            pressure_bound    = pressure_active - u[num_cells].u1 * (phi[num_cells+1] - phi[num_cells]) ;
            pressure_bound    = std::max(pressure_bound, 0.0) ;
            u[num_cells+1].u3 = pressure_bound / (gamma_adiabat[num_cells+1]-1.) + E_kinetic ;
            break ;
        case BoundaryType::reflecting:
            u[num_cells+1]     = u[num_cells]; 
            u[num_cells+1].u2 *= -1;
            phi[num_cells+1]   = phi[num_cells] ;
            break;
        case BoundaryType::fixed:
            u[num_cells+1]     = SHOCK_TUBE_UR;
            break;
        case BoundaryType::periodic:
            u[num_cells+1]     = u[1];
            phi[num_cells+1]   = phi[1] ;
            break;
    }
}

//
// Creates a wave at the inner boundary supposedly propagating upwards in the gravity well
//
void hydro_run::add_wave(double globalTime, double WAVE_PERIOD, double WAVE_AMPLITUDE)  {
    
    
    u[0].u2 = u[0].u1 * WAVE_AMPLITUDE * std::sin(12.*M_PI*globalTime/WAVE_PERIOD);
    u[1].u2 = u[0].u2;
    
}
hydro_run::~hydro_run() {
    
    cout<<"The destructor is called, but its empty so far."<<endl;
    
}

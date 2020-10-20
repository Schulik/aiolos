#include "aiolos.h"

//extern AOS* init_AOS(int num);

////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  CLASS SIMULATION
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c_Sim::c_Sim(string filename, string speciesfile, int debug) {

        this->debug = debug ;
        
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
        geometry = read_parameter_from_file<Geometry>(filename, "PARI_GEOMETRY", debug, Geometry::cartesian).value;
        order = read_parameter_from_file<IntegrationType>(filename, "PARI_ORDER", debug, IntegrationType::second_order).value;

        if (order == IntegrationType::first_order)
            num_ghosts = 1 ;
        else 
            num_ghosts = 2 ;

        if(type_of_grid == 0) {
            num_cells        = (int)((domain_max - domain_min)/dx0);
            cout<<"Domain specifics:  "<<domain_min<<" | . . . "<<num_cells<<" uniform cells . . . | "<<domain_max<<endl;
        }
        else {
            num_cells        = (int) ( (log10f(domain_max) - log10f(domain_min)) * cells_per_decade );
            dx0              = domain_min;
            cout<<"Domain specifics:  "<<domain_min<<" | . .  .  "<<num_cells<<" nonuniform cells .       .             .     | "<<domain_max<<endl;
        }
        
        //
        // Check that cell number looks fine. Shoutout if its not.
        //
        if( !(num_cells > 0 && num_cells < 999999) ) {
            std::stringstream err ;
            err<<"WARNING! Something seems wrong with the number of cells required. num_cells = "<<num_cells<<endl;
            throw std::invalid_argument(err.str()) ;
        }
            
        if(debug > 0) cout<<"Init: Finished reading grid parameters."<<endl;
        
        //dt          = read_parameter_from_file<double>(filename,"PARI_TIME_DT", debug).value;
        cflfactor   = read_parameter_from_file<double>(filename,"PARI_CFLFACTOR", debug, 1.).value;
        t_max       = read_parameter_from_file<double>(filename,"PARI_TIME_TMAX", debug).value;
        output_time = read_parameter_from_file<double>(filename,"PARI_TIME_OUTPUT", debug).value; 
        globalTime = 0.0;    
        timecount = 0;
        
        if(debug > 0) cout<<"Init: Finished reading time parameters."<<endl;

        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //    Control parameters for users
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        simname = filename;
        
        problem_number    = read_parameter_from_file<int>(filename,"PARI_PROBLEM_NUMBER", debug).value;
        use_self_gravity  = read_parameter_from_file<int>(filename,"PARI_SELF_GRAV_SWITCH", debug, 0).value;
        use_linear_gravity= read_parameter_from_file<int>(filename,"PARI_LINEAR_GRAV", debug, 0).value;
        use_rad_fluxes    = read_parameter_from_file<int>(filename,"PARI_USE_RADIATION", debug, 0).value;
        if(problem_number==2) {
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug, 1).value; //in Earth masses
            planet_position = read_parameter_from_file<double>(filename,"PARI_PLANET_POS", debug, 0.).value;  //inside the simulation domain
            rs              = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_LENGTH", debug, 0.).value; //Gravitational smoothing length in hill 
            rs_time         = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_TIME", debug, 0.).value; //Time until we reach rs starting at rs_at_moment
            rs_at_moment    = 0.2;
            
        } else {
            planet_mass = planet_position = rs = rs_time = 0 ;
            rs_at_moment = 0.2 ;
        }
        
        if(use_self_gravity)
            cout<<"WARNING: Self-gravity switched ON !!!"<<endl;
        
        if(debug > 0) cout<<"Init: Finished reading gravity parameters"<<endl;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Simulation data: Grid and variables
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        // Increment num cells so that we have enough space for the ghosts
        num_cells += 2*(num_ghosts-1) ;

        x_i        = np_zeros(num_cells+1);		//The cell boundaries
        x_i12      = np_zeros(num_cells+2);		//The cell mid positions
        surf       = np_zeros(num_cells+1);     // intercell surfaces
        vol        = np_zeros(num_cells+2);     // cell volumes
        dx         = np_zeros(num_cells+2);
        omegaplus  = np_zeros(num_cells+2);
        omegaminus = np_zeros(num_cells+2);
        source_pressure_prefactor_left    = np_zeros(num_cells+2);
        source_pressure_prefactor_right   = np_zeros(num_cells+2);
        enclosed_mass   = np_zeros(num_cells+2);
        phi             = np_zeros(num_cells+2);
        
        if(debug > 0) cout<<"Init: Setup grid memory"<<endl;
        
        //
        // Compute cell wall boundaries
        //         First and last two cells near boundaries are uniform
        //
        if(type_of_grid==1) {

            double dlogx = pow(10.,1./cells_per_decade);
            
            x_i[0] = domain_min / std::pow(dlogx, (num_ghosts-1)) ;
            for(int i=1; i<= num_cells; i++) {
                x_i[i]   = x_i[i-1] * dlogx;
            }
            
            //Assign the last boundary as domain maximum as long as nonuniform grid is in the test-phase
            domain_max = x_i[num_cells];
            cout<<"We have a changed DOMAIN MAX = "<<domain_max<<endl;
                
        }
        //Uniform grid
        else {
            x_i[0] = domain_min - dx0*(num_ghosts-1) ;
            for(int i=1; i<= num_cells; i++) {
                x_i[i] = x_i[i-1] + dx0;
            }
        }
       
        //
        // Surfaces and volumes for all geometries
        //
        
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
            x_i12[i] = 0.5 * (x_i[i] + x_i[i-1]);
        }
        
        //Differences
        for(int i=1; i<num_cells+1; i++) {
            dx[i]  = x_i[i]-x_i[i-1];
        }
        //Ghost dxs, those have only to be generated such that hydrostatic eq is preserved, they don't have to be physical
        dx[0] = dx[1]*dx[1]/dx[2] ;
        dx[num_cells+1] = dx[num_cells]*dx[num_cells]/dx[num_cells-1];
        
        //Ghost cells
        x_i12[0]           = x_i[0] - dx[0]/2 ;
        x_i12[num_cells+1] = x_i[num_cells] + dx[num_cells+1]/2; 
        
        //
        // Grid generation done. Now do a few simple checks on the generated grid values
        //TODO: Expand this with a few more sensible tests, NaN checks etc.
        
        int all_grid_ok = 1, broken_index;
        double broken_value;
        for(int i = 0; i <= num_cells+1; i++)
            if( dx[i] < 0 ) {
                all_grid_ok = 0;
                broken_index = i;
                broken_value = dx[i];
            }
        
        
        if(debug > 0) cout<<"Init: Finished grid setup"<<endl;
        if(!all_grid_ok) {
            std::stringstream err ;
            err <<"WARNING! Problem with grid values. Example broken value, dx["<<broken_index<<"]="<<broken_value<<endl;
            err<<" FIRST AND LAST DX = "<<dx[0]<<" / "<<dx[num_cells+1]<<endl;
            throw std::runtime_error(err.str()) ;
        }
            
        
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
        omegaminus[0] = 2*omegaminus[1] - omegaminus[2] ;
        omegaplus[0]  = 2*omegaplus[1] - omegaplus[2];
        
        omegaminus[num_cells+1] = 2*omegaminus[num_cells] - omegaminus[num_cells-1];
        omegaplus[num_cells+1]  = 2*omegaplus[num_cells] - omegaplus[num_cells-1];
        
        for(int i=1; i<num_cells+1; i++) {
            source_pressure_prefactor_left[i]  = (surf[i-1]/vol[i] - 1./dx[i]); 
            source_pressure_prefactor_right[i] = (surf[i]  /vol[i] - 1./dx[i]); 
        }
        
        if(debug > 0) {
                cout<<"DEBUGLEVEL 1: Samples of generated grid values"<<endl;
                cout<<"First cell coordinates: |<--"<<x_i[0]<<" // "<<x_i12[0]<<" // "<<x_i[1]<<"-->|"<<endl;
                cout<<"Last cell coordinates:  |<--"<<x_i[num_cells-1]<<" // "<<x_i12[num_cells-1]<<" // "<<x_i[num_cells]<<"-->|"<<endl;
        } 
        
        
        if(debug > 0) cout<<"Init: Setup metric factors"<<endl;
        
        //
        //
        // Initialize gravity first without any species. Hydrostatic construction needs gravity to start.
        //
        //
        
        init_grav_pot();
        
        //
        //
        // Create and initialize the species
        //
        //
        num_species = read_parameter_from_file<int>(filename,"PARI_NUM_SPECIES", debug, 1).value;
        
        if(debug > 0) cout<<"Init: About to setup species with num_species = "<<num_species<<endl;
        
        species.reserve(num_species);
        
        for(int s = 0; s < num_species; s++) {
            species.push_back(c_Species(this, filename, speciesfile, s, debug)); // Here we call the c_Species constructor
        }
        
        //
        //
        // After successful init, compute all primitive quantities and determine first timestep
        //
        //
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure(species[s].u);
        
        dt = get_cfl_timestep();
        
        if(debug > 0) cout<<"Init: Finished Init."<<endl;
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

c_Species::c_Species(c_Sim *base_simulation, string filename, string species_filename, int species_index, int debug) {
    
        base               = base_simulation;
        this_species_index = species_index;
        num_cells          = base->num_cells;
        this->debug        = debug;
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Begin"<<endl;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Physical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        boundary_left  = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_LEFT", debug, BoundaryType::fixed).value;
        boundary_right = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_RIGHT", debug, BoundaryType::fixed).value;
        const_T_space  = read_parameter_from_file<double>(filename,"PARI_CONST_TEMP", debug, 1.).value;
        TEMPERATURE_BUMP_STRENGTH = read_parameter_from_file<double>(filename,"TEMPERATURE_BUMP_STRENGTH", debug, 0.).value; 
        
        //
        // Determine which equation of states we are using
        //
        eos_pressure_type       = EOS_pressure_type::adiabatic;       //Read parameter from...
        eos_internal_energy_type = EOS_internal_energy_type::thermal;
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Finished reading boundaries."<<endl;
        if(debug > 0) cout<<"         Boundaries used in species["<<name<<"]: "<<boundary_left<<" / "<<boundary_right<<endl;

        //Read EOS specifications from file
        if(eos_pressure_type == EOS_pressure_type::tabulated) 
            read_tabulated_eos_data_pressure("eos_pressure.input");
            
        if(eos_internal_energy_type == EOS_internal_energy_type::tabulated) 
            read_tabulated_eos_data_eint("eos_internal_energy.input");
        
        if(base->use_rad_fluxes == 1)
            const_opacity   = read_parameter_from_file<double>(filename,"PARI_CONST_OPAC", debug, 1.).value;
        
        
        //
        // //This provides important stuff like, const_cv, const_gamma_adiabat for this species
        //
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Reading species data..."<<endl;
        read_species_data(species_filename, species_index); 
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Done reading species data."<<endl;
        
        
        u               = init_AOS(num_cells+2); //Conserved hyperbolic variables: density, mass flux, energy density
        dudt[0]         = init_AOS(num_cells+2);
        dudt[1]         = init_AOS(num_cells+2);
        source          = init_AOS(num_cells+2);  
        source_pressure = init_AOS(num_cells+2);
        flux            = init_AOS(num_cells+1);
                
        gamma_adiabat  = np_somevalue(num_cells+2, const_gamma_adiabat);
        cv             = np_somevalue(num_cells+2, const_cv);
        opacity        = np_somevalue(num_cells+2, const_opacity);
        temperature    = np_somevalue(num_cells+2, const_T_space);
        opticaldepth   = np_zeros(num_cells+2);
        radiative_flux = np_zeros(num_cells+1);
        
        pressure        = np_zeros(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        pressure_l      = np_zeros(num_cells+2); 
        pressure_r      = np_zeros(num_cells+2); 
        internal_energy = np_zeros(num_cells+2);
        speed           = np_zeros(num_cells+2);
        cs              = np_zeros(num_cells+2);
        
        timesteps    = np_zeros(num_cells+2);		    
        timesteps_cs = np_zeros(num_cells+2);	
        finalstep    = np_zeros(num_cells+2);
        
        
        
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Initialize/readin parameters
        //
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        //
        // Problem 1: two states, left and right, separated at a mid-position
        //
        if(base->problem_number == 1) {
            
            double u1l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1L", debug, 1.).value;
            double u2l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2L", debug, 0.).value;
            double u3l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3L", debug, 1.).value;
            
            double u1r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1R", debug, 1.).value;
            double u2r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2R", debug, 0.).value;
            double u3r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3R", debug, 0.1).value;
            
            SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug, 0.5).value;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            SHOCK_TUBE_UL = AOS(u1l, u1l*u2l, 0.5*u1l*u2l*u2l + u3l/(gamma_adiabat[0]-1.) );
            SHOCK_TUBE_UR = AOS(u1r, u1r*u2r, 0.5*u1r*u2r*u2r + u3r/(gamma_adiabat[num_cells]-1.));
            
            initialize_shock_tube_test(SHOCK_TUBE_UL, SHOCK_TUBE_UR);
            
            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 1."<<endl;
        }
        //
        // Problem 2: A constant background state, with a planet embedded into it, free to evolve as physics dictates
        //
        else if(base->problem_number == 2) {
            
            double u1 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1", debug, 1.).value;
            double u2 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2", debug, 0.).value;
            double u3 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3", debug, 1.).value;
            
            init_static_atmosphere = read_parameter_from_file<int>(filename,"PARI_INIT_STATIC", debug, 1).value; //Yesno, isothermal at the moment
            
            if(debug > 0) cout<<"        Species["<<species_index<<"] problem 2 init, pos 1."<<endl;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            //BACKGROUND_U = AOS(u1, u1*u2, 0.5*u1*u2*u2 + u3/(gamma_adiabat[0]-1.) );
            
            BACKGROUND_U = AOS(u1 * initial_fraction, u1*u2 * initial_fraction, 0.5*u1*u2*u2 * initial_fraction  + u3/(gamma_adiabat[0]-1.) );

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
            
            initialize_background(BACKGROUND_U);
            
            //
            // If we want to initialize a hydrostatic atmosphere, then we overwrite the so far given density/pressure data and set every velocity to 0
            //
            if(init_static_atmosphere == 1) {
                initialize_hydrostatic_atmosphere();
            }

            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 2."<<endl;
        }
        else 
            initialize_default_test();
        
        
        USE_WAVE = read_parameter_from_file<int>(filename,"WAVE_USE", debug, 0).value; 
        if(USE_WAVE==1) {
            WAVE_AMPLITUDE = read_parameter_from_file<double>(filename,"WAVE_AMPLITUDE", debug).value; 
            WAVE_PERIOD    = read_parameter_from_file<double>(filename,"WAVE_PERIOD", debug).value; 
        }
        
        // Apply boundary conditions
        apply_boundary_left(u) ;
        apply_boundary_right(u) ;
        
        if(debug > 0) cout<<"        Species["<<species_index<<"]: Init done."<<endl;
    
}



void c_Species::initialize_hydrostatic_atmosphere() {
    
    if(debug > 0) cout<<"            ATTENTION: Initializing hydrostatic construction for species "<<name<<" and overwriting prior initial values."<<endl;
    
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
            temperature[i] = - 1.0 * base->phi[i] / (cv[i] * gamma_adiabat[i]) + const_T_space;
            //temperature[i] = 100.;
        
            //Add temperature bumps and troughs
            temperature[i] += TEMPERATURE_BUMP_STRENGTH * 4.  * exp( - pow(base->x_i12[i] - 1.e-1 ,2.) / (0.1) );
            temperature[i] -= TEMPERATURE_BUMP_STRENGTH * 15. * exp( - pow(base->x_i12[i] - 7.e-3 ,2.) / (1.5e-3) );

    }
    
    //At this point, the right ghost cell is already initialized, so we can just build up a hydrostatic state from there
    for(int i=num_cells; i>=0; i--)  {
        
        //
        // Construct next density as to fulfil the hydrostatic condition
        //
        T_outer = temperature[i+1];
        T_inner = temperature[i];  

        factor_outer = (gamma_adiabat[i+1]-1.) * cv[i+1] * T_outer; //TODO: Replace with T_outer for non-isothermal EOS
        factor_inner = (gamma_adiabat[i]  -1.) * cv[i]   * T_inner; //TODO: Replace with T_inner for non-isothermal EOS
        metric_outer = (base->phi[i+1] - base->phi[i]) * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
        metric_inner = (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
        
        temp_rhofinal = u[i+1].u1 * (factor_outer + metric_outer) / (factor_inner - metric_inner);
        
        u[i] = AOS(temp_rhofinal, 0., cv[i] * temp_rhofinal * T_inner) ;
        
        //
        // Debug info
        // 
        //if( (i==20 || i== num_cells-20) && debug >= 0) {
        //if( i>0 ) {
        if(temp_rhofinal < 0) {
            char a;
            cout.precision(16);
            cout<<"NEGATIVE DENSITY IN INIT HYDROSTATIC i="<<i<<endl;
            if((factor_inner - metric_inner) < 0) {
                
                cout<<"     Negative denominator detected. Product (gamma-1)*cv*T_inner is too small for chosen discretization."<<endl;
                
            }
            cout<<endl;
            cout<<"     metric_outer = "<< metric_outer << " metric_inner = "<<metric_inner <<endl;
            cout<<"     factor_outer = "<< factor_outer << " factor_inner = "<<factor_inner <<endl;
            cout<<"     factor_outer+metric_outer = "<< (factor_outer + metric_outer) << " factor_inner-metric_inner = "<<( factor_inner - metric_inner) <<endl;
            cout<<"     RATIO = "<< ((factor_outer + metric_outer)/ (factor_inner - metric_inner)) <<endl;
            //cout<<"In hydostatic init: factor_dens = "<< (2.* factor_outer / delta_phi + 1.) / (2. * factor_inner / delta_phi - 1.) <<endl;
            cout<<"Ratio of densities inner/outer = "<< temp_rhofinal/u[i+1].u1 <<endl;
            cout<<"Ratio of temperatures inner/outer = "<<T_inner/T_outer<<" t_inner ="<<T_inner<<" t_outer ="<<T_outer<<endl;
            cout<<"Ratio of pressures inner/outer = "<<cv[i] * temp_rhofinal * T_inner /u[i+1].u3<<endl;
            //tempgrad = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            //tempgrad2 = ((gamma_adiabat-1.)*cv*T_outer + 0.5 * delta_phi) * u[i+1].u1 - ((gamma_adiabat-1.)*cv*T_inner  + 0.5 * delta_phi ) * u[i].u1 ;
            //residual = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            //cout<<"pressure diff "<<(gamma_adiabat-1.)*u[i+1].u3 - (gamma_adiabat-1.)*(u[i].u3)<<endl;
            //cout<<"pressure diff "<<( ((gamma_adiabat[i+1]-1.)*cv[i+1]*u[i+1].u1*T_outer - (gamma_adiabat[i]-1.)*cv[i]*u[i].u1*T_inner)/base->dx[i])<<endl;
            //cout<<"density sum with potential "<<(0.5 * (u[i].u1 + u[i+1].u1) * delta_phi/dx[i])<<endl;
            //cout<<"density diff "<<(u[i].u1 - u[i+1].u1)<<endl;
            //cout<<"density sum " <<(u[i].u1 + u[i+1].u1)<<endl;
            //cout<<"residual = "<<residual<<endl;
            
            //cout<<"residual2 = "<<residual<<endl;
            //cout<<"sum of hydrostatic gradients = "<<tempgrad<<endl;
            //cout<<"sum2 of hydrostatic gradients = "<<tempgrad2<<endl;
            cout<<"Resulting density == "<<temp_rhofinal<<endl;
            cout<<" density before "<<u[i+1].u1<<endl;
            cin>>a;
        }
    }
    
    
    cout<<"            Ended hydrostatic construction for species "<<name<<endl;
    //TODO:Give some measure if hydrostatic construction was also successful, and give warning if not.

}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void c_Species::initialize_default_test() {
    
    if(debug > 0) cout<<"            Initialized default simulation."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = AOS(1,1,1);
}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void c_Species::initialize_shock_tube_test(const AOS &leftval,const AOS &rightval) {
    
    if(debug > 0) cout<<"            Initialized shock tube test."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) {
        
        if(base->x_i12[i] < SHOCK_TUBE_MID) 
            u[i] = leftval;
        else
            u[i] = rightval;
    }
}

//
// Initial conditions for one constant background state
//
void c_Species::initialize_background(const AOS &background) {
    
    if(debug > 0) cout<<"            Initialized constant background state. state = "<<background.u1<<" "<<background.u2<<" "<<background.u3<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = background;
}


void c_Species::apply_boundary_left(std::vector<AOS>& u) {
    double E_kinetic, pressure_active, pressure_bound, dphi ;
    int num_ghosts = base->num_ghosts;
    int Ncell = num_cells - 2*(num_ghosts-1) ; // Correct for fact we increased num_cells
    switch(boundary_left) {
        case BoundaryType::user:
            user_boundary_left(u);
            break;
        case BoundaryType::open:
            for (int i=num_ghosts; i > 0; i--) {
                u[i-1]            = u[i]; 
                E_kinetic       = 0.5 * u[i].u2 * u[i].u2 / u[i].u1 ;
                pressure_active = (gamma_adiabat[i]-1.) * (u[i].u3 - E_kinetic);
                // Hydrostatic pressure extrapolation 
                dphi = (base->phi[i] - base->phi[i-1]) / (base->dx[i] + base->dx[i-1]) ;
                dphi *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                pressure_bound = pressure_active +  u[i].u1 * dphi ;
                pressure_bound = std::max(pressure_bound, 0.0) ;
                u[i-1].u3        = pressure_bound/(gamma_adiabat[i-1]-1.) + E_kinetic ;
            }
            break ;
        case BoundaryType::reflecting:
            for (int i=0; i < num_ghosts; i++) {
                int igh =  num_ghosts-1 -i;
                int iact = num_ghosts   +i;
                u[igh]     = u[iact]; 
                u[igh].u2 *= -1;
                base->phi[igh]   = base->phi[iact] ;
            }
            break;
        case BoundaryType::fixed:
            for (int i=0; i < num_ghosts; i++)
                u[i]     = SHOCK_TUBE_UL;
            break;
        case BoundaryType::periodic:
            for (int i=0; i < num_ghosts; i++) {
                int iact = Ncell + num_ghosts +i;
                u[i] = u[iact];
                base->phi[i]   = base->phi[iact] ;
            }
            break;
    }
}

void c_Species::apply_boundary_right(std::vector<AOS>& u) {
    double E_kinetic, pressure_active, pressure_bound, dphi ;
    int num_ghosts = base->num_ghosts;
    int Ncell = num_cells - 2*(num_ghosts-1) ;
    switch(boundary_right) {
        case BoundaryType::user:
            user_boundary_right(u);
            break;
        case BoundaryType::open:
            for (int i=Ncell+num_ghosts; i < Ncell+2*num_ghosts; i++) {
                u[i]    = u[i-1]; 
                E_kinetic         = 0.5 * u[i-1].u2 * u[i-1].u2 / u[i-1].u1 ;
                pressure_active   = (gamma_adiabat[i-1]-1.) * (u[i-1].u3 - E_kinetic);
                // Hydrostatic pressure extrapolation 
                dphi = (base->phi[i] - base->phi[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                dphi *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                pressure_bound = pressure_active -  u[i-1].u1 * dphi ;
                pressure_bound    = std::max(pressure_bound, 0.0) ;
                u[i].u3 = pressure_bound / (gamma_adiabat[i]-1.) + E_kinetic ;
            }
            break ;
        case BoundaryType::reflecting:
            for (int i=0; i < num_ghosts; i++) {
                int iact = Ncell + num_ghosts-1 -i ;
                int igh = Ncell + num_ghosts +i;
                u[igh]     = u[iact]; 
                u[igh].u2 *= -1;
                base->phi[igh]   = base->phi[iact] ; //TODO: This will be called many times, too many!
            }
            break;
        case BoundaryType::fixed:
            for (int i=Ncell+ num_ghosts; i < Ncell+2*num_ghosts; i++)
                u[i]     = SHOCK_TUBE_UR;
            break;
        case BoundaryType::periodic:
            for (int i=Ncell+num_ghosts; i < Ncell+2*num_ghosts; i++) {
                int iact = num_ghosts + i;
                u[i]     = u[iact];
                base->phi[i]   = base->phi[iact] ; //TODO: This will be called too many times!
            }
            break;
    }
}

//
// Creates a wave at the inner boundary supposedly propagating upwards in the gravity well
//
void c_Species::add_wave(std::vector<AOS>& u, double globalTime)  {
    
    u[0].u2 = u[0].u1 * WAVE_AMPLITUDE * std::sin(12.*M_PI*globalTime/WAVE_PERIOD);
    u[1].u2 = u[0].u2;
    
}
c_Sim::~c_Sim() {
    //Empty
}

c_Species::~c_Species() {
    //Empty
}

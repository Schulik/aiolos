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

    
        if(debug > 0) cout<<"Init position 0."<<endl;
        
        this->debug = debug ;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //  Numerical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double dx0;
        type_of_grid     = read_parameter_from_file<int>(filename,"PARI_GRID_TYPE", debug, 0).value;
        domain_min       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MIN", debug).value;
        domain_max       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MAX", debug).value;
        geometry = read_parameter_from_file<Geometry>(filename, "PARI_GEOMETRY", debug, Geometry::cartesian).value;
        order = read_parameter_from_file<IntegrationType>(filename, "PARI_ORDER", debug, IntegrationType::second_order).value;

        if(debug > 0) cout<<"Using integration order "<<order<<" while second order would be "<<IntegrationType::second_order<<endl;
        
        if (order == IntegrationType::first_order)
            num_ghosts = 1 ;
        else 
            num_ghosts = 2 ;

        if(type_of_grid == 0) {
            dx0              = read_parameter_from_file<double>(filename,"PARI_DOMAIN_DX", debug).value;
            num_cells        = (int)((domain_max - domain_min)/dx0);
            cout<<"Domain specifics:  "<<domain_min<<" | . . . "<<num_cells<<" uniform cells . . . | "<<domain_max<<endl;
        }
        else {
            cells_per_decade = read_parameter_from_file<double>(filename,"PARI_CELLS_PER_DECADE", debug).value;
            num_cells        = (int) ( (log10f(domain_max) - log10f(domain_min)) * cells_per_decade );
            dx0              = domain_min;
            cout<<"Domain specifics:  "<<domain_min<<" | . .  .  "<<num_cells<<" nonuniform cells .       .             .     | "<<domain_max<<endl;
        }
        
        //
        // Check that cell number looks fine. Shoutout if its not.
        //
        if( !(num_cells > 0 && num_cells < 999999) ) {
            std::stringstream err ;
            err<<"    WARNING! Something seems wrong with the number of cells required. num_cells = "<<num_cells<<endl;
            throw std::invalid_argument(err.str()) ;
        }
            
        if(debug > 0) cout<<"Init: Finished reading grid parameters."<<endl;
        
        
        cflfactor   = read_parameter_from_file<double>(filename,"PARI_CFLFACTOR", debug, 1.).value;
        t_max       = read_parameter_from_file<double>(filename,"PARI_TIME_TMAX", debug).value;
        output_time = read_parameter_from_file<double>(filename,"PARI_TIME_OUTPUT", debug).value; 
        monitor_time = read_parameter_from_file<double>(filename,"PARI_TIME_DT", debug).value;
        globalTime = 0.0;    
        timecount = 0;
        
        if(t_max / monitor_time > 1e6) {
            char b;
            cout<<"    WARNING! Specified t_max and t_dt will result in "<<(t_max / monitor_time)<<" entries in the monitor file."<<endl<<"    This might create a large file and/or impact performance significantly. Press the any key to continue."<<endl;
            cin>>b;
        }
        
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
        init_wind         = read_parameter_from_file<int>(filename,"PARI_INIT_WIND", debug, 0).value;
        alpha_collision   = read_parameter_from_file<double>(filename,"PARI_ALPHA_COLL", debug, 0).value;
        collision_model   = read_parameter_from_file<char>(filename,"PARI_COLL_MODEL", debug, 'C').value;
        friction_solver   = read_parameter_from_file<int>(filename,"FRICTION_SOLVER", debug, 0).value;
        
        if(problem_number == 2)
            monitor_output_index = num_cells/2; //TODO: Replace with index of sonic radius for dominant species?
        else
            monitor_output_index = 1;
        
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

        x_i        = np_zeros(num_cells+1);		// The cell boundaries
        x_i12      = np_zeros(num_cells+2);		// The cell mid positions
        x_iVC      = np_zeros(num_cells+2);		// The cell volumetric centres
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
        // Differences, surfaces and volumes for all geometries
        //
        
        //Differences
        for(int i=1; i<num_cells+1; i++) {
            dx[i]  = x_i[i]-x_i[i-1];
        }
        // Ghost dxs: valid for linear or log grids
        dx[0] = dx[1]*dx[1]/dx[2] ;
        dx[num_cells+1] = dx[num_cells]*dx[num_cells]/dx[num_cells-1];

        // Surface areas
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
        
        //Compute shell volumes / voluemtric centres
        for(int i=0; i<=num_cells+1; i++) {
            double xl, xr ;
            if (i < num_cells)
                xr = x_i[i];
            else
                xr = x_i[i-1] + dx[i];
            
            if (i > 0)
                xl = x_i[i-1];
            else
                xl = x_i[0] - dx[0] ;

            switch (geometry) {          
            case Geometry::cartesian:
                vol[i] = xr - xl ;
                x_iVC[i] = (xr*xr - xl*xl) / (2*vol[i]) ;
                break;
            case Geometry::cylindrical:
                vol[i] = M_PI * (xr*xr - xl*xl);
                x_iVC[i] = 2 * M_PI * (xr*xr*xr - xl*xl*xl) / (3*vol[i]) ;
                break;
            case Geometry::spherical:
                vol[i] = (4*M_PI/3) * (xr*xr*xr - xl*xl*xl);
                x_iVC[i] = (4*M_PI) * (xr*xr*xr*xr - xl*xl*xl*xl) / (4*vol[i]) ;
                break;
            }
        }
        
        //Compute cell mid positions. Ghost cells also have mid positions in order to balance their pressure gradients
        // but need to be calculated after this loop
        for(int i=1; i<num_cells+1; i++) {
            x_i12[i] = 0.5 * (x_i[i] + x_i[i-1]);
        }
    
        
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
        // !!!!IMPORTANT CHANGE HERE TO NUM_SPECIES_ACT
        //
        num_species = read_parameter_from_file<int>(filename,"PARI_NUM_SPECIES", debug, 1).value;

        if (NUM_SPECIES != Eigen::Dynamic && num_species != NUM_SPECIES) {
            std::stringstream err ;
            err << "Error: You compiled aiolois with a different number of species requested"
                << " to that requested in the parameter file.\n You must recompile with either: "
                << " 1) a dynamic number of species (default) or 2) the correct number of species";
            throw std::invalid_argument(err.str()) ;
        }


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
        //
        // Matrix init via Eigen
        //

        if(debug > 0) cout<<"Init: Setting up friction workspace."<<endl;


        alphas_sample  = Eigen::VectorXd::Zero(num_cells+2);
        
        if(num_species > 1) {
            friction_matrix_T     = Matrix_t::Zero(num_species, num_species);
            friction_coefficients = Matrix_t::Zero(num_species, num_species);
            friction_coeff_mask   = Matrix_t::Ones(num_species, num_species);
            friction_vec_input    = Vector_t(num_species);
            friction_vec_output   = Vector_t(num_species);
            friction_dEkin        = Vector_t(num_species);
            dens_vector           = Vector_t(num_species);
            numdens_vector        = Vector_t(num_species);
            mass_vector           = Vector_t(num_species);

            identity_matrix       = Matrix_t::Identity(num_species, num_species);
            unity_vector          = Vector_t::Ones(num_species);
            LU                    = decltype(LU)(num_species) ;
        }
        else {
            friction_coeff_mask = Matrix_t::Ones(num_species, num_species);
        }
        
        //
        // Punch holes in the coefficient mask for the dust species, so they don't feel each other
        //
        for(int si=0; si < num_species; si++)
            for(int sj=0; sj < num_species; sj++) 
                if(species[si].is_dust_like == 1 && species[sj].is_dust_like == 1) {
                    friction_coeff_mask(si,sj) = 0.;
                    friction_coeff_mask(sj,si) = 0.;
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
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Finished reading boundaries."<<endl;
        if(debug > 0) cout<<"         Boundaries used in species["<<name<<"]: "<<boundary_left<<" / "<<boundary_right<<endl;

        if(base->use_rad_fluxes == 1)
            const_opacity   = read_parameter_from_file<double>(filename,"PARI_CONST_OPAC", debug, 1.).value;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // //Readin species file, const_cv, const_gamma_adiabat for this species
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Reading species data..."<<endl;
        read_species_data(species_filename, species_index); 
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Done reading species data."<<endl;
        
        u               = init_AOS(num_cells+2); //Conserved hyperbolic variables: density, mass flux, energy density
        dudt[0]         = init_AOS(num_cells+2);
        dudt[1]         = init_AOS(num_cells+2);
        source          = init_AOS(num_cells+2);  
        source_pressure = init_AOS(num_cells+2);
        flux            = init_AOS(num_cells+1);

        //
        // Determine which equation of states we are using, and assign thermodynamic variables
        //
        
        //Gas species
        if(is_dust_like == 0) {
            cv            = 0.5 * degrees_of_freedom * Rgas_fake / mass_amu; 
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom;
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu) ;
        }
        //Dust species
        else {
            cv            = 1.5 * Rgas_fake / mass_amu ; // 3 Degrees of freedom per atom (high T limit of debye theory)
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom; // Total degrees of freedom per grain
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu) ;
        }
        
        if(debug >= 0) cout<<"        Species["<<species_index<<"] got a gamma_adiabatic = "<<gamma_adiabat<<" and cv = "<<cv<<endl;

        opacity        = np_somevalue(num_cells+2, const_opacity);
        opticaldepth   = np_zeros(num_cells+2);
        radiative_flux = np_zeros(num_cells+1);
        u_analytic     = np_zeros(num_cells+2);

        prim   = std::vector<AOS_prim>(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        prim_l = std::vector<AOS_prim>(num_cells+2);
        prim_r = std::vector<AOS_prim>(num_cells+2);
        
        timesteps    = np_zeros(num_cells+2);		    
        timesteps_cs = np_zeros(num_cells+2);	
        finalstep    = np_zeros(num_cells+2);
        
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Initialize/readin scenario parameters
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
            double u3r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3R", debug, 1.).value;
            
            SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug, 0.5).value;
            
            u1l = initial_fraction * u1l;
            u1r = initial_fraction * u1r;
            
            if (is_dust_like) { 
                u3l *= 2 / mass_amu ;
                u3r *= 2 / mass_amu ;
            } 
            
            if(base->init_wind) {
                
                u2l = density_excess;
                u2r = density_excess;
            }

            // Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            AOS_prim pl(u1l, u2l, u3l) ;
            eos->compute_conserved(&pl, &SHOCK_TUBE_UL, 1) ;

            AOS_prim pr(u1r, u2r, u3r) ;
            eos->compute_conserved(&pr, &SHOCK_TUBE_UR, 1) ;
            
            initialize_shock_tube_test(SHOCK_TUBE_UL, SHOCK_TUBE_UR);
            
            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 1. with left primitives ="<<pl.density<<" /"<<pl.speed<<" init atmo "<<init_static_atmosphere  <<" density excess"<<density_excess<<endl;
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
            AOS_prim p(u1 * initial_fraction, u2, u3);
            eos->compute_conserved(&p, &BACKGROUND_U, 1);
            
            initialize_background(BACKGROUND_U);
            
            //
            // If we want to initialize a hydrostatic atmosphere, then we overwrite the so far given density/pressure data and set every velocity to 0
            //
            if(init_static_atmosphere == 1) {
                initialize_hydrostatic_atmosphere();
            }

            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 2."<<endl;
        }
        else if (base->problem_number == 3) {
            initialize_sound_wave() ;
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
    
    double temp_rhofinal;
    double factor_inner, factor_outer;
    double T_inner;
    double T_outer;
    double metric_inner;
    //double residual;
    double metric_outer;
    //
    // Start with the outermost cell and build up a hydrostatic atmosphere
    // Fulfilling KKM16, Eqn. 17
    //
    //long double residual  = 0.;
    
    //
    // First, initialize (adiabatic) temperature
    //
    for(int i=num_cells+1; i>=0; i--) {
            
            prim[i].temperature = - 1.0 * base->phi[i] / (cv * gamma_adiabat) + const_T_space;
            //prim[i].temperature = const_T_space;
            
            //Add temperature bumps and troughs
            prim[i].temperature += TEMPERATURE_BUMP_STRENGTH * 4.  * exp( - pow(base->x_i12[i] - 1.e-1 ,2.) / (0.1) );
            prim[i].temperature -= TEMPERATURE_BUMP_STRENGTH * 15. * exp( - pow(base->x_i12[i] - 7.e-3 ,2.) / (1.5e-3) );

    }
    
    //At this point, the right ghost cell is already initialized, so we can just build up a hydrostatic state from there
    for(int i=num_cells; i>=0; i--)  {
        
        //
        // Construct next density as to fulfil the hydrostatic condition
        //
        T_outer = prim[i+1].temperature ;
        T_inner = prim[i].temperature ;

        //factor_outer = (gamma_adiabat-1.) * cv * T_outer; 
        //factor_inner = (gamma_adiabat-1.) * cv * T_inner; 
        
        eos->get_p_over_rho_analytic(&T_outer, &factor_outer);
        eos->get_p_over_rho_analytic(&T_inner, &factor_inner);
        
        metric_outer = (base->phi[i+1] - base->phi[i]) * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
        metric_inner = (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
        
        temp_rhofinal = u[i+1].u1 * (factor_outer + metric_outer) / (factor_inner - metric_inner);
        
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner) ;
        
        //
        // Debug info
        // 
        //if( (i==20 || i== num_cells-20) && debug >= 0) {
        //if( i==300 ) {
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
            cout<<"     Ratio of densities inner/outer = "<< temp_rhofinal/u[i+1].u1 <<endl;
            cout<<"     Ratio of temperatures inner/outer = "<<T_inner/T_outer<<" t_inner ="<<T_inner<<" t_outer ="<<T_outer<<endl;
            cout<<"     Ratio of pressures inner/outer = "<<cv * temp_rhofinal * T_inner /u[i+1].u3<<endl;
            cout<<"     Resulting density == "<<temp_rhofinal<<endl;
            cout<<"     density before "<<u[i+1].u1<<endl;
            cin>>a;
        }
    }
    
    if(u[2].u1 > 1e40) {
        cout<<"    ---WARNING---"<<endl;
        cout<<"    IN CONSTRUCT HYDROSTATIC:"<<endl;
        cout<<"    species["<<name<<"] has reached a very high initial density  of > 1e40 at the inner boundary."<<endl;
        cout<<"    The well-balanced scheme seems to be unable to compensate minor velocity fluctuations under those circumstances."<<endl;
        cout<<"    Be advised, your solution is probably about to unphysically explode."<<endl;
    }
        
  
    compute_pressure(u);
    
    //
    // Apply multiplicators for over/underdensities just below and above the sonic point
    //
    if(base->init_wind==1) {
        double csi              = prim[num_cells-1].sound_speed;
        double smoothed_factor;
        bondi_radius = base->planet_mass / (pow(csi,2.)); //Most restrictive sonic point, as we are not quite isothermal
    
        for( int i = 1; i < num_cells + 1; i++) {
            smoothed_factor = density_excess * (1./(1.+std::exp(-(base->x_i12[i] - bondi_radius)/(0.01 * bondi_radius )) ) ); 
            
            u[i].u1 = u[i].u1 * (1. + smoothed_factor);
            
            u[i].u3 = u[i].u3 * (1. + smoothed_factor);
            if(i==1)
                cout<<"In initwind: First smoothed factor = "<<smoothed_factor<<endl;
        }
        
        cout<<"            Ended hydrostatic construction for species "<<name<<", last smoothed factor was "<<smoothed_factor<<" while  density_excess="<<density_excess<<" r_b="<<bondi_radius<<endl;
    }
    else
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

//
void c_Species::initialize_sound_wave() {

    double amp = 1e-6 ;
    double rho0 = 1.0 * initial_fraction ;
    double cs = 1.0 * (base->domain_max - base->domain_min) ;
    double p0 = rho0*cs*cs/(gamma_adiabat) ;

    // Assume gas is H2
    if (is_dust_like)
        p0 *= 2 / mass_amu ;

    for(int i=0; i<=num_cells+1; i++) {
        double kx = 2*M_PI * (i+0.5 - base->num_ghosts) / double(num_cells - 2*(base->num_ghosts-1)) ;
        double f = amp * std::sin(kx) ;

        AOS_prim p (rho0 *(1 + f), cs*f, p0*(1 + gamma_adiabat*f)) ;
        eos->compute_conserved(&p, &u[i], 1) ;
    }
}


void c_Species::apply_boundary_left(std::vector<AOS>& u) {
    int num_ghosts = base->num_ghosts;
    int Ncell = num_cells - 2*(num_ghosts-1) ; // Correct for fact we increased num_cells
    switch(boundary_left) {
        case BoundaryType::user:
            user_boundary_left(u);
            break;
        case BoundaryType::open:
            for (int i=num_ghosts; i > 0; i--) {
                AOS_prim prim ;
                eos->compute_primitive(&u[i],&prim, 1) ;
                double dphi = (base->phi[i] - base->phi[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                dphi       *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                prim.pres   = prim.pres + prim.density * dphi ;
                prim.pres   = std::max( prim.pres, 0.0) ;
                eos->compute_conserved(&prim, &u[i-1], 1) ;
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
                int iact = Ncell + i;
                u[i] = u[iact];
                base->phi[i]   = base->phi[iact] ;
            }
            break;
    }
}

void c_Species::apply_boundary_right(std::vector<AOS>& u) {
    int num_ghosts = base->num_ghosts;
    int Ncell = num_cells - 2*(num_ghosts-1) ;
    switch(boundary_right) {
        case BoundaryType::user:
            user_boundary_right(u);
            break;
        case BoundaryType::open:
            for (int i=Ncell+num_ghosts; i < Ncell+2*num_ghosts; i++) {
                AOS_prim prim ;
                eos->compute_primitive(&u[i-1],&prim, 1) ;
                double dphi = (base->phi[i] - base->phi[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                dphi *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                prim.pres = prim.pres -  prim.density * dphi ;
                prim.pres    = std::max( prim.pres, 0.0) ;
                eos->compute_conserved(&prim, &u[i], 1) ;
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
                int iact = i - Ncell;
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
    // Delete EOSs, note here, not c_Species.
    for (int i=0; i < num_species; i++) 
        delete species[i].eos ;
}

c_Species::~c_Species() {
    //Empty
}

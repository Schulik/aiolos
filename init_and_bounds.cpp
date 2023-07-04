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

/**
 * Simulation class constructor and initializer of all the physics. Reads all simulation parameters from the *.par file, constructs all species using data from the *spc and linked
 * *opa files.
 * 
 * Computes the spatial grid, areas, volumes, wavelength grid. Initiates matrix solver objects,  species data (around line 580) and radiative quantities.
 * 
 * @param[in] filename_solo Parameter file without any other appendages, directory.
 * @param[in] speciesfile_solo Species file without any other appendages, directory.
 * @param[in] workingdir Working directory string
 * @param[in] debug debug level from command line or param file
 * @param[in] debug_cell Get detailed information for a specific cell (To be implemented by user..)
 * @param[in] debug_steps Get detailed information for a specific timestep (To be implemented by user..)
 */
c_Sim::c_Sim(string filename_solo, string speciesfile_solo, string workingdir, int debug, int debug_cell, int debug_steps) {

        if(debug > 0) cout<<"Init position 0."<<endl;
        
        steps = -1;
        this->debug      = debug ;
        this->debug_cell = debug_cell;
        this->debug_steps= debug_steps;
        simname          = filename_solo;
        this->workingdir = workingdir;
        string filename    = workingdir + filename_solo;
        parfile            = workingdir + filename_solo;
        
        
        if(speciesfile_solo.compare("default.spc")==0)
            speciesfile_solo = read_parameter_from_file<string>(filename,"SPECIES_FILE", debug, "default.spc").value; //The program always needs to know some information about the species used.
        
        cout<<" Starting construction of simulation. Attempting to find required speciesfile = "<<speciesfile_solo<<endl;
        
        speciesfile = workingdir + speciesfile_solo;
        
        cout<<" Opening parameter file "<<filename<<endl;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //  Numerical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double dx0;
        type_of_grid     = read_parameter_from_file<int>(filename,"PARI_GRID_TYPE", debug, 0).value; //Spacing type of grid. Linear, log, bi-log etc. More in enum.h
        domain_min       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MIN", debug).value; //Smallest value of 1-D domain in cm 
        domain_max       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MAX", debug).value; //Largest value of 1-D domain   in cm
        geometry         = read_parameter_from_file<Geometry>(filename, "PARI_GEOMETRY", debug, Geometry::cartesian).value; //Cartesian, Polar, Spherical. Determines differentials. More in enum.h
        order            = read_parameter_from_file<IntegrationType>(filename, "PARI_ORDER", debug, IntegrationType::second_order).value; //Spatial order of differentials.
        
        num_bands_in     = read_parameter_from_file<int>(filename,"PARI_NUM_BANDS", debug, 1).value; //Number of instellation bands. Low performance impact.
        num_bands_out    = read_parameter_from_file<int>(filename,"NUM_BANDS_OUT", debug, num_bands_in).value; //Number of outgoing radiation bands. Enormous performance impact.
        lambda_min_in       = read_parameter_from_file<double>(filename,"PARI_LAM_MIN", debug, 1e-1).value; //Minimum ingoing, non-zero limit of the wavelength domain in microns.
        lambda_max_in       = read_parameter_from_file<double>(filename,"PARI_LAM_MAX", debug, 10.).value;  //Maximum ingoing, non-inf limit of the wavelength domain in microns. 
                             //Note how by default there is another band between 0 and lambda_min_in, and another band lamda_max_in and infty, if the band numbers are high enough
        lambda_per_decade_in= read_parameter_from_file<double>(filename,"PARI_LAM_PER_DECADE", debug, 10.).value; //Log-equidistant spacing in wavelngth
        lambda_min_out       = read_parameter_from_file<double>(filename,"PARI_LAM_MIN_OUT", debug, lambda_min_in).value;
        lambda_max_out       = read_parameter_from_file<double>(filename,"PARI_LAM_MAX_OUT", debug, lambda_max_in).value;
                            //Note how by default there is another band between 0 and lambda_min_out, and another band lamda_max_out and infty, if the band numbers are high enough
        lambda_per_decade_out= read_parameter_from_file<double>(filename,"PARI_LAM_PER_DECADE_OUT", debug, lambda_per_decade_in).value;
        fluxfile             = read_parameter_from_file<string>(filename,"FLUX_FILE", debug, "---").value;  //File to read in the wavelength-dependent data of the
        fluxmultiplier       = read_parameter_from_file<double>(filename,"FLUX_MULTIPLIER", debug, 1.).value; //Wavelength-constant multiplier for all top-of-atmosphere band fluxes
        
        star_mass        =  read_parameter_from_file<double>(filename,"PARI_MSTAR", debug, 1.).value; //Mass of the star in solar masses
        //star_mass        *= msolar;
        init_star_mass   = msolar * read_parameter_from_file<double>(filename,"INIT_MSTAR", debug, star_mass).value;
	star_mass        *= msolar;
        ramp_star_mass_t0 = read_parameter_from_file<double>(filename,"RAMP_MSTAR_T0", debug, 0.).value;
        ramp_star_mass_t1 = read_parameter_from_file<double>(filename,"RAMP_MSTAR_T1", debug, 1.).value;

        T_star           = read_parameter_from_file<double>(filename,"PARI_TSTAR", debug, 5777.).value; //tau=2/3 temperature of the star, i.e. Teff in K
        R_star           = read_parameter_from_file<double>(filename,"PARI_RSTAR", debug, 0.).value;    //radius of the star in solar radii
        UV_star          = read_parameter_from_file<double>(filename,"PARI_UVSTAR", debug, 0.).value;   // UV luminosity of the star, used if no flux file is given, in erg/s
        X_star           = read_parameter_from_file<double>(filename,"PARI_XSTAR", debug, 0.).value;    //X-ray luminosity of the star, not used currently
        Lyalpha_star     = read_parameter_from_file<double>(filename,"PARI_LYASTAR", debug, 0.).value;  // Not used currently
        planet_semimajor= read_parameter_from_file<double>(filename,"PARI_PLANET_DIST", debug, 1.).value; //Distance to primary black-body in AU
        
        R_other          = read_parameter_from_file<double>(filename,"R_OTHER", debug, 0.).value;  //Second bolometric black-body source (i.e. second star or giant planet host) in stellar radii
        T_other          = read_parameter_from_file<double>(filename,"T_OTHER", debug, 0.).value;  //Temperature of other black-body in K
        d_other          = read_parameter_from_file<double>(filename,"D_OTHER", debug, 1.).value;  // Distance to other black-body in AU
       
        T_int            = read_parameter_from_file<double>(filename,"PARI_TINT", debug, 0.).value; //Internal planetary temperature at inner boundary. T_int in Guillot2010. in K.
        // old notation.
        if (T_int==0) 
            T_int        = read_parameter_from_file<double>(filename,"PARI_TPLANET", debug, 0.).value; //Internal planetary temperature at inner boundary. T_int in Guillot2010. in K.

        use_planetary_temperature = read_parameter_from_file<int>(filename,"USE_PLANET_TEMPERATURE", debug, 0).value;//Add the T_int heating to the lower boundary cell?
        core_cv          = read_parameter_from_file<double>(filename,"PARI_CORE_CV", debug, 0).value; //Heat capacity of the core 
        T_planet         = read_parameter_from_file<double>(filename,"PARI_TPLANET_INIT", debug, 0).value; //To be unused?
        //Note that using T_planet might require smaller timesteps (if they are not already restricted by CFL), in those cases use max_timestep_change to limit the dt increase
        
        radiation_matter_equilibrium_test = read_parameter_from_file<int>(filename,"RAD_MATTER_EQUI_TEST", debug, 0).value; //Unused/only for tests
        radiation_diffusion_test_linear   = read_parameter_from_file<int>(filename,"RAD_DIFF_TEST_LIN", debug, 0).value;    //Unused/only for tests
        radiation_diffusion_test_nonlinear= read_parameter_from_file<int>(filename,"RAD_DIFF_TEST_NLIN", debug, 0).value;   //Unused/only for tests
        couple_J_into_T = read_parameter_from_file<int>(filename,"COUPLE_J_INTO_T", debug, 1).value;   //Unused/only for tests
        
        if(debug > 0) cout<<"Using integration order "<<order<<" while second order would be "<<IntegrationType::second_order<<endl;
        
        if (order == IntegrationType::first_order)
            num_ghosts = 1 ;
        else 
            num_ghosts = 2 ;

        if(type_of_grid == 0) {
            dx0              = read_parameter_from_file<double>(filename,"PARI_DOMAIN_DX", debug).value; //For uniform grid, resolution element size
            num_cells        = (int)((domain_max - domain_min)/dx0);
            cout<<"Domain specifics:  "<<domain_min<<" | . . . "<<num_cells<<" uniform cells . . . | "<<domain_max<<endl;
        }
        else {
            cells_per_decade = read_parameter_from_file<double>(filename,"PARI_CELLS_PER_DECADE", debug).value; //For logarithmic grid, number cells per decade in resolution
            num_cells        = (int) ( (log10f(domain_max) - log10f(domain_min)) * cells_per_decade );
            dx0              = domain_min;
            
            if(type_of_grid==2) {
                    grid2_transition       = read_parameter_from_file<double>(filename,"GRID2_TRANSITION", debug, 99.e99).value; //For bi-logarithmic grid (two different resolution regions), define transition boundary
                    grid2_cells_per_decade = read_parameter_from_file<double>(filename,"GRID2_CELLS_PER_DECADE", debug, 10).value; // For bi-logarithmic grid, number cells per decade in resolution in second region at large r (first region sits at small r)
                    
                    num_cells          = (int) ( (log10f(grid2_transition) - log10f(domain_min)) * cells_per_decade );
                    grid2_transition_i = num_cells + num_ghosts;
                    num_cells         += (int) ( (log10f(domain_max) - log10f(grid2_transition)) * grid2_cells_per_decade );
                    
            }
            cout<<"Domain specifics:  "<<domain_min<<" | . .  .  "<<num_cells<<" nonuniform cells .       .             .     | "<<domain_max<<endl;
        }
        reverse_hydrostat_constrution = read_parameter_from_file<int>(filename,"REVERSE_HYDROSTAT_CONSTRUCTION", debug, 0).value; //If false, or zero, construct density from large r to small r.
                                                                                                                   //Otherwise reversed. First density is the PARI_INIT_DATA_U1parameter.
        init_wind         = read_parameter_from_file<int>(filename,"PARI_INIT_WIND", debug, 0).value; //Initialize a step in density at init_sonic_radius of variable magnitude. 
        //Magnitude of density jump is set by species dependent parameter in *spc file in column 7 "initial density excess". If init_sonic_radius is negative, then do nothing. 
        init_sonic_radius = read_parameter_from_file<double>(filename,"INIT_SONIC_RADIUS", debug, -1e10).value;
        
        //
        // Check that cell number looks fine. Shoutout if its not.
        //
        if( !(num_cells > 0 && num_cells < 999999) ) {
            std::stringstream err ;
            err<<"    WARNING! Something seems wrong with the number of cells required. num_cells = "<<num_cells<<endl;
            throw std::invalid_argument(err.str()) ;
        }
            
        if(debug > 0) cout<<"Init: Finished reading grid parameters."<<endl;
        
        cflfactor   = read_parameter_from_file<double>(filename,"PARI_CFLFACTOR", debug, 0.9).value;  //Multiplier on the cfl timestep length
        t_max       = read_parameter_from_file<double>(filename,"PARI_TIME_TMAX", debug, 1e0).value;  // Simulate until t=t_max in seconds.
        dt_max      = read_parameter_from_file<double>(filename,"PARI_DTMAX", debug, 1e99).value;     // Limit the largest possible timestep size in seconds.
        max_timestep_change = read_parameter_from_file<double>(filename,"MAX_TIMESTEP_CHANGE", debug, 1.1).value; //Max. change of timestep per following timestep in s.
        dt_min_init         = read_parameter_from_file<double>(filename,"DT_MIN_INIT", debug, 1e-20).value;       //Initial dt in s
        output_time = read_parameter_from_file<double>(filename,"PARI_TIME_OUTPUT", debug, 1e99).value;           //Create an output every xxx simulated seconds.
        output_time_offset = read_parameter_from_file<double>(filename,"TIME_OUTPUT_OFFSET", debug, 0.).value;    //Create outputs every PARI_TIME_OUTPUT but only starting after offset, in s
        monitor_time = read_parameter_from_file<double>(filename,"PARI_TIME_DT", debug).value;                    //Put measurements into the monitor file every xx s
        CFL_break_time = read_parameter_from_file<double>(filename,"CFL_BREAK_TIME", debug, std::numeric_limits<double>::max()).value ; //Use PARI_CFLFACTOR if t<CLF_break_time. Otherwise, set cflfactor to 0.9
        energy_epsilon = read_parameter_from_file<double>(filename,"ENERGY_EPSILON", debug, 0.01).value; //Relative allowed change of internal energy in any single grid cell. Limits timestep size additional to the CFL condition.
        
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
        
        problem_number    = read_parameter_from_file<int>(filename,"PARI_PROBLEM_NUMBER", debug).value;     //0==gas tube, 1==planet run
        use_self_gravity  = read_parameter_from_file<int>(filename,"PARI_SELF_GRAV_SWITCH", debug, 0).value;//Gravitational field of gas in domain added to potential?
        use_tides         = read_parameter_from_file<int>(filename,"USE_TIDES", debug, 0).value;            //Use the tidal field from the first central object (i.e. star_mass, R_star etc)?
        use_linear_gravity= read_parameter_from_file<int>(filename,"PARI_LINEAR_GRAV", debug, 0).value;     //Linear acceleration ~r instead of 1/r^2
        use_rad_fluxes    = read_parameter_from_file<int>(filename,"PARI_USE_RADIATION", debug, 0).value;   //Switch to turn on the radiation module
        use_convective_fluxes = read_parameter_from_file<int>(filename,"USE_CONVECTION", debug, 0).value;   //Switch to turn on convective energy transport in the radiation module
        K_zz_init = read_parameter_from_file<double>(filename,"KZZ_INIT", debug, 0.).value;                 //Initial atmospheric mixing parameter in cm^2/s
        convect_boundary_strength = read_parameter_from_file<double>(filename,"CONVECT_BOUNDARY_STRENGTH", debug, 1.1).value; //Unused currently.
        
        use_collisional_heating = read_parameter_from_file<int>(filename,"PARI_USE_COLL_HEAT", debug, 1).value; //Switch on collisional energy exchange between species
        use_drag_predictor_step = read_parameter_from_file<int>(filename, "PARI_SECONDORDER_DRAG", debug, 0).value; //Switch on drag predictor substep
        alpha_collision        = read_parameter_from_file<double>(filename,"PARI_ALPHA_COLL", debug, 1.).value; //Multiplier for collision alphas 
        alpha_collision_ions   = read_parameter_from_file<double>(filename,"PARI_ALPHA_IONS", debug, 1.).value; //Multiplier for collision alphas for ions
        max_mdot              = read_parameter_from_file<double>(filename,"MAX_MDOT", debug, -1.).value;        //Max negative mdot through outer boundary in g/s
        
        rad_energy_multiplier=read_parameter_from_file<double>(filename,"PARI_RAD_MULTIPL", debug, 1.).value; //Only for tests. Unused currently.
        collision_model   = read_parameter_from_file<char>(filename,"PARI_COLL_MODEL", debug, 'P').value;     //Collision model. P=physical collision rates (i.e. Schunk&Nagy). C=constant collision rate of value PARI_ALPHAS_COLL
        opacity_model     = read_parameter_from_file<char>(filename,"PARI_OPACITY_MODEL", debug, 'C').value; //Opacity model. List in opacities.cpp line 90.
        use_avg_temperature = read_parameter_from_file<int>(filename,"USE_AVG_T", debug, 0).value; //Compute avg c_v rho T  and use that as T_avg
	avg_temperature_t0  = read_parameter_from_file<double>(filename,"avg_temperature_t0", debug, 1e99).value; //All T_s = T_avg before t0. Between t0 and t1 a ramp-down begins
	avg_temperature_t1  = read_parameter_from_file<double>(filename,"avg_temperature_t1", debug, 1e99).value; //All T_s remain T_s after t1. T_avg is ignored after t1

        cout<<" CHOSEN OPACITY MODEL = "<<opacity_model<<endl;
        if(opacity_model == 'M')
            init_malygin_opacities();
        if(opacity_model == 'K') {
                cout<<" Hybrid opacities selected. NOTE! that the input file expects in this mode two file names to be given for each species, one *aiopa for kR and kP, and *opa for kS."<<endl;
        }
        
        no_rad_trans               = read_parameter_from_file<double>(filename,"NO_RAD_TRANS", debug, 1.).value; //Multiplier for strength for thermal radiative losses in radiation transport. Set to 1e-100 to emulate perfect energy-limited escape.
        solve_for_j                = read_parameter_from_file<int>(filename,"SOLVE_FOR_J", debug, 1).value; //Debugging parameter. Switch to zero for decoupling of T and J in simple rad transport
        photocooling_multiplier    = read_parameter_from_file<double>(filename,"PHOTOCOOL_MULTIPLIER", debug, 1.).value; //Multiplier for non-thermal cooling rates
        photocooling_expansion     = read_parameter_from_file<double>(filename,"PHOTOCOOL_EXPANSION", debug, 1.).value; //Multiplier for non-thermal second order cooling rates
        radiation_rampup_time      = read_parameter_from_file<double>(filename,"RAD_RAMPUP_TIME", debug, 0.).value; //Ramp up the irradiation in all bands smoothly over xxx seconds.
        init_radiation_factor      = read_parameter_from_file<double>(filename,"INIT_RAD_FACTOR", debug, 0.).value; //Unused
        //radiation_solver           = read_parameter_from_file<int>(filename,"RADIATION_SOLVER", debug, 0).value; //replaced by use_rad_fluxes
        closed_radiative_boundaries = read_parameter_from_file<int>(filename,"PARI_CLOSED_RADIATIVE_BOUND", debug, 0).value; //Reflect thermal radiation at outer boundaries?
        const_opacity_solar_factor = read_parameter_from_file<double>(filename,"CONSTOPA_SOLAR_FACTOR", debug, 1.).value;    //Multiplier for all stellar opacities (independent of opacity model)
        const_opacity_rosseland_factor = read_parameter_from_file<double>(filename,"CONSTOPA_ROSS_FACTOR", debug, 1.).value; //Multiplier for all rosseland opacities (")
        const_opacity_planck_factor = read_parameter_from_file<double>(filename,"CONSTOPA_PLANCK_FACTOR", debug, 1.).value;  //Multiplier for all planck opacities (")
        temperature_model = read_parameter_from_file<char>(filename,"INIT_TEMPERATURE_MODEL", debug, 'P').value; //Initialize temperature as adiabatic+constant (P) or constant (C)
        friction_solver   = read_parameter_from_file<int>(filename,"FRICTION_SOLVER", debug, 0).value; //0: no friction, 1: analytic (only for two species), 2: numerical
        update_coll_frequently = read_parameter_from_file<int>(filename,"UPDATE_COLL_FREQUENTLY", debug, 1).value; //Switches on or off another update of collision rates after temperature and before temperature exchange update for rad solver==2. Important for code performance with many species
        do_hydrodynamics  = read_parameter_from_file<int>(filename,"DO_HYDRO", debug, 1).value;  //Switch hydrodynamics, including its CFL condition on/off
        start_hydro_time  = read_parameter_from_file<double>(filename,"START_HYDRO_TIME", debug, -1.).value; //Hydro is started once globalTime > start_hydro_time and do_hydro == 1
        photochemistry_level = read_parameter_from_file<int>(filename,"PHOTOCHEM_LEVEL", debug, 0).value; //1==C2Ray solver, 2==general photo and thermochemistry solver
        dust_to_gas_ratio = read_parameter_from_file<double>(filename,"DUST_TO_GAS", debug, 0.).value;    //dust-to-gas ratio for Semenov/Malygin opacities
        temperature_floor = read_parameter_from_file<double>(filename,"TEMPERATURE_FLOOR", debug, 0.).value;  //Limit the temperature to a minimum in radiation solver (negative T crashes can still occur due to negative pressures, whichever is found first)
        max_temperature   = read_parameter_from_file<double>(filename,"TEMPERATURE_MAX", debug, 9e99).value;  //Limit the temperature to a maximum
        use_chemistry      = read_parameter_from_file<int>(filename,"DO_CHEM", debug, 0).value;               //Unused currently.
        use_total_pressure = read_parameter_from_file<int>(filename,"USE_TOTAL_PRESSURE", debug, 0).value;    //Compute total pressure of all species? Unused currently.
        coll_rampup_time      = read_parameter_from_file<double>(filename,"COLL_RAMPUP_TIME", debug, 0.).value; //Ramp up collision rates gently from INIT_COLL_FACTOR to PARI_ALPHA_COLL in the rampup time
        init_coll_factor      = read_parameter_from_file<double>(filename,"INIT_COLL_FACTOR", debug, 0.).value;
        
        chemistry_precision         = read_parameter_from_file<double>(filename,"CHEM_PRECISION", debug, 1e-2).value; //Relative preicion for the number densities to be accepted by solver
        chemistry_numberdens_floor         = read_parameter_from_file<double>(filename,"CHEM_FLOOR", debug, 1e-20).value; //Lower number densities floor per cell, normalized to total cell dens
        chemistry_maxiter           = read_parameter_from_file<int>(filename,"CHEM_MAXITER", debug, 4).value; //Number of max chem solver iterations per timestep. Increase as 2^int
        chemistry_miniter           = read_parameter_from_file<int>(filename,"CHEM_MINITER", debug, 4).value; //Number of min chem solver iterations per timestep. 
        ion_precision         = read_parameter_from_file<double>(filename,"ION_PRECISION", debug, 1e-12).value; //Precision of brent solver for C2Ray number densities
        ion_heating_precision = read_parameter_from_file<double>(filename,"ION_HEATING_PRECISION", debug, 1e-12).value; //Precision of brent solver for C2Ray heating rates
        ion_maxiter                  = read_parameter_from_file<int>(filename,"ION_MAXITER", debug, 128).value;         //Maxiter of brent solver for C2Ray number densities
        ion_heating_maxiter          = read_parameter_from_file<int>(filename,"ION_HEATING_MAXITER", debug, 128).value;//Maxiter of brent solver for C2Ray heating rates
        intermediate_chemfloor_check = read_parameter_from_file<int>(filename,"INTERMEDIATE_CHEM_CHECK", debug, 0).value; //Check chem_floor in between iterations?
        chem_momentum_correction     = read_parameter_from_file<int>(filename,"CHEM_MOMENTUM_CORR", debug, 1).value;      //Compute momentum transport due to ndot
        chem_ekin_correction         = read_parameter_from_file<int>(filename,"CHEM_EKIN_CORR", debug, 1).value;          //Compute energy transport due to ndot
        dt_skip_ichem                = read_parameter_from_file<int>(filename,"CHEM_DT_SKIP", debug, 1).value;            //Skip chemistry update every xxx timesteps in main loop. Accumulate timesteps until next chem solver call. Currently buggy.
        
        right_extrap_press_multiplier =read_parameter_from_file<double>(filename,"BOUND_EXTRAP_MUL", debug, 1.0).value; //Hydrodynamic pressure extrapolation into ghost cell multiplier
        dt_skip_dchem = 0.; 
        
        rad_solver_max_iter = read_parameter_from_file<int>(filename,"MAX_RAD_ITER", debug, 1).value;  //Currently unused
        xi_rad              = read_parameter_from_file<double>(filename,"XI_RAD", debug, 2.).value;    //Radiation Xi factor, see Fig. 5 in code paper.
        bond_albedo       = read_parameter_from_file<double>(filename,"BOND_ALBEDO", debug, 0.).value; //Bond albedo, Number between 0. and 1.
        use_init_discont_smoothing = read_parameter_from_file<int>(filename,"INIT_DISCONT_SMOOTHING", debug, 0).value;  //Smooths out discontinuities in the initial density profiles with a powerlaw. Applied before init_wind discontinuity is put in.
        
        if(problem_number == 2)
            monitor_output_index = num_cells/2; 
        else
            monitor_output_index = 1;
        
        if(problem_number==2) {
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug, 1).value; //in Earth masses
            planet_mass     *= mearth;
            
            planet_position = read_parameter_from_file<double>(filename,"PARI_PLANET_POS", debug, 0.).value;  //In spatial coordinates. 
            rs              = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_LENGTH", debug, 0.).value; //Gravitational smoothing length in rhill 
            rs_time         = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_TIME", debug, 0.).value; //Ramp-up time for smoothing length from 0 till PARI_SMOOTHING_LENGTH
            rs_at_moment    = 0.2;
            
        } else {
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug, 0).value; //in Earth masses
            planet_mass     *= mearth;
            
            planet_position = read_parameter_from_file<double>(filename,"PARI_PLANET_POS", debug, 0.).value;  //inside the simulation domain
            rs              = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_LENGTH", debug, 0.).value; //Gravitational smoothing length in hill 
            rs_time         = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_TIME", debug, 0.).value; //Time until we reach rs starting at rs_at_moment
            rs_at_moment    = 0.2;
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
        enclosed_mass     = np_zeros(num_cells+2);
        enclosed_mass_tmp = np_zeros(num_cells+2);
        phi               = np_zeros(num_cells+2);
        //total_pressure    = np_zeros(num_cells+2);
        total_press_l            = std::vector<double>(num_cells+2);
        total_press_r            = std::vector<double>(num_cells+2);
        total_adiabatic_index   = std::vector<double>(num_cells+2);
        
        if(num_bands_in < 1) {
            cout<<" In INIT RADIATION, invalid num_bands_in = "<<num_bands_in<<" changing to num_bands_in = 1."<<endl;
            num_bands_in = 1;
        }
        if(num_bands_out < 1) {
            cout<<" In INIT RADIATION, invalid num_bands_out = "<<num_bands_out<<" changing to num_bands_out = 1."<<endl;
            num_bands_out = 1;
        }
        if(use_rad_fluxes == 2 && use_convective_fluxes == 1 ) {
            cout<<" Simple radiative solver used, but with convection. Note that convection in this case is not implemented yet!"<<endl;
        }
        
        l_i_in        = np_zeros(num_bands_in+1);    // Wavelenght bin boundaries
        l_i12_in      = np_zeros(num_bands_in);      // Wavelength bin midpoints or averages
        l_i_out       = np_zeros(num_bands_out+1);    // Wavelenght bin boundaries
        l_i12_out     = np_zeros(num_bands_out);      // Wavelength bin midpoints or averages
        BAND_IS_HIGHENERGY = inp_zeros(num_bands_in); //We save wether a band is a highenergy band. This is important for which solver is activated.
        
        l_i_in[0]              = lminglobal/100.;     //Transform wavelengths into microns
        l_i_in[num_bands_in]   = lmaxglobal/100.;
        l_i_out[0]             = lminglobal/100.;     //Transform wavelengths into microns
        l_i_out[num_bands_out] = lmaxglobal/100.;
        
        /////////////////////////////////////////////////////////
        ////////////////////////////////Bands in
        /////////////////////////////////////////////////////////
        if(num_bands_in == 2) {
            l_i_in[1]   = 0.5*(lambda_max_in + lambda_min_in);
            l_i12_in[0] = lambda_min_in;
            l_i12_in[1] = lambda_max_in;
        }
        else if(num_bands_in > 2)
        {
            l_i_in[1]                = lambda_min_in;
            l_i_in[num_bands_in-1]   = lambda_max_in;
            double dlogl2      = pow(lambda_max_in/lambda_min_in, 1./(num_bands_in-2));
            
            for(int b=2; b<num_bands_in-1; b++) {
                l_i_in[b]      = l_i_in[b-1] * dlogl2;
                cout<<" in NUM_BANDS>2, b = "<<b<<" l_i_in[b] = "<<l_i_in[b]<<endl;
            }
            
            for(int b=0; b<num_bands_in; b++) {
                l_i12_in[b]  = pow( 10., 0.5 * (std::log10(l_i_in[b]) + std::log10(l_i_in[b+1])));   
                
            }
        }
        /////////////////////////////////////////////////////////
        ////////////////////////////////Bands out
        /////////////////////////////////////////////////////////
        if(num_bands_out == 2) {
            l_i_out[1]   = 0.5*(lambda_max_out + lambda_min_out);
            l_i12_out[0] = lambda_min_out;
            l_i12_out[1] = lambda_max_out;
        }
        else if(num_bands_out > 2)
        {
            l_i_out[1]             = lambda_min_out;
            l_i_out[num_bands_out-1]   = lambda_max_out;
            double dlogl2      = pow(lambda_max_out/lambda_min_out, 1./(num_bands_out-2));
            //l_i[0] = lambda_min;
            
            for(int b=2; b<num_bands_out-1; b++) {
                l_i_out[b]      = l_i_out[b-1] * dlogl2;
                cout<<" in NUM_BANDS>2, b = "<<b<<" l_i_out[b] = "<<l_i_out[b]<<endl;
            }
            
            for(int b=0; b<num_bands_out; b++) {
                l_i12_out[b]  = pow( 10., 0.5 * (std::log10(l_i_out[b]) + std::log10(l_i_out[b+1])));   
                
            }
        }
        
        //
        // Assign high and low energy band switches. Those will later decide about the computation of dS and ionization
        // Low-energy bands directly input solar radiation into the thermal energy. High-energy bands ionize and are therefore more complicated to treat.
        //
        num_he_bands = 0;
        for(int b=0; b<num_bands_in; b++) {
            BAND_IS_HIGHENERGY[b] = 0;
            
            if(photochemistry_level > 0)
                if (l_i_in[b + 1] < 0.09116) {
                    BAND_IS_HIGHENERGY[b] = 1;
                    num_he_bands++;
                }
        }
        photon_energies = np_somevalue(num_bands_in, 20. * ev_to_K * kb);
        
        for(int b=0; b<num_bands_in; b++) {
            photon_energies[b] = 1.24/( l_i_in[b + 1] ) * ev_to_K * kb ;
            cout<<"Assigned to highenergy band b = "<<b<<" a photon energy of "<<photon_energies[b]/ev_to_K/kb<<" eV "<<endl;
        }

        cout<<" WAVELENGTH GRID FINISHED. num_bands_in / num_he_bands / num_bands_out = "<<num_bands_in<<" / "<<num_he_bands<<" / "<<num_bands_out<<endl;

        /*
        cout<<"   ATTENTION: Bands go from "<<l_i[0]<<" till "<<l_i[num_bands]<<", while lmin/lmax was given as "<<lambda_min<<"/"<<lambda_max<<endl;
        cout<<"   Please confirm this is correct."<<endl;
        char tmpchar;
        cin>>tmpchar;
        */
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
        else if(type_of_grid==2) { //power-law grid with higher log-resolution near the planet
            
            double dlogx  = pow(10.,1./cells_per_decade);
            double dlogx2 = pow(10.,1./grid2_cells_per_decade);
            double templogx = dlogx;
            //int grid2_transition_i = 0;
            
            x_i[0] = domain_min / std::pow(dlogx, (num_ghosts-1)) ;
            cout<<" grid_transition_i = "<<grid2_transition_i<<endl;
            
            for(int i=1; i<= num_cells; i++) {
                
                if(i <= grid2_transition_i) {
                    
                    templogx = dlogx;
                } else if (i >= (grid2_transition_i + 5)) {
                    templogx = dlogx2;
                } else {
                    templogx = dlogx + double(i-grid2_transition_i)/5. * (dlogx2 - dlogx);
                }
                    
                x_i[i]   =  x_i[i-1] * templogx;
            }
            
            //Assign the last boundary as domain maximum as long as nonuniform grid is in the test-phase
            //grid2_transition = domain_max;
            domain_max       = x_i[num_cells];
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

        R_core = x_i[1];
        
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
        
        int all_grid_ok = 1, broken_index;
        double broken_value;
        for(int i = 0; i <= num_cells+1; i++) {
            if( dx[i] < 0 ) {
                all_grid_ok = 0;
                broken_index = i;
                broken_value = dx[i];
            }
        }
            
        
        if(debug > 0) cout<<"Init: Finished grid setup"<<endl;
        if(!all_grid_ok) {
            std::stringstream err ;
            err <<"WARNING! Problem with grid values. Example broken value, dx["<<broken_index<<"]="<<broken_value<<endl;
            err<<" FIRST AND LAST DX = "<<dx[0]<<" / "<<dx[num_cells+1]<<endl;
            throw std::runtime_error(err.str()) ;
        }
        cout<<endl;
        
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
        K_zz            = np_somevalue(num_cells+2, K_zz_init);
        
        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        ////
        //
        // Create and initialize the species
        //
        ////
        ///////////////////////////////////////////////////////////////////////// 
        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        num_species = read_parameter_from_file<int>(filename,"PARI_NUM_SPECIES", debug, 1).value; //Read the number of species
        if(num_species<1) {
            std::stringstream err;
            err<<" Number of species has to be >= 1!"<<endl;
            throw std::invalid_argument(err.str()) ;
        }

        if (NUM_SPECIES != Eigen::Dynamic && num_species != NUM_SPECIES) {
            std::stringstream err ;
            err << "Error: You compiled aiolos with a different number of species requested"
                << " to that requested in the parameter file.\n You must recompile with either: "
                << " 1) a dynamic number of species (default) or 2) the correct number of species";
            throw std::invalid_argument(err.str()) ;
        }

        init_T_temp   = read_parameter_from_file<double>(filename,"INIT_T_TEMP", debug, 1.).value; //Overwrite initial temperature with this value. Use discouraged. 
        if(debug > 0) cout<<"Init: About to setup species with num_species = "<<num_species<<endl; 
        
        species.reserve(num_species);
        
        for(int s = 0; s < num_species; s++) {
            species.push_back(c_Species(this, filename, speciesfile, s, debug)); // Here we call the c_Species constructor
        }
        cout<<endl;
        
        if(use_self_gravity) {
            for(int cnt=0; cnt<50; cnt++) {
                //cout<<"Before update mass and pot"<<endl;
                update_mass_and_pot();
                
                //cout<<"Before second run in update dens"<<endl;
                for(int s = 0; s < num_species; s++) {
                    //species[s].update_kzz_and_gravpot(s); 
                    
                    for(int i=0; i<num_cells+2; i++) {
                        species[s].phi_s[i] = phi[i];
                    }
                    
                    species[s].initialize_hydrostatic_atmosphere(filename);
                }
            }
        }
        
        //Finally, look for the species indices for some important species. Stored in the utility integers e_idx, hnull_idx, hplus_idx etc.
        init_highenergy_cooling_indices();
        
        ///////////////////////////////////////////////////////////////////////// 
        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        ////
        //
        // After successful init, compute all primitive quantities and determine first timestep
        //
        ////
        ///////////////////////////////////////////////////////////////////////// 
        /////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure(species[s].u);
        
        dt = t_max ; // Old value needed by get_cfl_timestep()
        dt = get_cfl_timestep();
        
        if(debug > 0) cout<<"Init: Finished Init. Got initial dt = "<<dt<<" This is only temporary dt_init, not used for the first timestep."<<endl;
        //
        // Matrix init via Eigen
        //

        if(debug > 0) cout<<"Init: Setting up friction workspace."<<endl;
        
        alphas_sample   = Eigen::VectorXd::Zero(num_cells+2);
        friction_sample = Eigen::VectorXd::Zero(num_cells+2);
        
        reaction_matrix = Matrix_t::Zero(num_species, num_species);
        reaction_b      = Vector_t(num_species);
        
        reaction_matrix_ptr = new Eigen::MatrixXd[omp_get_max_threads()];
        reaction_b_ptr      = new Eigen::VectorXd[omp_get_max_threads()];
        LUchem_ptr          = new Eigen::PartialPivLU<Matrix_t>[omp_get_max_threads()];
        
        for(int i =0; i<omp_get_max_threads(); i++ ) {
            reaction_matrix_ptr[i] = Matrix_t::Zero(num_species, num_species);
            reaction_b_ptr[i]      = Vector_t(num_species);
            //LUchem_ptr[i]          = Eigen::PartialPivLU<Matrix_t>;
        }
        //n_news          = Vector_t(num_species);
        
        for(int s = 0; s < num_species; s++) {
            species[s].update_kzz_and_gravpot(s);
        }
        
        if(num_species >= 1) {
            friction_matrix_T     = Matrix_t::Zero(num_species, num_species);
            friction_coefficients = Matrix_t::Zero(num_species, num_species);
            friction_coeff_mask   = Matrix_t::Ones(num_species, num_species);
            inv_totmasses         = Matrix_t::Zero(num_species, num_species);
            resonant_pair_matrix  = Matrix_t::Ones(num_species, num_species);
            friction_vec_input    = Vector_t(num_species);
            friction_vec_output   = Vector_t(num_species);
            friction_dEkin        = Vector_t(num_species);
            dens_vector           = Vector_t(num_species);
            numdens_vector        = Vector_t(num_species);
            mass_vector           = Vector_t(num_species);
            temperature_vector    = Vector_t(num_species);
            temperature_vector_augment    = Vector_t(num_species);

            identity_matrix       = Matrix_t::Identity(num_species, num_species);
            unity_vector          = Vector_t::Ones(num_species);
            LU                    = decltype(LU)(num_species) ;
        }
        else {
            friction_coeff_mask = Matrix_t::Ones(num_species, num_species);
        }
        
        //
        // Dust species dont feel each other: Punch holes in the coefficient mask for the dust species, so they don't feel each other
        // Also set the total inv mass matrix values
        //
        for(int si=0; si < num_species; si++)
            for(int sj=0; sj < num_species; sj++) {
                inv_totmasses(si,sj) = 1./(species[si].mass_amu*amu + species[sj].mass_amu*amu);
                
                if(species[si].is_dust_like == 1 && species[sj].is_dust_like == 1) {
                    friction_coeff_mask(si,sj) = 0.;
                    friction_coeff_mask(sj,si) = 0.;
                }
                
            }
                
    
    radiation_matrix_T   = Matrix_t::Zero(num_species, num_species);
    radiation_matrix_M   = Matrix_t::Zero(num_species, num_species);
    radiation_vec_input  = Vector_t(num_species);
    radiation_vec_output = Vector_t(num_species);
    radiation_cv_vector  = Vector_t(num_species);
    radiation_T3_vector  = Vector_t(num_species);
    
    previous_monitor_J = std::vector<double>(num_bands_out) ;
    previous_monitor_T = std::vector<double>(num_species) ;
    
    solar_heating = Eigen::VectorXd::Zero(num_bands_in,  1);
    solar_heating_final = Eigen::VectorXd::Zero(num_bands_in,  1);
    S_band        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    dS_band       = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    dS_band_special=Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    dS_band_zero  = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    planck_matrix = Eigen::MatrixXd::Zero(num_plancks, 2);
    
    //
    // Compute Planck matrix based on loggrid and 100 K blackbody
    //
    double sum_planck = 0;
    double dsum_planck;
    double lmin, lmax;
    double dlogx = pow(lmaxglobal/lminglobal, 1./((double)num_plancks));
    lT_spacing     = dlogx;
    //double lT_crit = 14387.770; //=hc/k_b, in units of mum K
    double lumi100 = sigma_rad * pow(100.,4.) / pi;
    
    for(int p = 0; p < num_plancks; p++) {

        lmin      = lminglobal * pow(dlogx,(double)p);
        lmax      = lmin * dlogx;
        
        dsum_planck = compute_planck_function_integral2(lmin, lmax, 100.) / lumi100;
        sum_planck += dsum_planck;
        
        planck_matrix(p,0) = lmax*100.;  // Wavelength
        planck_matrix(p,1) = sum_planck; // Planck integral until this matrix element
        
    }
 
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 1."<<endl;
    for(int s=0; s<num_species; s++) {
        species[s].dS = Eigen::VectorXd::Zero(num_cells+2,  1);
        species[s].dG = Eigen::VectorXd::Zero(num_cells+2,  1);
        species[s].dGdT = Eigen::VectorXd::Zero(num_cells+2,  1);
    }
    
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 2."<<endl;
    
    double templumi = 0;
    
    if(true) {
        
            for(int b=0; b<num_bands_in; b++) {
                if(debug > 0) cout<<"Init: Assigning stellar luminosities. b ="<<b<<endl;
                
                if(num_bands_in == 1) {
                    solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.);
                    solar_heating(b) += sigma_rad * pow(T_other,4.) * pow(R_other*rsolar,2.)/pow(d_other*au,2.);
                    solar_heating(b) += UV_star/(4.*pi*pow(planet_semimajor*au,2.));
                    solar_heating(b) += X_star/(4.*pi*pow(planet_semimajor*au,2.));
                    templumi += solar_heating(b);
                    
                    solar_heating_final(b) = solar_heating(b);
                    cout<<" Solar heating is "<<solar_heating(b)<<" T_star = "<<T_star<<" pow(R_star*rsolar,2.) "<<pow(R_star*rsolar,2.)<<" pow(planet_semimajor*au,2.) "<<planet_semimajor<<endl;
                }
                else{
                    cout<<"SOLAR HEATING in bin "<<b;
                    cout<<" from/to lmin/lmax"<<l_i_in[b];
                    cout<<"/"<<l_i_in[b+1]<<" with frac = "<<compute_planck_function_integral4(l_i_in[b], l_i_in[b+1], T_star);
                    
                    solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.) * compute_planck_function_integral4(l_i_in[b], l_i_in[b+1], T_star);
                    solar_heating(b) += sigma_rad * pow(T_other,4.) * pow(R_other*rsolar,2.)/pow(d_other*au,2.) * compute_planck_function_integral4(l_i_in[b], l_i_in[b+1], T_other);;
                    templumi         += solar_heating(b);
                    
                    //if (BAND_IS_HIGHENERGY[b] == 1) {
                    if(l_i_in[b+1] <= 0.09161) { //Detect the EUV band 
                            
                        if(l_i_in[b+1] <= 0.0161 ) { //Detect the X ray band
                            solar_heating(b) += X_star/(4.*pi*pow(planet_semimajor*au,2.));
                            cout<<endl<<"        band "<<b<<" detected as X band. Currently, this is treated as thermal band. Assigning X flux "<<X_star/(4.*pi*pow(planet_semimajor*au,2.))<<" based on X lumi "<<X_star;
                            BAND_IS_HIGHENERGY[b] = 0;
                        }
                        else {
                            solar_heating(b) += UV_star/(4.*pi*pow(planet_semimajor*au,2.));
                            cout<<endl<<"        band "<<b<<" detected as UV band. Assigning UV flux "<<UV_star/(4.*pi*pow(planet_semimajor*au,2.))<<" based on UV lumi "<<UV_star;
                        }
                    }
                    
                    solar_heating_final(b) = solar_heating(b);
                }
                
                cout<<" is F = "<<solar_heating(b)<<endl;
                
        }
        
    } 
    
    if(fluxfile.compare("---")!=0) {//else {
        
        for(int b = 0; b<num_bands_in; b++) {
            
            double tempenergy = 0.;
            int cnt = 0;
            
            double flux;
            double lam;
            
            ifstream file(fluxfile);
            string line;
            if(!file) {
                cout<<"Couldnt open flux spectrum file "<<filename<<"!!!!!!!!!!1111"<<endl;
            }
            
            while(std::getline( file, line )) {
    
                std::vector<string> stringlist = stringsplit(line," ");
                
                lam   = std::stod(stringlist[0]);
                flux  = std::stod(stringlist[1]);
                
                if(lam > l_i_in[b] && lam < l_i_in[b+1]) {
                    tempenergy += lam*flux;
                    cnt++;
                }
             }
             file.close();
             
             if(cnt > 0) {
                 tempenergy /= (double)cnt;
                 solar_heating(b) += tempenergy * fluxmultiplier;
                 templumi += solar_heating(b);
                 solar_heating_final(b) += solar_heating(b);
            } else {
                solar_heating(b) = 0.;
            }
            
            cout<<"SOLAR HEATING read from file in bin "<<b;
            cout<<" from/to lmin/lmax"<<l_i_in[b];
            cout<<" is F = "<<solar_heating(b)<<endl;
        }
        
    }
    
    cout<<"TOTAL SOLAR HEATING / Flux = "<<templumi<<" Luminosity = "<<(templumi*4.*pi*rsolar*rsolar*pi)<<" | T_irr  = "<<pow(templumi/sigma_rad,0.25) <<" T_eq = "<<pow(templumi/sigma_rad/4,0.25) <<" Teq/2**0.25 = "<<pow(templumi/sigma_rad/8,0.25)<<endl;
    
    double totallumi = 0;
    double totallumiincr = 0;
    for(int b=0; b < num_bands_in; b++) {
        
        totallumiincr = compute_planck_function_integral4(l_i_in[b], l_i_in[b+1], T_star);
        totallumi     += totallumiincr;
        if(debug > 0 ) cout<<" INIT BANDS, b = "<<b<<" l_i_in[b] = "<<l_i_in[b]<<" l_i[b+1] = "<<l_i_in[b+1]<<" l_i12[b] = "<<l_i12_in[b]<<" fraction = "<<totallumiincr<<" fraction for T=1430 K = "<<compute_planck_function_integral3(l_i_in[b], l_i_in[b+1], 1430.)<<endl;
    }
    cout<<" 2nd Luminosity check, L / sigma T^4/pi is = "<<totallumi/compute_planck_function_integral4(l_i_in[0], l_i_in[num_bands_in-1], T_star)<<endl;
    

    // Initialize the planet's temperature if it has not already been done:
    if (T_planet == 0 && use_planetary_temperature) {
        T_planet = sigma_rad * (T_int*T_int) * (T_int*T_int) ;
        for (int b=0; b < num_bands_in; b++)
            T_planet += 0.25*solar_heating(b) ;      
        T_planet = pow(T_planet/sigma_rad, 0.25) ;
    }

    total_opacity        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    cell_optical_depth   = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    radial_optical_depth = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);

    total_opacity_twotemp        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    cell_optical_depth_twotemp   = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    cell_optical_depth_highenergy= Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    radial_optical_depth_twotemp = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    highenergy_switch            = Eigen::MatrixXd::Ones(num_species, num_bands_in);
    
    Etot_corrected = Eigen::MatrixXd::Zero(num_cells+2, 1);
    
    Jrad_FLD       = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    Jrad_init      = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    if(use_rad_fluxes == 1)
        tridiag        = BlockTriDiagSolver<Eigen::Dynamic>(num_cells+2, num_bands_out + num_species) ;
    else if(use_rad_fluxes == 2)
        tridiag        = BlockTriDiagSolver<Eigen::Dynamic>(num_cells+2, num_bands_out) ;
        
    double rade_l = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERL", debug, 1.).value; //Erad_left for radiative shock test
    double rade_r = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERR", debug, 1.).value; //Erad_right for radiative shock test
    init_J_factor = read_parameter_from_file<double>(filename,"INIT_J_FACTOR", debug, -1.).value; //Discontinued
    
    
    for(int j = 0; j <= num_cells+1; j++) {
        
        if(num_bands_out == 1)  {
            for(int s = 0; s< num_species; s++){
                
                Jrad_init(j,0) = rad_energy_multiplier * c_light /4. /pi;
                
                if(radiation_matter_equilibrium_test == 0)
                    Jrad_FLD(j,0)  = rad_energy_multiplier * sigma_rad*pow(species[s].prim[j].temperature,4) / pi;
                    
                else if(radiation_matter_equilibrium_test == 1)
                    Jrad_FLD(j,0) = Jrad_init(j,0);
                
                else if(radiation_matter_equilibrium_test == 3)
                {
                    if(j == (num_cells/2+1)) {
                        Jrad_FLD(j,0) = Jrad_init(j,0);
                        cout<<" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RADMATTER TEST 3, I am in cell "<<j<<" and it got J = "<<Jrad_FLD(j,0)<<endl;
                    }
                    
                    else
                        Jrad_FLD(j,0) = 0.;
                }
                
                else if(radiation_matter_equilibrium_test >= 4) {//Nonlinear diffusion test, and proper RHD shock tests
                    
                    double SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug, 0.5).value;//Location of discontinuity between left and right values in shock tube
                    
                    if(x_i12[j] < SHOCK_TUBE_MID) 
                        Jrad_FLD(j,0) = rade_l;
                    else
                        Jrad_FLD(j,0) = rade_r;
                    
                }
                
                if(debug > 3) {
                //if(debug > 1 && j==5) {
                    cout<<" Jrad("<<j<<","<<0<<") = "<<Jrad_FLD(j,0)<<" dJrad_species["<<s<<"] = "<<rad_energy_multiplier * compute_planck_function_integral3(l_i_out[0], l_i_out[1], species[s].prim[j].temperature)<<" T_rad = "<<pow(pi*Jrad_FLD(j,0)/sigma_rad,0.25)<<" at Tgas = "<<species[s].prim[j].temperature<<endl;
                }
                 
            }

        } else {
            
            for(int b = 0; b < num_bands_out; b++) {
                for(int s = 0; s< num_species; s++){
                    Jrad_FLD(j,b)  = rad_energy_multiplier * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], species[s].prim[j].temperature) * sigma_rad*pow(species[s].prim[j].temperature,4) / pi;
                    if(Jrad_FLD(j,b) < 1e-50)
                        Jrad_FLD(j,b) = 1e-50;
                    
                    Jrad_init(j,b) = Jrad_FLD(j,b);
                    if(debug > 1 && j==5) {
                        cout<<" Jrad("<<j<<","<<b<<") = "<<Jrad_FLD(j,b)<<" dJrad_species["<<s<<"] = "<<rad_energy_multiplier * compute_planck_function_integral3(l_i_out[b], l_i_out[b+1], species[s].prim[j].temperature)<<" at T ="<<species[s].prim[j].temperature<<endl; 
                     }
                }
            }
        }
        
        //
        // Set initial gradient for stability: FLD apparently doesn't like to be initialized in thermal equilibrium in the optically thin regions
        //
        if(j>10 && !closed_radiative_boundaries)
           Jrad_FLD(j,0) = Jrad_FLD(10,0) *  (x_i12[10]*x_i12[10]) / (x_i12[j]*x_i12[j]);  

    }
    cout<<"photochem level = "<<photochemistry_level<<endl;
    n_init = np_zeros(num_species);
    n_tmp  = np_zeros(num_species);
    
    for (int s = 0; s < num_species; s++) {
        n_init[s] = species[s].prim[num_cells-1].number_density;
        n_tmp[s]  = n_init[s];
    }
    
    if(photochemistry_level == 2) {
        cout<<" Init chemistry."<<endl<<endl;
        
        try {
            
                init_reactions(1);
        }
        catch(int i ) {
            
            if(i==0)
                cout<<" FATAL ERROR caught in init_chemistry: Reagular reaction educt index exceeds num_species! "<<endl;
            else if(i==1)
                cout<<" FATAL ERROR caught in init_chemistry: Reagular reaction product index exceeds num_species! "<<endl;
            else if(i==2)
                cout<<" FATAL ERROR caught in init_chemistry: Photochem reaction educt index exceeds num_species! "<<endl;
            else if(i==3)
                cout<<" FATAL ERROR caught in init_chemistry: Photochem reaction product index exceeds num_species! "<<endl;
            else if(i==4)
                cout<<" FATAL ERROR caught in init_chemistry: Photochem reaction band index exceeds num_bands! "<<endl;
            else
                cout<<"Misc integer error code caught in init_chemistry."<<endl;
            
        } 
        catch(...) {
            
            cout<<"Misc error caught in init_chemistry."<<endl;
        }
        
    }
    
    //c_Sim::c_Sim(string filename_solo, string speciesfile_solo, string workingdir, int debug, int debug_cell, int debug_steps) {
    write_into_execution_log(workingdir, filename_solo, speciesfile_solo);
    
    cout<<" Finished init. Returning to main now."<<endl<<endl;
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
 * Species class constructor.
 * 
 * Sets up the boundaries to be used with this boundary, equation of state, initialises primitive and conservative arrays.
 * Starts hydrostatic construction of temperature and density structures and initializes opacity data.
 * Most of the base data for this needs to be setup in the *spc file and linked via the opacity files.
 * 
 * @param[in] base_simulation Pointer to the c_Sim object containing this species class. This way one can quickly access important numbers by using base->important_number
 * @param[in] filename Parameter file string, including directory substrings
 * @param[in] species_filename Species file string, including directory substrings
 * @param[in] species_index Own species number to be saved within this object
 * @param[in] debug Debug level set either on command line or in param file.
 */
c_Species::c_Species(c_Sim *base_simulation, string filename, string species_filename, int species_index, int debug) {
    
        base               = base_simulation;
        this_species_index = species_index;
        num_cells          = base->num_cells;
        num_bands_in       = base->num_bands_in;
        num_bands_out      = base->num_bands_out;
        this->workingdir   = base->workingdir;
        this->debug        = debug;
        
        if(debug >= 0) cout<<"        Species["<<species_index<<"] = "<<speciesname<<" Beginning init."<<endl;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Physical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        boundary_left  = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_LEFT", debug, BoundaryType::fixed).value; //Which boundaries to use. See enum.h or this file, lines 1800
        boundary_right = read_parameter_from_file<BoundaryType>(filename,"PARI_BOUND_TYPE_RIGHT", debug, BoundaryType::fixed).value;

        const_T_space  = read_parameter_from_file<double>(filename,"PARI_CONST_TEMP", debug, 1.).value; //Temperature at boundary. Use depending on INIT_TEMPERATURE_MODEL
        TEMPERATURE_BUMP_STRENGTH    = read_parameter_from_file<double>(filename,"TEMPERATURE_BUMP_STRENGTH", debug, 0.).value; //Discontinued
        pressure_broadening_factor   = read_parameter_from_file<double>(filename,"PRESSURE_BROADENING", debug, 0.).value;  //For opacity model==P, weaken/strengthen pressure broadening
        pressure_broadening_exponent = read_parameter_from_file<double>(filename,"BROADENING_EXP", debug, 1.).value;       //For opacity model==P, opacity = (1+broad_factor * (p/1e5)^exponent)
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Finished reading boundaries."<<endl;
        if(debug > 0) cout<<"         Boundaries used in species["<<speciesname<<"]: "<<boundary_left<<" / "<<boundary_right<<endl;

        const_opacity   = read_parameter_from_file<double>(filename,"PARI_CONST_OPAC", debug, 1.).value; //Value of constant opacity, modified by individual multipliers for const opa model
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // //Readin species file, const_cv, const_gamma_adiabat for this species
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Reading species data..."<<endl;
        read_species_data(species_filename, species_index); 
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Done reading species data."<<endl;
        
        u               = init_AOS(num_cells+2); //Conserved hyperbolic variables: density, mass flux, energy density
        u0              = init_AOS(num_cells+2); //Conserved hyperbolic variables: density, mass flux, energy density
        dudt[0]         = init_AOS(num_cells+2);
        dudt[1]         = init_AOS(num_cells+2);
        source          = init_AOS(num_cells+2);  
        source_pressure = init_AOS(num_cells+2);
        flux            = init_AOS(num_cells+1);
        lconvect        =std::vector<double>(num_cells+1);

        //
        // Determine which equation of states we are using, and assign thermodynamic variables: heat capacity, adiabatic ratio, and equation of state
        //
        
        //Gas species
        if(is_dust_like == 0) {
            cv            = 0.5 * degrees_of_freedom * Rgas / mass_amu; 
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom;
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu*amu) ;
        }
        //Dust species
        else {
            cv            = 1.5 * Rgas / mass_amu; // 3 Degrees of freedom per atom (high T limit of debye theory)
            mass_amu      = mass_amu * degrees_of_freedom/3 ; // Convert to mass of grain in amu
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom; // Total degrees of freedom per grain
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu*amu) ;
        }
        
        if(debug >= 0) cout<<"        Species["<<species_index<<"] got a gamma_adiabatic = "<<gamma_adiabat<<" and cv = "<<cv<<" and a charge of "<<static_charge<<endl;

        u_analytic     = np_zeros(num_cells+2);

        prim   = std::vector<AOS_prim>(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        primlast = std::vector<AOS_prim>(num_cells+2); 
        de_e   = std::vector<double>(num_cells+2);
        prim_l = std::vector<AOS_prim>(num_cells+2);
        prim_r = std::vector<AOS_prim>(num_cells+2);
        temp_temperature = std::vector<double>(num_cells+2);
        phi_s            = std::vector<double>(num_cells+2);
        
        timesteps     = np_zeros(num_cells+2);		    
        timesteps_cs  = np_zeros(num_cells+2);	
        timesteps_rad = np_zeros(num_cells+2);
        finalstep     = np_zeros(num_cells+2);
        timesteps_de  = np_zeros(num_cells+2);
        
        K_zzf = std::vector<double>(num_cells+2);  
        this->update_kzz_and_gravpot(species_index);
        
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Initialize/readin scenario parameters
        //
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        //
        // Problem 1: two states, left and right, separated at a mid-position
        //
        if(base->problem_number == 1) {
            
            double u1l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1L", debug, 1.).value; //density left
            double u2l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2L", debug, 0.).value; //momentum left
            double u3l = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3L", debug, 1.).value; //total energy left
            
            double u1r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1R", debug, 1.).value; 
            double u2r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2R", debug, 0.).value;
            double u3r = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3R", debug, 1.).value;
            
            SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug, 0.5).value; //Position of discontinuity
            
            u1l *= initial_fraction;
            u1r *= initial_fraction;
            
            u3l *= initial_fraction;
            u3r *= initial_fraction;
            
            if (is_dust_like) { 
                u3l *= 3 / degrees_of_freedom ;
                u3r *= 3 / degrees_of_freedom ;
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
            
            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 1. with left rho/v/p/T ="<<pl.density<<"/"<<pl.speed<<"/"<<pl.pres <<"/"<<pl.temperature<<" init atmo "<<init_static_atmosphere  <<" density excess"<<density_excess<<endl;
        }
        //
        // Problem 2: A constant background state, with a planet embedded into it, free to evolve as physics dictates
        //
        else if(base->problem_number == 2) {
            
            double u1 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U1", debug, 1.).value; //Starting value for density construction.
            double u2 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U2", debug, 0.).value; //initial momentum. Should be 0, hence mostly irrelevant parameter.
            double u3 = read_parameter_from_file<double>(filename,"PARI_INIT_DATA_U3", debug, 1.).value; //Initial total energy. With 0 momentum and T given by construction, essentialy irrelevant parameter
            
            init_static_atmosphere = read_parameter_from_file<int>(filename,"PARI_INIT_STATIC", debug, 1).value; //Yesno, if yes, densities are constructed in hydrostatic manner.
            
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
                initialize_hydrostatic_atmosphere(filename);
                //TODO: PUTIN short Phi- (T,rho)-iteration to allow hydrostatic starts with significant self-garvit
            }
            else if (init_static_atmosphere == 2) {
                const_rho_scale = read_parameter_from_file<double>(filename,"CONST_RHO_SCALE", debug, 550.e5).value; //Const scale height if exponential atmosphere is constructed
                initialize_exponential_atmosphere();
            }
            
            //Initialize non-zero velocity
            if(base->init_sonic_radius > 0. && base->init_wind == 2){
                init_analytic_wind_solution();
            }
                

            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 2."<<endl;
        }
        else if (base->problem_number == 3) {
            initialize_sound_wave() ;
        }
        else if (base->problem_number == 0) {
            user_initial_conditions() ;
        } 
        else 
            initialize_default_test();
        
        
        USE_WAVE = read_parameter_from_file<int>(filename,"WAVE_USE", debug, 0).value; //Wave test
        if(USE_WAVE==1) {
            WAVE_AMPLITUDE = read_parameter_from_file<double>(filename,"WAVE_AMPLITUDE", debug).value; 
            WAVE_PERIOD    = read_parameter_from_file<double>(filename,"WAVE_PERIOD", debug).value; 
        }
        
        mass_reservoir = 0.;
        t_evap = 50.;
        latent_heat = 3.34e5*1e7/1e3; //3.34e5 J/kg
        p_sat = 1e6;
        
        // Apply boundary conditions
        apply_boundary_left(u) ;
        apply_boundary_right(u) ;
        
        if(debug > 0) cout<<"        Species["<<species_index<<"]: Init done."<<endl;
        
        opacity                = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out); //num_cells * num_bands_out
        opacity_planck         = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out); //num_cells * num_bands_out
        opacity_twotemp        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in); //num_cells * num_bands_in
        fraction_total_solar_opacity = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in); //num_cells * num_bands_in
        
        //NOTE: If the chosen opacity model is such that a *.opa file is read-in, then this initialization can be found in io.cpp, in void c_Species::read_species_data(string filename, int species_index)
    
}


/**
 * Initializer of a well-balanced hydrostatic density profile.
 * 
 * Starts with a guess of the temperature profile, either isothermal, adiabatic, or a combination thereof. With this, we compute the density structure, such that
 * it fulfills the hydrostatic pressure gradient, given the species' EOS. Integration happens either from outside in or inside out.
 * After the initial density profile is solved for, we apply a number of density cuts: For quasi-isothermal EOS we need a density jump to initialize a wind, 
 * for adiabatic+radiative EOS we need to cut the density to above ~1e-20, so that the rad solver doesn't crash. For simulations with ionisation/chemistry we 
 * initialize the setups as neutral and then have to force the electrons to follow the proton distribution, etc.
 * Another cut is applied to prevent heavy shocks in the case of multi-species simulations with ionisation: Beyond the sonic point, assume density drops off with sensible powerlaw
 * such that the tau=1 radius does not sit in the outer simulation domain (which would have a different density profile in a self-consistent solution).
 * 
 * @param[in] filename Parameter file filename, including the folder substrings.
 */
void c_Species::initialize_hydrostatic_atmosphere(string filename) {
    
    if(debug > 0) cout<<"            ATTENTION: Initializing hydrostatic construction for species "<<speciesname<<" and overwriting prior initial values."<<endl;
    
    double temp_rhofinal;
    double factor_inner, factor_outer;
    double T_inner;
    double T_outer;
    double metric_inner;
    //double residual;
    double metric_outer;
    
    double dphi_factor = read_parameter_from_file<double>(filename,"PARI_DPHI_FACTOR", debug, 1.).value; //Modifier on initial entropy
        
    if(base->photochemistry_level > 0 && std::abs(static_charge) > 0 )
        base->density_floor = read_parameter_from_file<double>(filename,"ION_FLOOR", debug, 1.e-20).value; 
    else
        base->density_floor = read_parameter_from_file<double>(filename,"DENSITY_FLOOR", debug, 1.e-20).value; //Caps minimum density at construction. Thermal radiation transport likes to crash at densities < 1e-20
    //
    // Start with the outermost cell and build up a hydrostatic atmosphere
    // Fulfilling KKM16, Eqn. 17
    //
    //long double residual  = 0.;
    
    //
    // First, initialize (adiabatic) temperature
    //
    for(int i=num_cells+1; i>=0; i--) {
            
        int mode;
        
        if(base->temperature_model == 'P') {
            mode = 1;
            if(i==num_cells+2) {
                prim[i].temperature = const_T_space;
            } else {
                mode = 11;
                
                if(base->x_i12[i] < 1.2e9) {
                    prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat) + const_T_space; 
                    //prim[i].temperature = const_T_space; 
                } else {
                    prim[i].temperature = - dphi_factor * base->get_phi_grav(base->x_i12[i], base->planet_mass) / (cv * gamma_adiabat) + const_T_space;
                    //prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat);
                    
                }
                    
            }
            
        } else if(base->temperature_model == '0') {
            mode = 2;
            
            if(this_species_index == 0) {
                mode = 22;
                
                if(base->x_i12[i] < 1.e99) {
                    mode = 221;
                    prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat) + const_T_space; //get_phi_grav(x_i12[i], enclosed_mass[0]);
                } else {
                    mode = 222;
                    prim[i].temperature = - dphi_factor * base->get_phi_grav(base->x_i12[i], base->planet_mass) / (cv * gamma_adiabat) + const_T_space;
                    //prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat) + const_T_space; get_phi_grav(x_i12[i], enclosed_mass[0]);
                }
            } else {
                mode = 23;
                prim[i].temperature = base->species[0].prim[i].temperature;
                
            }
            
        }
        else {
            mode = 3;
            //prim[i].temperature = const_T_space;
            //if(i>=num_cells-1)
                //prim[i].temperature = 1.0*const_T_space - 1e-4 * base->phi[i] / (cv * gamma_adiabat); //Need initial temperature gradient for FLD
                prim[i].temperature = const_T_space; // - 1e-5 * std::pow(base->x_i12[i], 0.5); //Need initial temperature gradient for FLD

        }
        //prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat) + const_T_space;
        
        if(debug>=2) {
            if(mode == 23) 
                cout<<"Just assigned T ="<<prim[i].temperature<<" to i = "<<i<<" with phi_grav = "<<base->phi[i]<<" in mode = "<<mode<<" T[0] = "<< base->species[0].prim[i].temperature<<endl;
            else
                cout<<"Just assigned T ="<<prim[i].temperature<<" to i = "<<i<<" with phi_grav = "<<base->phi[i]<<" in mode = "<<mode<<endl;
            
            //char a;
            //cin>>a;
            
        }
    }
    
    //cout<<" POS1 dens[2] = "<<u[2].u1<<" temp[2] = "<<prim[2].temperature<<" num_cells+2 = "<<num_cells+2<<endl;
    
    //At this point, the right ghost cell is already initialized, so we can just build up a hydrostatic state from there
    int iter_start;
    if(base->type_of_grid==-2)
        iter_start = base->grid2_transition_i;
    else
        iter_start = num_cells;
    
    if(base->reverse_hydrostat_constrution == 0) {
        
            for(int i=iter_start; i>=0; i--)  {
            
            //
            // Construct next density as to fulfil the hydrostatic condition
            //
            T_outer = prim[i+1].temperature ;
            T_inner = prim[i].temperature ;

            //factor_outer = (gamma_adiabat-1.) * cv * T_outer; 
            //factor_inner = (gamma_adiabat-1.) * cv * T_inner; 
            
            eos->get_p_over_rho_analytic(&T_outer, &factor_outer); //sets factor = k*T/mu in the case of an ideal gas EOS
            eos->get_p_over_rho_analytic(&T_inner, &factor_inner);
            
            //double dphi =  (base->phi[i+1] * K_zzf[i+1] - base->phi[i] * K_zzf[i]);
            double dphi = phi_s[i+1] - phi_s[i];
            
            metric_outer = dphi * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
            metric_inner = dphi * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
            
            temp_rhofinal = u[i+1].u1 * (factor_outer + metric_outer) / (factor_inner - metric_inner);
            
            u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner);
            
            if(temp_rhofinal < 0) {
                
                char a;
                cout.precision(16);
                cout<<"NEGATIVE DENSITY IN INIT HYDROSTATIC i="<<i<<"/"<<num_cells<<endl;
                if((factor_inner - metric_inner) < 0) {
                    
                    cout<<"     Negative denominator detected. Product (gamma-1)*cv*T_inner is too small for chosen discretization."<<endl;
                    
                }
                cout<<endl;
                cout<<"     dphi debug: phi["<<i+1<<"] = "<<base->phi[i+1] * K_zzf[i+1]<<" phi["<<i<<"] = "<<base->phi[i]*K_zzf[i]<<" num_cells = "<<num_cells<<endl;
                cout<<"     metric_inner debug: dPhi = "<<(base->phi[i+1]*K_zzf[i+1] - base->phi[i]*K_zzf[i])<<" dx[i+1]+dx[i] = "<<(base->dx[i+1] + base->dx[i])<<" omega*dx  = "<<base->omegaminus[i]  * base->dx[i]<<endl;
                cout<<"     factor_inner debug: gamma-1 = "<<(gamma_adiabat-1.)<<" cv = "<<cv<<" T_inner = "<<T_inner<<endl;
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
        
    } else {//REVERSE_HYDROSTAT_CONSTRUCTION
            int negdens = 0; 
            
            //Ghost cells
            for(int i =0; i<2; i++)
                u[i] = AOS(u[i].u1, 0., cv * u[i].u1 *prim[i].temperature);
            
            for(int i=1; i<=iter_start; i++)  {
            
            //
            // Construct next density as to fulfil the hydrostatic condition
            //
            T_outer = prim[i+1].temperature ;
            T_inner = prim[i].temperature ;
            
            eos->get_p_over_rho_analytic(&T_outer, &factor_outer);
            eos->get_p_over_rho_analytic(&T_inner, &factor_inner);
            
            double dphi = phi_s[i+1] - phi_s[i];
            //double dphi = base->phi[i+1]*K_zzf[i+1] - base->phi[i]*K_zzf[i];
            //double dphi = base->phi[i+1] - base->phi[i];
            
            metric_outer = dphi * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
            metric_inner = dphi * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
            
            temp_rhofinal = u[i].u1 *  (factor_inner - metric_inner)/(factor_outer + metric_outer);
            
            if(dphi < 0.)
                temp_rhofinal = u[i].u1;
            
            //cout<<" s/i = "<<speciesname<<"/"<<i<<"   metric_inner debug: dPhi = "<<dphi<<" rhonew/rhoold = "<<temp_rhofinal<<"/"<<u[i].u1<<endl;
            
            //if(speciesname.compare("S2")==0) { //Force light electrons to be the same number as protons
            if(mass_amu < 0.90) {    
            //if(speciesname.compare("e-")==0 ) { //Force light electrons to be the same number as protons
                
                //int p_index = base->get_species_index("H+");
                int p_index = base->get_species_index("S1 H+ p+ p", 0);
                if(p_index < 0) {
                    cout<<" IN HYDROSTAT CONSTRUCTION for e-: Protons (H+) species not found! "<<endl;
                }
                temp_rhofinal = base->species[p_index].u[i+1].u1 * initial_fraction / base->species[p_index].initial_fraction;
                
            }
            
            //////////////////////////////////////////////////////////////////////////////////////
            //double floor = base->density_floor / mass_amu * std::pow(base->x_i12[i]/base->x_i12[1], -4.);
            
            double floor          = base->density_floor * std::pow(base->x_i12[i]/base->x_i12[1], -4.) * initial_fraction;
            //double lastval_scaled = u[i].u1 * std::pow(base->x_i12[i]/base->x_i12[i+1], 8.);
            double lastval_scaled = u[i].u1 * u[i].u1/u[i-1].u1;
            
            if(temp_rhofinal < 0.) {
                
                if(lastval_scaled > floor && base->use_init_discont_smoothing)
                    temp_rhofinal = lastval_scaled;
                else
                    temp_rhofinal = lastval_scaled; //this was formerly triggered when  //if( (temp_rhofinal < floor) || ( base->x_i12[i] > 0.5*base->rhill  && base->use_tides == 2       ) ) {
                
                negdens = 1;
            }
            /*
            if(base->use_init_discont_smoothing && temp_rhofinal < 0. && temp_rhofinal > floor) {
                //if(temp_rhofinal < lastval_scaled)
                    temp_rhofinal = lastval_scaled;
            }
            
            //if( (temp_rhofinal < floor) || ( base->x_i12[i] > 0.5*base->rhill  && base->use_tides == 2       ) ) {
            if(false){
                
                if(temp_rhofinal < 0.)
                    negdens = 1;
                
                temp_rhofinal = floor;
            }*/
            //////////////////////////////////////////////////////////////////////////////////////
            
            //if(this_species_index == 2) {
            //    cout<<" electrons, i/rho = "<<i<<"/"<<temp_rhofinal<<endl;
            //}
            
            //cout<<" s/i = "<<speciesname<<"/"<<i<<"   metric_inner debug: dPhi = "<<(base->phi[i+1]*K_zzf[i+1] - base->phi[i]*K_zzf[i])<<" rhonew/rhoold = "<<temp_rhofinal<<"/"<<u[i].u1<<" kzz = "<<K_zzf[i+1]<<"/"<<K_zzf[i]<<" To/Ti = "<<T_outer<<"/"<<T_inner<<endl;
            
            u[i+1] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_outer);
            
            if(debug > 2) {
            //if(i < 3) {
                
                char a;
                cout.precision(16);
                cout<<"NEGATIVE DENSITY IN INIT HYDROSTATIC i="<<i<<"/"<<num_cells<<endl;
                if((factor_inner - metric_inner) < 0) {
                    
                    cout<<"     Negative denominator detected. Product (gamma-1)*cv*T_inner is too small for chosen discretization."<<endl;
                    
                }
                cout<<endl;
                cout<<"     dphi debug: phi["<<i+1<<"] = "<<base->phi[i+1]*K_zzf[i+1]<<" phi["<<i<<"] = "<<base->phi[i]*K_zzf[i]<<" num_cells = "<<num_cells<<endl;
                cout<<"     metric_inner debug: dPhi = "<<(base->phi[i+1]*K_zzf[i+1] - base->phi[i]*K_zzf[i])<<" dx[i+1]+dx[i] = "<<(base->dx[i+1] + base->dx[i])<<" omega*dx  = "<<base->omegaminus[i]  * base->dx[i]<<endl;
                cout<<"     factor_inner debug: gamma-1 = "<<(gamma_adiabat-1.)<<" cv = "<<cv<<" T_inner = "<<T_inner<<endl;
                cout<<"     metric_outer = "<< metric_outer << " metric_inner = "<<metric_inner <<endl;
                cout<<"     factor_outer = "<< factor_outer << " factor_inner = "<<factor_inner <<endl;
                cout<<"     factor_outer+metric_outer = "<< (factor_outer + metric_outer) << " factor_inner-metric_inner = "<<( factor_inner - metric_inner) <<endl;
                cout<<"     RATIO = "<< ((factor_outer + metric_outer)/ (factor_inner - metric_inner)) <<endl;
                //cout<<"In hydostatic init: factor_dens = "<< (2.* factor_outer / delta_phi + 1.) / (2. * factor_inner / delta_phi - 1.) <<endl;
                cout<<"     Ratio of densities inner/outer = "<< temp_rhofinal/u[i+1].u1 <<endl;
                cout<<"     Ratio of temperatures inner/outer = "<<T_inner/T_outer<<" t_inner ="<<T_inner<<" t_outer ="<<T_outer<<endl;
                cout<<"     Ratio of pressures inner/outer = "<<cv * temp_rhofinal * T_inner /u[i+1].u3<<endl;
                cout<<"     Resulting density == "<<temp_rhofinal<<endl;
                cout<<"     floor ="<<floor<<endl;
                cout<<"     density before "<<u[i+1].u1<<endl;
                cin>>a;
            }
        }
        
        if(negdens) 
            cout<<" Negative densities in init_hydrostatic. Check conditions for hydrostatic construction and/or debug."<<endl;
        
        //u[1] = AOS(u[2].u1 , 0., u[2].u3);
        //u[0] = AOS(u[3].u1 , 0., u[3].u3);
        //prim[1] = prim[2];
        //prim[0] = prim[3];
    }
    
    //cout<<" POS2 dens[2] = "<<u[2].u1<<" temp[2] = "<<prim[2].temperature<<" rhoe1, rhoe2 = "<<u[1].u3<<"/"<<u[2].u3<<" cv ="<<cv<<endl;
    //cout<<"Assigned densities in num_cell-1 = "<<u[num_cells-1].u1<<endl;
    //cout<<"Assigned densities in num_cell = "<<  u[num_cells].u1<<endl;
    //cout<<"Assigned densities in num_cell+1 = "<<u[num_cells+1].u1<<endl;
    
    if(base->type_of_grid == -2) {
        
        for(int i = iter_start+1; i<num_cells+1; i++) {
            double rr = pow(base->x_i[i]/base->x_i[iter_start],5.);
            u[i] = AOS(u[iter_start].u1 / rr, 0., cv * u[iter_start].u1 / rr * prim[i].temperature) ;
        }
    }
    
    for(int i = 0; i<num_cells+1; i++) {
            //double floor = base->density_floor / mass_amu * std::pow(base->x_i12[i]/base->x_i12[1], -4.);
            double floor = base->density_floor * std::pow(base->x_i12[i]/base->x_i12[1], -4.) * initial_fraction;
            if(u[i].u1 < floor) {
                //floor = 1e6 * mass_amu*amu / (kb * prim[i].temperature); 
                //cout<<"FIXING SMALL DENSITIES in i ="<<i<<" for species = "<<this_species_index<<endl;
                u[i] = AOS(floor, 0., cv * floor * prim[i].temperature) ;
            }
    }
    
    if(base->init_T_temp > 2.7) {
        
            if(this_species_index == 0)
                cout<<" In INIT: Overwriting T with user-defined temperature:"<<base->init_T_temp;
                
            for(int i=num_cells; i>=0; i--)  {
                
                if(base->x_i12[i]>1.0e10) {
                    
                    //double temprho = u[i].u1;
                
                    //u[i] = AOS(temprho, 0., cv * temprho * base->init_T_temp) ;
                }
                
                //prim[i].temperature = base->init_T_temp;
            }
        }
        //else
        //    cout<<" In INIT: Did NOT overwrite T, because init_T_temp < 2.7= = "<<base->init_T_temp;
    
    
    if(u[2].u1 > 1e40) {
        cout<<"    ---WARNING---"<<endl;
        cout<<"    IN CONSTRUCT HYDROSTATIC:"<<endl;
        cout<<"    species["<<speciesname<<"] has reached a very high initial density  of > 1e40 at the inner boundary. Density = "<<u[2].u1<<endl;
        cout<<"    The well-balanced scheme seems to be unable to compensate minor velocity fluctuations under those circumstances."<<endl;
        cout<<"    Be advised, your solution is probably about to unphysically explode."<<endl;
    }
    cout<<endl;
    
    compute_pressure(u);

    
    for(int i=num_cells; i>=0; i--)  {
        primlast[i].internal_energy = prim[i].internal_energy;
    }
    //
    // Apply multiplicators for over/underdensities just below and above the sonic point
    //
    
    if(base->init_wind==1) {
        double csi              = prim[num_cells-1].sound_speed;
        double smoothed_factor;
        
        
        if(base->init_sonic_radius>0)
            bondi_radius = base->init_sonic_radius;
        else
            bondi_radius = base->planet_mass*G / (pow(csi,2.)); //Most restrictive sonic point, as we are not quite isothermal
        
        for( int i = 1; i < num_cells + 1; i++) {
            smoothed_factor = density_excess * (1./(1.+std::exp(-(base->x_i12[i] - bondi_radius)/(0.01 * bondi_radius )) ) ); 
            
            u[i].u1 = u[i].u1 * (1. + smoothed_factor);
            
            u[i].u3 = u[i].u3 * (1. + smoothed_factor);
            if(i==1 && debug > 0)
                cout<<"In initwind: First smoothed factor = "<<smoothed_factor<<endl;
        }
        
        cout<<"......... Ended hydrostatic construction for species "<<speciesname<<", last density was "<<u[num_cells+1].u1<<" while  density_excess="<<density_excess<<" r_b="<<bondi_radius<<endl<<endl;
    }
    else
        if(debug > 0)
            cout<<"@@@ Ended hydrostatic construction for species "<<speciesname<<endl<<endl;
    
    //TODO:Give some measure if hydrostatic construction was also successful, and give warning if not (e.g. when resolution too low and floors were reached).

    //cout<<"End of hydrostat"<<endl;
}

/**
 * Historical initialiser for analytic isothermal density profiles, unused now.
 */
void c_Species::initialize_exponential_atmosphere() {
    
    for(int i=num_cells+1; i>=0; i--) {
            prim[i].temperature = const_T_space;
    }
    
    for(int i=num_cells; i>=0; i--)  {
        double temp_rhofinal = u[i].u1*std::exp(-1./const_rho_scale*(base->x_i[i] - base->x_i[num_cells]));
        
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * prim[i].temperature) ;
    }
    cout<<endl<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<endl<<endl;
}


/**
 * Iterator skelleton for self-consistent density-temperature profiles. For users who wish to use fully self-consistent 
 * initial conditions as given by e.g. the Guillot2010 model.
 * 
 * @param[in] filename Parameter file filename, including the folder substrings.
 */
void c_Species::initialize_hydrostatic_atmosphere_iter(string filename) {
    cout<<" Initializing hydrostatic atmosphere based on filename = "<<filename<<endl;
}


/**
 * Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
 */
void c_Species::initialize_default_test() {
    
    if(debug > 0) cout<<"            Initialized default simulation."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = AOS(1,1,1);
}

/**
 * Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
 * 
 * @param[in] leftval Left hydrodynamic state (rho_L, mom_L, Energy_L)
 * @param[in] rightval Right hydrodynamic state (rho_R, mom_R, Energy_R)
 */
void c_Species::initialize_shock_tube_test(const AOS &leftval,const AOS &rightval) {
    
    if(debug > 0) cout<<"            Initialized shock tube test."<<endl;
    
    for(int i=0; i<=num_cells+1; i++) {
        
        if(base->x_i12[i] < SHOCK_TUBE_MID) 
            u[i] = leftval;
        else
            u[i] = rightval;
    }
}

/**
 * Initial conditions for one constant background state
 * 
 * @param[in] background Hydrodynamic state to be initialized in the entire domain (rho, mom, Energy)
 */
void c_Species::initialize_background(const AOS &background) {
    
    if(debug > 0) cout<<"            Initialized constant background state. state = "<<background.u1<<" "<<background.u2<<" "<<background.u3<<endl;
    
    for(int i=0; i<=num_cells+1; i++) 
        u[i] = background;
}

/**
 * Sound wave propagation test, to be compared against analytic model.
 */
void c_Species::initialize_sound_wave() {

    double amp = 1e-6 ;
    double rho0 = 1.0 * initial_fraction ;
    double cs = 1.0 * (base->domain_max - base->domain_min) ;
    double p0 = rho0*cs*cs/(gamma_adiabat) ;

    // Make dust / gas temperatures similar
    if (is_dust_like)
        p0 *= 3 / degrees_of_freedom ;

    for(int i=0; i<=num_cells+1; i++) {
        double kx = 2*M_PI * (i+0.5 - base->num_ghosts) / double(num_cells - 2*(base->num_ghosts-1)) ;
        double f = amp * std::sin(kx) ;

        AOS_prim p (rho0 *(1 + f), cs*f, p0*(1 + gamma_adiabat*f)) ;
        eos->compute_conserved(&p, &u[i], 1) ;
    }
}

/**
 * Build first order hydrodynamic boundary conditions
 * user = 0, open = 1, reflecting = 2, fixed = 3, periodic = 4, giantplanet = 5
 * 
 * @param[in] u Hydrodynamic state in last ghost cell before active domain begins
 */
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
                double dphi = (base->phi[i]*K_zzf[i] - base->phi[i-1]*K_zzf[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                dphi       *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                prim.pres   = prim.pres + prim.density * dphi ;
                prim.pres   = std::max( prim.pres, 0.0) ;
                prim.speed = u[i].u2/u[i].u1;
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
            for (int i=0; i < num_ghosts; i++) {
                
                if(base->problem_number == 1)
                    u[i] = SHOCK_TUBE_UL ;
                else {
                    AOS_prim prim;
                    prim.density = BACKGROUND_U.u1;
                    if(this->prim[2].speed > 0)
                        prim.speed   = this->prim[2].speed;
                    else
                        prim.speed   = 0.;
                    prim.temperature = const_T_space;
                    
                    eos->update_eint_from_T(&(prim), 1);
                    eos->update_p_from_eint(&(prim), 1);
                    eos->compute_conserved(&(prim), &(u[i]), 1) ; //Ghost cell u is fixed after init, only need to update p in prim
                }
            }
            
            break;
        case BoundaryType::periodic:
            for (int i=0; i < num_ghosts; i++) {
                int iact = Ncell + i;
                u[i] = u[iact];
                base->phi[i]   = base->phi[iact] ;
            }
            break;
        case BoundaryType::giantplanet: //This is currently identical to the left reflecting boundaries, assuming the left boundary is a hard core
            for (int i=0; i < num_ghosts; i++) {
                int igh =  num_ghosts-1 -i;
                int iact = num_ghosts   +i;
                u[igh]     = u[iact]; 
                u[igh].u2 *= -1;
                base->phi[igh]   = base->phi[iact] ;
            }
            break;
            break;
    }
}

/**
 * Build first order hydrodynamic boundary conditions
 * user = 0, open = 1, reflecting = 2, fixed = 3, periodic = 4, giantplanet = 5
 * 
 * @param[in] u Hydrodynamic state in last ghost cell after active domain ends
 */
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
                
                double dphi = (base->phi[i]*K_zzf[i] - base->phi[i-1]*K_zzf[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                dphi       *= (base->omegaplus[i]*base->dx[i] + base->omegaminus[i-1]*base->dx[i-1]) ;
                prim.pres = prim.pres - base->right_extrap_press_multiplier * prim.density * dphi ;    
                
                //Older variant, should be more precise but has weird pressure slope
                //double dphi2 = (phi_s[i] - phi_s[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                //dphi2 *= (prim.density * base->omegaplus[i]*base->dx[i] + this->prim[i].density * base->omegaminus[i-1]*base->dx[i-1]) ;
                //prim.pres = prim.pres - dphi2 ; 
                
                prim.pres = std::max( prim.pres, 1e-1*prim.pres) ; //TODO: Replace 1e-3 with an estimate for the max pressure jump in a adiabatic shock
                //prim.pres = std::max( prim.pres, 0.0) ;          //27.10.2021: Not in use anymore, due to this causing problems with negative temperatures. p=0 -> E = 0 -> T = 0  and negative after a bit of hydro
                
                eos->compute_conserved(&prim, &u[i], 1) ;
            }
            break ;
        case BoundaryType::reflecting:
            for (int i=0; i < num_ghosts; i++) {
                int iact = Ncell + num_ghosts-1 -i ;
                int igh = Ncell + num_ghosts +i;

                u[igh]     = u[iact]; 
                u[igh].u2 *= -1;
                base->phi[igh]   = base->phi[iact] ;
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
                base->phi[i]   = base->phi[iact] ; 
            }
            break;
        case BoundaryType::giantplanet:
            
            for (int i=Ncell+num_ghosts; i < Ncell+2*num_ghosts; i++) {
                AOS_prim prim ;
                eos->compute_primitive(&u[i-1],&prim, 1) ;
                
                //double dphi = (phi_s[i] - phi_s[i-1]) / (base->dx[i-1] + base->dx[i]) ;
                double r =base->x_i12[i];
                double mdot      = 4.*3.141592*prim.density*prim.speed*r*r;
                double freefallv = -std::sqrt(2.*G*base->planet_mass/r);
                    
                double mdotlimit = base->max_mdot; //1.9e+18 g/s = 1e-2 mearth/yr
                if(this->is_dust_like)
                    mdotlimit = 1e-4 * mdotlimit;
                    
                string cases;
                if(mdot < mdotlimit) {
                    prim.pres = prim.pres;// - dphi * (prim.density * base->omegaplus[i]*base->dx[i] + this->prim[i].density * base->omegaminus[i-1]*base->dx[i-1]); //Limiting case
                    prim.speed = freefallv;
                    prim.density = mdotlimit/freefallv/r/r;
                }
                else {
                    prim.pres = std::min(prim.pres * 2., 1.);// - prim.density * dphi ; //Forcing inflow
                }
                
                eos->compute_conserved(&prim, &u[i], 1) ;
                eos->compute_auxillary(&prim, 1);
            }
            break ;
    }
}

/**
 * Creates a wave at the inner boundary supposedly propagating upwards in the gravity well
 * 
 * @param[out] u Hydrodynamic state in boundary
 * @param[in] globalTime Time
 */
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

/**
 * Setting the base pointer for reaction objects and computing the total product masses
 */
void c_reaction::set_base_pointer(c_Sim *base_simulation) {
    base = base_simulation;
    
    for(int& pj : products) {
            this->products_total_mass += p_stoch[pj] * base->species[pj].mass_amu;
    }
}

/**
 * Setting the base pointer for reaction objects and computing the total product masses as well as the photoenergy split factors.
 */
void c_photochem_reaction::set_base_pointer(c_Sim *base_simulation) {
    base = base_simulation;
    
    energy_split_factor = 1.;
    double denom = 0.;
    for(int& pj : products) {
            this->products_total_mass += p_stoch[pj] * base->species[pj].mass_amu;
            
            //cout<<"In photochem setup, added species number pj = "<<pj<<" mass = "<<base->species[pj].mass_amu<<" and stoch = "<<p_stoch[pj]<<endl;
            
            energy_split_factor *= base->species[pj].mass_amu;
            double temp = 1.;
            for(int& pk : products) {
                    if(pk != pj)
                        temp *= base->species[pk].mass_amu;
            }
            //cout<<" energy split denom sum "<<temp<<endl;
            denom += temp;
    }
    //cout<<" nom = "<<energy_split_factor<<endl;
    energy_split_factor /= denom;
    
    //cout<<" energy split factor = "<<energy_split_factor<<" / 1/energy_split_factor = "<<1./energy_split_factor<<" denom = "<<denom<<endl;
    
    base->highenergy_switch(educts[0], this->band) = 0.;
}


/**
 * Constructor wrapper for the reaction class. If reaction rates are given via a,b,c, we ignore b and c and just pass a on as first guess of the reaction rate to the actual constructor.
 * Switches arguments in real constructor in a safe way in case this is the reverse reaction of another reaction.
 */
c_reaction::c_reaction(bool is_reverse, bool mtype_reaction, int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double a, double b, double c) {
    this->is_reverse_reac = is_reverse;
    this->mtype_reaction  = mtype_reaction;
    this->reac_a = a;
    this->reac_b = b;
    this->reac_c = c;
    
    if(!is_reverse_reac)
        c_reaction_real(num_species, e_indices, p_indices, e_stoch, p_stoch, reac_a);
    else
        c_reaction_real(num_species, p_indices, e_indices, p_stoch, e_stoch, reac_a);
    
} 

/**
 * Constructor wrapper for the reaction class. 
 * Switches arguments in real constructor in a safe way in case this is the reverse reaction of another reaction.
 */
c_reaction::c_reaction(bool is_reverse, bool mtype_reaction, int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double reaction_rate) 
{
    this->is_reverse_reac = is_reverse;
    this->mtype_reaction  = mtype_reaction;
    
    if(!is_reverse_reac)
        c_reaction_real(num_species, e_indices, p_indices, e_stoch, p_stoch, reaction_rate);
    else
        c_reaction_real(num_species, p_indices, e_indices, p_stoch, e_stoch, reaction_rate);
            
}

/**
 * Constructor for the reaction class, which takes care of thermochemistry. Writes all data about the reaction so that the reaction matrix can later be built via iteration over all reactions.
 * 
 * @param[in] num_species Needs to know the global number of species
 * @param[in] e_indices list of global species indices participating as reactants (educts)
 * @param[in] p_indices list of global species indices participating as products
 * @param[in] e_stoch list of educt stochiometric factors
 * @param[in] p_stoch list of product stochiometric factors 
 * @param[in] reaction_rate initial guess for the reaction rate (exact value doesn't matter as long as its not NaN or -inf), it's recomputed via update_reaction_rate(temperature) at every timestep.
 */
void c_reaction::c_reaction_real(int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double reaction_rate) {
    
    for(long unsigned int e=0; e<e_indices.size(); e++) {
        if(e_indices[e] >= num_species) {
            
            cout<<" FATAL ERROR caught in init_chemistry: Target reactant index exceeds num_species! num_species = "<<num_species<<endl;
            char wait;
            cin>>wait;
            throw 0;
        }
            
    }
    for(long unsigned int p=0; p<p_indices.size(); p++) {
        if(p_indices[p] >= num_species) {
            
            cout<<" FATAL ERROR caught in init_chemistry: Target product index exceeds num_species! num_species = "<<num_species<<endl;
            char wait;
            cin>>wait;
            throw 1;
        }
            
    }
        
    this->educts   = e_indices;
    this->products = p_indices;
            
    this->e_num   = e_indices.size();
    this->p_num   = p_indices.size();
    
    this->dndts   = std::vector<double>(num_species);
    
    this->e_stoch_i = std::vector<int>(num_species);//  e_stoch;
    this->p_stoch_i = std::vector<int>(num_species); // p_stoch;
        
    this->e_stoch = std::vector<double>(num_species); //e_stoch;
    this->p_stoch = std::vector<double>(num_species); //p_stoch;
    this->reac_a  = reaction_rate;
    this->r       = reaction_rate;
    this->products_total_mass = 0.; // Needs to be set after the base pointer is set.        

    for(int s=0; s<num_species; s++) {
        this->e_stoch[s] = 0;
        this->p_stoch[s] = 0;
    }
        
    delta_stoch = 0.;
    for(int e=0; e<e_num; e++) {
        this->e_stoch[e_indices[e]]   = e_stoch[e];
        this->e_stoch_i[e_indices[e]] = (int) (e_stoch[e]*1.01);
        delta_stoch += e_stoch[e];
    }
                
    for(int p=0; p<p_num; p++) {
        this->p_stoch[p_indices[p]]   = p_stoch[p];
        this->p_stoch_i[p_indices[p]] = (int) (p_stoch[p]*1.01);
        delta_stoch -= p_stoch[p];
    }
}

/**
 * Constructor for the photochem_reaction class which takes care of photochemistry. Writes all data about the reaction so that the reaction matrix can later be built via iteration over all reactions.
 * 
 * @param[in] num_species Needs to know the global number of species
 * @param[in] num_bands Global number of bands
 * @param[in] band      Band number at which this reaction starts to be active, i.e. it's active for all b<= num_bands (this is the implementation of the ionisation threshold!)
 * @param[in] e_indices list of global species indices participating as reactants (educts)
 * @param[in] p_indices list of global species indices participating as products
 * @param[in] e_stoch list of educt stochiometric factors
 * @param[in] p_stoch list of product stochiometric factors 
 * @param[in] branching Branching ratio. A number <= 1, multiplier on the photoionisation rate
 * @param[in] threshold This is the ionisation threshold energy, user needs to make sure that this is consistent with the wavelength limit l_i_in[band + 1].
 */
c_photochem_reaction::c_photochem_reaction(int num_species, int num_bands, int band, std::vector<int> e_indices, std::vector<int> p_indices, std::vector<double> e_stoch, std::vector<double> p_stoch, double branching, double threshold)
{
        //cout<<"in init photochem reaction "<<endl;
        for(long unsigned int e=0; e<e_indices.size(); e++) {
            //cout<<" e = "<<e_indices[e]<<endl;
            if(e_indices[e] >= num_species) {
                cout<<" FATAL ERROR caught in init_chemistry: Target photochem reactant index exceeds num_species! "<<endl;
                char wait;
                cin>>wait;
                throw 2;
            }
        }
        for(long unsigned int p=0; p<p_indices.size(); p++) {
            //cout<<" p ="<<p_indices[p]<<endl;
            if(p_indices[p] >= num_species) {
                cout<<" FATAL ERROR caught in init_chemistry: Target photochem product index exceeds num_species! "<<endl;
                char wait;
                cin>>wait;
                throw 3;
            }
                
        }
        //cout<<"e_indices = "<<e_indices<<" p_indices = "<<p_indices<<endl;
        //cout<<"pos1"<<endl;
        //this->reaction_number = reacnum;
        this->educts   = e_indices;
        this->products = p_indices;
        //cout<<"pos2"<<endl;
        this->e_num   = e_indices.size();
        this->p_num   = p_indices.size();
        //cout<<"pos3, enum/pnum = "<<this->e_num<<" "<<this->p_num<<endl;
        this->e_stoch = std::vector<double>(num_species);//  e_stoch;
        this->p_stoch = std::vector<double>(num_species); // p_stoch;
        //cout<<"pos4"<<endl;
        
        for(int s=0; s<num_species; s++) {
            this->e_stoch[s] = 0.;
            this->p_stoch[s] = 0.;
        }
        
        for(int e=0; e<e_num; e++) {
            this->e_stoch[e_indices[e]] = e_stoch[e];
            //cout<<" set e_stoch = "<<e_stoch[e]<<endl;
            //cout<<" set e_indices = "<<e_indices[e]<<endl;
        }
        
        count_p = 0;
        for(int p=0; p<p_num; p++) {
            this->p_stoch[p_indices[p]] = p_stoch[p];
            count_p ++;
            //cout<<" set p_stoch = "<<p_stoch[p]<<endl;
            //cout<<" set p_indices = "<<p_indices[p]<<endl;
        }
            
        this->products_total_mass = 0.; // Needs to be set after the base pointer is set.
        
        //cout<<"pos6"<<endl;
        if(band >= num_bands) {
            cout<<" FATAL ERROR IN SETTING UP PHOTOCHEM: band = "<<band<<" >= num_bands = "<<num_bands<<endl;
            throw 4;
        }
            
        
        this->band = band; //TODO: Check that band <= num_bands
        this->branching_ratio = branching;
        this->threshold_energy = threshold * ev_to_K * kb ;
        //cout<<"pos7"<<endl;
}

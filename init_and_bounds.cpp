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

c_Sim::c_Sim(string filename_solo, string speciesfile_solo, string workingdir, int debug) {

        if(debug > 0) cout<<"Init position 0."<<endl;
        
        steps = -1;
        this->debug      = debug ;
        simname          = filename_solo;
        this->workingdir = workingdir;
        
        string filename    = workingdir + filename_solo;
        string speciesfile = workingdir + speciesfile_solo;
        
        cout<<" Opening parameter file "<<filename<<endl;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //  Numerical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double dx0;
        type_of_grid     = read_parameter_from_file<int>(filename,"PARI_GRID_TYPE", debug, 0).value;
        domain_min       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MIN", debug).value;
        domain_max       = read_parameter_from_file<double>(filename,"PARI_DOMAIN_MAX", debug).value;
        geometry         = read_parameter_from_file<Geometry>(filename, "PARI_GEOMETRY", debug, Geometry::cartesian).value;
        order            = read_parameter_from_file<IntegrationType>(filename, "PARI_ORDER", debug, IntegrationType::second_order).value;
        
        lambda_min       = read_parameter_from_file<double>(filename,"PARI_LAM_MIN", debug, 1e-1).value;
        lambda_max       = read_parameter_from_file<double>(filename,"PARI_LAM_MAX", debug, 10.).value;
        lambda_per_decade= read_parameter_from_file<double>(filename,"PARI_LAM_PER_DECADE", debug, 10.).value;
        T_star           = read_parameter_from_file<double>(filename,"PARI_TSTAR", debug, 5777.).value;
        R_star           = read_parameter_from_file<double>(filename,"PARI_RSTAR", debug, 1.).value;
        UV_star          = read_parameter_from_file<double>(filename,"PARI_UVSTAR", debug, 1.).value;
        num_bands        = read_parameter_from_file<int>(filename,"PARI_NUM_BANDS", debug, 1).value;
        T_core           = read_parameter_from_file<double>(filename,"PARI_TPLANET", debug, 200.).value;
        use_planetary_temperature = read_parameter_from_file<int>(filename,"USE_PLANET_TEMPERATURE", debug, 0).value;
        core_cv           = read_parameter_from_file<double>(filename,"PARI_CORE_CV", debug, 1.e9).value;
        
        radiation_matter_equilibrium_test = read_parameter_from_file<int>(filename,"RAD_MATTER_EQUI_TEST", debug, 0).value;
        radiation_diffusion_test_linear   = read_parameter_from_file<int>(filename,"RAD_DIFF_TEST_LIN", debug, 0).value;
        radiation_diffusion_test_nonlinear= read_parameter_from_file<int>(filename,"RAD_DIFF_TEST_NLIN", debug, 0).value;
        
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
            
            if(type_of_grid==2) {
                    grid2_transition       = read_parameter_from_file<double>(filename,"GRID2_TRANSITION", debug, 99.e99).value;
                    grid2_cells_per_decade = read_parameter_from_file<double>(filename,"GRID2_CELLS_PER_DECADE", debug, 10).value;
                    
                    num_cells          = (int) ( (log10f(grid2_transition) - log10f(domain_min)) * cells_per_decade );
                    grid2_transition_i = num_cells + num_ghosts;
                    num_cells         += (int) ( (log10f(domain_max) - log10f(grid2_transition)) * grid2_cells_per_decade );
                    
            }
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
        CFL_break_time = read_parameter_from_file<double>(filename,"CFL_BREAK_TIME", debug, 
                                                          std::numeric_limits<double>::max()).value ;
        
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
        
        problem_number    = read_parameter_from_file<int>(filename,"PARI_PROBLEM_NUMBER", debug).value;
        use_self_gravity  = read_parameter_from_file<int>(filename,"PARI_SELF_GRAV_SWITCH", debug, 0).value;
        use_linear_gravity= read_parameter_from_file<int>(filename,"PARI_LINEAR_GRAV", debug, 0).value;
        use_rad_fluxes    = read_parameter_from_file<int>(filename,"PARI_USE_RADIATION", debug, 0).value;
        init_wind         = read_parameter_from_file<int>(filename,"PARI_INIT_WIND", debug, 0).value;
        alpha_collision   = read_parameter_from_file<double>(filename,"PARI_ALPHA_COLL", debug, 0).value;
        rad_energy_multiplier=read_parameter_from_file<double>(filename,"PARI_RAD_MULTIPL", debug, 1.).value;
        collision_model   = read_parameter_from_file<char>(filename,"PARI_COLL_MODEL", debug, 'C').value;
        opacity_model     = read_parameter_from_file<char>(filename,"PARI_OPACITY_MODEL", debug, 'C').value;
        temperature_model = read_parameter_from_file<char>(filename,"INIT_TEMPERATURE_MODEL", debug, 'P').value;
        friction_solver   = read_parameter_from_file<int>(filename,"FRICTION_SOLVER", debug, 0).value;
        do_hydrodynamics  = read_parameter_from_file<int>(filename,"DO_HYDRO", debug, 1).value;
        rad_solver_max_iter = read_parameter_from_file<int>(filename,"MAX_RAD_ITER", debug, 1).value;
        
        if(problem_number == 2)
            monitor_output_index = num_cells/2; //TODO: Replace with index of sonic radius for dominant species?
        else
            monitor_output_index = 1;
        
        if(problem_number==2) {
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug, 1).value; //in Earth masses
            planet_mass     *= mearth;
            planet_semimajor= read_parameter_from_file<double>(filename,"PARI_PLANET_DIST", debug, 1.).value; //in Earth masses
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
        
        if(num_bands < 1) {
            cout<<" In INIT RADIATION, invalid num_bands = "<<num_bands<<" changing to num_bands = 1."<<endl;
            
            num_bands = 1;
        }
        
        l_i       = np_zeros(num_bands+1);    // Wavelenght bin boundaries
        l_i12     = np_zeros(num_bands);      // Wavelength bin midpoints or averages
        
        l_i[0]         = lminglobal/100.;     //Transform wavelengths into microns
        l_i[num_bands] = lmaxglobal/100.;
        
        if(num_bands == 2) {
            l_i[1]   = 0.5*(lambda_max + lambda_min);
            l_i12[0] = lambda_min;
            l_i12[1] = lambda_max;
        }
        else if(num_bands > 2)
        {
            //double dlogl = pow(10., 1./lambda_per_decade);
            //double dlogl2 = pow(lambda_max/lambda_min, 1./(num_bands-2));
            l_i[1]             = lambda_min;
            l_i[num_bands-1]   = lambda_max;
            double dlogl2      = pow(lambda_max/lambda_min, 1./(num_bands-2));
            //l_i[0] = lambda_min;
            
            for(int b=2; b<num_bands-1; b++) {
                l_i[b]      = l_i[b-1] * dlogl2;
                cout<<" in NUM_BANDS>2, b = "<<b<<" l_i[b] = "<<l_i[b]<<endl;
            }
            
            for(int b=0; b<num_bands; b++) {
                l_i12[b]  = pow( 10., 0.5 * (std::log10(l_i[b]) + std::log10(l_i[b+1])));   
                
            }
            
        }
        
        
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
            
            x_i[0] = domain_min / std::pow(dlogx, (num_ghosts-1)) ;
            for(int i=1; i<= num_cells; i++) {
                
                if(x_i[i-1] < grid2_transition)
                    x_i[i]   =  x_i[i-1] * dlogx;
                else
                    x_i[i]   =  x_i[i-1] * dlogx2;
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

        R_core = dx[1];
        
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
            //source_pressure_prefactor_left[i]  = (surf[i-1]/vol[i] - 1./dx[i]); 
            //source_pressure_prefactor_right[i] = (surf[i]  /vol[i] - 1./dx[i]); 
            
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
        
        if(debug > 0) cout<<"Init: Finished Init. Got initial dt = "<<dt<<" This is only temporary dt_init, not used for the first timestep."<<endl;
        //
        // Matrix init via Eigen
        //

        if(debug > 0) cout<<"Init: Setting up friction workspace."<<endl;


        alphas_sample  = Eigen::VectorXd::Zero(num_cells+2);
        
        if(num_species >= 1) {
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
    
    radiation_matrix_T   = Matrix_t::Zero(num_species, num_species);
    radiation_matrix_M   = Matrix_t::Zero(num_species, num_species);
    radiation_vec_input  = Vector_t(num_species);
    radiation_vec_output = Vector_t(num_species);
    radiation_cv_vector  = Vector_t(num_species);
    radiation_T3_vector  = Vector_t(num_species);
    
    previous_monitor_J = std::vector<double>(num_bands) ;
    previous_monitor_T = std::vector<double>(num_species) ;
    
    solar_heating = Eigen::VectorXd::Zero(num_bands,  1);
    S_band        = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    dS_band       = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
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
        
        
        //if(debug > 2)
        //   cout<<" in planck_matrix, lT = "<<planck_matrix(p,0)<<" percentile = "<<planck_matrix(p,1)<<endl;
    }
    
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 1."<<endl;
    for(int s=0; s<num_species; s++) {
        species[s].dS = Eigen::VectorXd::Zero(num_cells+2,  1);
    }
    
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 2."<<endl;
    
    double templumi = 0;
    for(int b=0; b<num_bands; b++) {
        if(debug > 0) cout<<"Init: Assigning stellar luminosities. b ="<<b<<endl;
        
        if(num_bands == 1) {
            solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.)  + UV_star * 1.;
            templumi += solar_heating(b);
            
            cout<<" Solar heating is "<<solar_heating(b)<<endl;
        }
        else{
            cout<<"SOLAR HEATING in bin "<<b;
            cout<<" from/to lmin/lmax"<<l_i[b];
            cout<<"/"<<l_i[b+1]<<" with T_star = "<<T_star;
            
            solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.) * compute_planck_function_integral3(l_i[b], l_i[b+1], T_star)  + UV_star * 1.;
            templumi += solar_heating(b);
            
            cout<<" is "<<solar_heating(b)<<endl;
            
        }
    }
    cout<<"TOTAL SOLAR HEATING / Lumi = "<<templumi<<" lumi = "<<(templumi*4.*pi*rsolar*rsolar*pi)<<endl;
    
    double totallumi = 0;
    double totallumiincr = 0;
    for(int b=0; b < num_bands; b++) {
        
        totallumiincr = compute_planck_function_integral3(l_i[b], l_i[b+1], T_star);
        totallumi     += totallumiincr;
        cout<<" INIT BANDS, b = "<<b<<" l_i[b] = "<<l_i[b]<<" l_i[b+1] = "<<l_i[b+1]<<" l_i12[b] = "<<l_i12[b]<<" fraction = "<<totallumiincr<<endl;
    }
    cout<<" Total lumi / sigma T^4/pi is = "<<totallumi<<endl;
    
    total_opacity        = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    cell_optical_depth   = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    radial_optical_depth = Eigen::MatrixXd::Zero(num_cells+2, num_bands);

    T_FLD          = Eigen::MatrixXd::Zero(num_cells+2, num_species);
    T_FLD2         = Eigen::MatrixXd::Zero(num_cells+2, num_species);
    T_FLD3         = Eigen::MatrixXd::Zero(num_cells+2, num_species);
    
    Etot_corrected = Eigen::MatrixXd::Zero(num_cells+2, 1);
    
    Jrad_FLD       = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    Jrad_FLD2      = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    Jrad_FLD3      = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    Jrad_init      = Eigen::MatrixXd::Zero(num_cells+2, num_bands);
    Jrad_FLD_total = Eigen::VectorXd::Zero(num_cells+2,  1);
    Jrad_FLD_total2= Eigen::VectorXd::Zero(num_cells+2,  1);
    tridiag        = BlockTriDiagSolver<Eigen::Dynamic>(num_cells+2, num_bands + num_species) ;
    
    double rade_l = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERL", debug, 1.).value;
    double rade_r = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERR", debug, 1.).value;
    
    for(int j = 0; j < num_cells+1; j++) {
        
        if(num_bands == 1)  {
            for(int s = 0; s< num_species; s++){
                
                //Jrad_FLD(j,0) += rad_energy_multipier * sigma_rad*pow(species[s].prim[j].temperature,4) / pi;
                //cout<<" Assigning value j = "<<j<<" species = "<<s<<" with rad_energy_multiplier = "<<rad_energy_multiplier<<" and T = "<<species[s].prim[j].temperature<<endl;
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
                    
                    double SHOCK_TUBE_MID = read_parameter_from_file<double>(filename,"PARI_INIT_SHOCK_MID", debug, 0.5).value;
                    
                    if(x_i12[j] < SHOCK_TUBE_MID) 
                        Jrad_FLD(j,0) = rade_l;
                    else
                        Jrad_FLD(j,0) = rade_r;
                    
                }
                 
                
                if(debug > 1) {
                //if(debug > 1 && j==5) {
                    cout<<" Jrad("<<j<<","<<0<<") = "<<Jrad_FLD(j,0)<<" dJrad_species["<<s<<"] = "<<rad_energy_multiplier * compute_planck_function_integral3(l_i[0], l_i[1], species[s].prim[j].temperature)<<" T_rad = "<<pow(pi*Jrad_FLD(j,0)/sigma_rad,0.25)<<" at Tgas = "<<species[s].prim[j].temperature<<endl;
                }
                    
            }
        } else {
            
            for(int b = 0; b < num_bands; b++) {
                for(int s = 0; s< num_species; s++){
                    Jrad_FLD(j,b)  = rad_energy_multiplier * compute_planck_function_integral3(l_i[b], l_i[b+1], species[s].prim[j].temperature);
                    if(Jrad_FLD(j,b) < 1e-50)
                        Jrad_FLD(j,b) = 1e-50;
                    
                    Jrad_init(j,b) = Jrad_FLD(j,b);
                    //Jrad_FLD(j,b) += rad_energy_multipier * compute_planck_function_integral(l_i[b], l_i[b+1], species[s].prim[j].temperature);
                    //if(debug >= 0 && j==1)
                     if(debug > 1 && j==5) {
                        cout<<" Jrad("<<j<<","<<b<<") = "<<Jrad_FLD(j,b)<<" dJrad_species["<<s<<"] = "<<rad_energy_multiplier * compute_planck_function_integral3(l_i[b], l_i[b+1], species[s].prim[j].temperature)<<" at T ="<<species[s].prim[j].temperature<<endl; 
                     }
                        
                }
            }
            
            
        }
        
    }
    
    cout<<" IN INIT: FINISHED INIT. RETURNING TO MAIN NOW."<<endl;
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
        num_bands          = base->num_bands;
        this->workingdir   = base->workingdir;
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
        TEMPERATURE_BUMP_STRENGTH  = read_parameter_from_file<double>(filename,"TEMPERATURE_BUMP_STRENGTH", debug, 0.).value; 
        pressure_broadening_factor = read_parameter_from_file<double>(filename,"PRESSURE_BROADENING", debug, 0.).value; 
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Finished reading boundaries."<<endl;
        if(debug > 0) cout<<"         Boundaries used in species["<<speciesname<<"]: "<<boundary_left<<" / "<<boundary_right<<endl;

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
        // Determine which equation of states we are using, and assign thermodynamic variables: heat capacity, adiabatic ratio, and equation of state
        //
        
        //Gas species
        if(is_dust_like == 0) {
            cv            = 0.5 * degrees_of_freedom * Rgas / mass_amu; 
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom;
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu) ;
        }
        //Dust species
        else {
            cv            = 1.5 * Rgas / mass_amu; // 3 Degrees of freedom per atom (high T limit of debye theory)
            gamma_adiabat = (degrees_of_freedom + 2.)/ degrees_of_freedom; // Total degrees of freedom per grain
            eos           = new IdealGas_EOS(degrees_of_freedom, cv, mass_amu) ;
        }
        
        if(debug >= 0) cout<<"        Species["<<species_index<<"] got a gamma_adiabatic = "<<gamma_adiabat<<" and cv = "<<cv<<endl;

        u_analytic     = np_zeros(num_cells+2);

        prim   = std::vector<AOS_prim>(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        primlast = std::vector<AOS_prim>(num_cells+2); 
        de_e   = std::vector<double>(num_cells+2);
        prim_l = std::vector<AOS_prim>(num_cells+2);
        prim_r = std::vector<AOS_prim>(num_cells+2);
        temp_temperature = std::vector<double>(num_cells+2);
        
        timesteps    = np_zeros(num_cells+2);		    
        timesteps_cs = np_zeros(num_cells+2);	
        timesteps_rad = np_zeros(num_cells+2);
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
            
            u1l *= initial_fraction;
            u1r *= initial_fraction;
            
            u3l *= initial_fraction;
            u3r *= initial_fraction;
            
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
            
            if(debug > 0) cout<<"        Species["<<species_index<<"] initialized problem 1. with left rho/v/p/T ="<<pl.density<<"/"<<pl.speed<<"/"<<pl.pres <<"/"<<pl.temperature<<" init atmo "<<init_static_atmosphere  <<" density excess"<<density_excess<<endl;
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
                initialize_hydrostatic_atmosphere(filename);
            }
            else if (init_static_atmosphere == 2) {
                initialize_exponential_atmosphere();
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
        
        opacity                = Eigen::MatrixXd::Zero(num_cells+2, num_bands); //num_cells * num_bands
        opacity_planck         = Eigen::MatrixXd::Zero(num_cells+2, num_bands); //num_cells * num_bands
        fraction_total_opacity = Eigen::MatrixXd::Zero(num_cells+2, num_bands); //num_cells * num_bands
    
}



void c_Species::initialize_hydrostatic_atmosphere(string filename) {
    
    if(debug > 0) cout<<"            ATTENTION: Initializing hydrostatic construction for species "<<speciesname<<" and overwriting prior initial values."<<endl;
    
    double temp_rhofinal;
    double factor_inner, factor_outer;
    double T_inner;
    double T_outer;
    double metric_inner;
    //double residual;
    double metric_outer;
    
    double dphi_factor = read_parameter_from_file<double>(filename,"PARI_DPHI_FACTOR", debug, 1.).value;
    //
    // Start with the outermost cell and build up a hydrostatic atmosphere
    // Fulfilling KKM16, Eqn. 17
    //
    //long double residual  = 0.;
    
    //
    // First, initialize (adiabatic) temperature
    //
    for(int i=num_cells+1; i>=0; i--) {
            
        if(base->temperature_model == 'P')
            prim[i].temperature = - dphi_factor * base->phi[i] / (cv * gamma_adiabat) + const_T_space;
        else
            prim[i].temperature = const_T_space;
            
            //Add temperature bumps and troughs
            //prim[i].temperature += TEMPERATURE_BUMP_STRENGTH * 4.  * exp( - pow(base->x_i12[i] - 1.e-1 ,2.) / (0.1) );
            //prim[i].temperature -= TEMPERATURE_BUMP_STRENGTH * 15. * exp( - pow(base->x_i12[i] - 7.e-3 ,2.) / (1.5e-3) );
    }
    
    //At this point, the right ghost cell is already initialized, so we can just build up a hydrostatic state from there
    int iter_start;
    if(base->type_of_grid==2)
        iter_start = base->grid2_transition_i;
    else
        iter_start = num_cells;
    
    for(int i=iter_start; i>=0; i--)  {
        
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
        
        //if(i==num_cells || i==num_cells-1)
        //    metric_inner = dphi_factor * (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
        //else
        metric_inner = (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
        
        temp_rhofinal = u[i+1].u1 * (factor_outer + metric_outer) / (factor_inner - metric_inner);
        
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner) ;
        
        //cout<<"    HYDROSTATIC: i/r[i]/rho/T = "<<i<<"/"<<base->x_i[i]<<"/"<<temp_rhofinal<<"/"<<prim[i].temperature<<endl;
        //
        // Debug info
        // 
        //if( (i==20 || i== num_cells-20) && debug >= 0) {
        //if( i==300 ) {
        if(temp_rhofinal < 0) {
            
            char a;
            cout.precision(16);
            cout<<"NEGATIVE DENSITY IN INIT HYDROSTATIC i="<<i<<"/"<<num_cells<<endl;
            if((factor_inner - metric_inner) < 0) {
                
                cout<<"     Negative denominator detected. Product (gamma-1)*cv*T_inner is too small for chosen discretization."<<endl;
                
            }
            cout<<endl;
            cout<<"     dphi debug: phi["<<i+1<<"] = "<<base->phi[i+1]<<" phi["<<i<<"] = "<<base->phi[i]<<" num_cells = "<<num_cells<<endl;
            cout<<"     metric_inner debug: dPhi = "<<(base->phi[i+1] - base->phi[i])<<" dx[i+1]+dx[i] = "<<(base->dx[i+1] + base->dx[i])<<" omega*dx  = "<<base->omegaminus[i]  * base->dx[i]<<endl;
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
    
    if(base->type_of_grid == 2) {
        
        for(int i = iter_start+1; i<num_cells+1; i++) {
            double rr = pow(base->x_i[i]/base->x_i[iter_start],3.);
            u[i] = AOS(u[iter_start].u1 / rr, 0., cv * u[iter_start].u1 / rr * const_T_space) ;
            
        }
    }
    
    if(u[2].u1 > 1e40) {
        cout<<"    ---WARNING---"<<endl;
        cout<<"    IN CONSTRUCT HYDROSTATIC:"<<endl;
        cout<<"    species["<<speciesname<<"] has reached a very high initial density  of > 1e40 at the inner boundary."<<endl;
        cout<<"    The well-balanced scheme seems to be unable to compensate minor velocity fluctuations under those circumstances."<<endl;
        cout<<"    Be advised, your solution is probably about to unphysically explode."<<endl;
    }
        
    
  
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
        bondi_radius = base->planet_mass*G / (pow(csi,2.)); //Most restrictive sonic point, as we are not quite isothermal
    
        for( int i = 1; i < num_cells + 1; i++) {
            smoothed_factor = density_excess * (1./(1.+std::exp(-(base->x_i12[i] - bondi_radius)/(0.01 * bondi_radius )) ) ); 
            
            u[i].u1 = u[i].u1 * (1. + smoothed_factor);
            
            u[i].u3 = u[i].u3 * (1. + smoothed_factor);
            if(i==1)
                cout<<"In initwind: First smoothed factor = "<<smoothed_factor<<endl;
        }
        
        cout<<"            Ended hydrostatic construction for species "<<speciesname<<", last smoothed factor was "<<smoothed_factor<<" while  density_excess="<<density_excess<<" r_b="<<bondi_radius<<endl;
    }
    else
        cout<<"            Ended hydrostatic construction for species "<<speciesname<<endl;
    
    
    //TODO:Give some measure if hydrostatic construction was also successful, and give warning if not.

}

void c_Species::initialize_exponential_atmosphere() {
    
    for(int i=num_cells+1; i>=0; i--) {
            prim[i].temperature = const_T_space;
    }
    
    
    for(int i=num_cells; i>=0; i--)  {
        double temp_rhofinal = u[i].u1*std::exp(-1./550e5*(base->x_i[i] - base->x_i[num_cells]));
        
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * prim[i].temperature) ;
    }
    cout<<endl<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<endl<<endl;
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

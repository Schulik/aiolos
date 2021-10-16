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

c_Sim::c_Sim(string filename_solo, string speciesfile_solo, string workingdir, int debug, int debug_cell, int debug_steps) {

        if(debug > 0) cout<<"Init position 0."<<endl;
        
        steps = -1;
        this->debug      = debug ;
        this->debug_cell = debug_cell;
        this->debug_steps= debug_steps;
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
        
        num_bands_in     = read_parameter_from_file<int>(filename,"PARI_NUM_BANDS", debug, 1).value;
        num_bands_out    = read_parameter_from_file<int>(filename,"NUM_BANDS_OUT", debug, num_bands_in).value;
        lambda_min_in       = read_parameter_from_file<double>(filename,"PARI_LAM_MIN", debug, 1e-1).value;
        lambda_max_in       = read_parameter_from_file<double>(filename,"PARI_LAM_MAX", debug, 10.).value;
        lambda_per_decade_in= read_parameter_from_file<double>(filename,"PARI_LAM_PER_DECADE", debug, 10.).value;
        lambda_min_out       = read_parameter_from_file<double>(filename,"PARI_LAM_MIN_OUT", debug, lambda_min_in).value;
        lambda_max_out       = read_parameter_from_file<double>(filename,"PARI_LAM_MAX_OUT", debug, lambda_max_in).value;
        lambda_per_decade_out= read_parameter_from_file<double>(filename,"PARI_LAM_PER_DECADE_OUT", debug, lambda_per_decade_in).value;
        
        star_mass        =  read_parameter_from_file<double>(filename,"PARI_MSTAR", debug, 1.).value;
        star_mass        *= msolar;
        T_star           = read_parameter_from_file<double>(filename,"PARI_TSTAR", debug, 5777.).value;
        R_star           = read_parameter_from_file<double>(filename,"PARI_RSTAR", debug, 0.).value;
        UV_star          = read_parameter_from_file<double>(filename,"PARI_UVSTAR", debug, 0.).value;
        X_star           = read_parameter_from_file<double>(filename,"PARI_XSTAR", debug, 0.).value;
        Lyalpha_star     = read_parameter_from_file<double>(filename,"PARI_LYASTAR", debug, 0.).value;
        
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
        reverse_hydrostat_constrution = read_parameter_from_file<int>(filename,"REVERSE_HYDROSTAT_CONSTRUCTION", debug, 0).value;
        
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
        max_timestep_change = read_parameter_from_file<double>(filename,"MAX_TIMESTEP_CHANGE", debug, 1.1).value;
        dt_min_init         = read_parameter_from_file<double>(filename,"DT_MIN_INIT", debug, 1e-20).value;
        output_time = read_parameter_from_file<double>(filename,"PARI_TIME_OUTPUT", debug).value; 
        monitor_time = read_parameter_from_file<double>(filename,"PARI_TIME_DT", debug).value;
        CFL_break_time = read_parameter_from_file<double>(filename,"CFL_BREAK_TIME", debug, std::numeric_limits<double>::max()).value ;
        energy_epsilon = read_parameter_from_file<double>(filename,"ENERGY_EPSILON", debug, 0.01).value;
        
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
        use_tides         = read_parameter_from_file<int>(filename,"USE_TIDES", debug, 0).value;
        use_linear_gravity= read_parameter_from_file<int>(filename,"PARI_LINEAR_GRAV", debug, 0).value;
        use_rad_fluxes    = read_parameter_from_file<int>(filename,"PARI_USE_RADIATION", debug, 0).value;
        use_convective_fluxes = read_parameter_from_file<int>(filename,"USE_CONVECTION", debug, 0).value;
        convect_boundary_strength = read_parameter_from_file<double>(filename,"CONVECT_BOUNDARY_STRENGTH", debug, 1.1).value;
        
        use_collisional_heating = read_parameter_from_file<int>(filename,"PARI_USE_COLL_HEAT", debug, 1).value;
        use_drag_predictor_step = read_parameter_from_file<int>(filename, "PARI_SECONDORDER_DRAG", debug, 0).value;
        init_wind         = read_parameter_from_file<int>(filename,"PARI_INIT_WIND", debug, 0).value;
        alpha_collision   = read_parameter_from_file<double>(filename,"PARI_ALPHA_COLL", debug, 0).value;
        //init_mdot              = read_parameter_from_file<double>(filename,"PARI_MDOT", debug, -1.).value;
        init_sonic_radius = read_parameter_from_file<double>(filename,"INIT_SONIC_RADIUS", debug, -1e10).value;
        rad_energy_multiplier=read_parameter_from_file<double>(filename,"PARI_RAD_MULTIPL", debug, 1.).value;
        collision_model   = read_parameter_from_file<char>(filename,"PARI_COLL_MODEL", debug, 'C').value;
        opacity_model     = read_parameter_from_file<char>(filename,"PARI_OPACITY_MODEL", debug, 'C').value;
        if(opacity_model == 'M')
            init_malygin_opacities();
        cout<<" CHOSEN OPACITY MODEL = "<<opacity_model<<endl;
        
        no_rad_trans               = read_parameter_from_file<double>(filename,"NO_RAD_TRANS", debug, 1.).value;
        photocooling_multiplier    = read_parameter_from_file<double>(filename,"PHOTOCOOL_MULTIPLIER", debug, 1.).value;
        //radiation_solver           = read_parameter_from_file<int>(filename,"RADIATION_SOLVER", debug, 0).value;
        closed_radiative_boundaries = read_parameter_from_file<int>(filename,"PARI_CLOSED_RADIATIVE_BOUND", debug, 0).value;
        const_opacity_solar_factor = read_parameter_from_file<double>(filename,"CONSTOPA_SOLAR_FACTOR", debug, 1.).value;
        const_opacity_rosseland_factor = read_parameter_from_file<double>(filename,"CONSTOPA_ROSS_FACTOR", debug, 1.).value;
        const_opacity_planck_factor = read_parameter_from_file<double>(filename,"CONSTOPA_PLANCK_FACTOR", debug, 1.).value;
        temperature_model = read_parameter_from_file<char>(filename,"INIT_TEMPERATURE_MODEL", debug, 'P').value;
        friction_solver   = read_parameter_from_file<int>(filename,"FRICTION_SOLVER", debug, 0).value;
        do_hydrodynamics  = read_parameter_from_file<int>(filename,"DO_HYDRO", debug, 1).value;
        photochemistry_level = read_parameter_from_file<int>(filename,"PHOTOCHEM_LEVEL", debug, 0).value;
        dust_to_gas_ratio = read_parameter_from_file<double>(filename,"DUST_TO_GAS", debug, 0.).value;
        
        ion_precision         = read_parameter_from_file<double>(filename,"ION_PRECISION", debug, 1e-12).value;
        ion_heating_precision = read_parameter_from_file<double>(filename,"ION_HEATING_PRECISION", debug, 1e-12).value;
        ion_maxiter                  = read_parameter_from_file<int>(filename,"ION_MAXITER", debug, 128).value;
        ion_heating_maxiter          = read_parameter_from_file<int>(filename,"ION_HEATING_MAXITER", debug, 128).value;
        
        cout<<" ION_HEATING_PRECISION set to "<<ion_heating_precision<<" heating+maxiters = "<<ion_heating_maxiter<<endl;
        cout<<" ION_PRECISION set to "<<ion_precision<<" ion_maxiter = "<<ion_maxiter<<endl;
        
        rad_solver_max_iter = read_parameter_from_file<int>(filename,"MAX_RAD_ITER", debug, 1).value;
        bond_albedo       = read_parameter_from_file<double>(filename,"BOND_ALBEDO", debug, 0.).value; 
        planet_semimajor= read_parameter_from_file<double>(filename,"PARI_PLANET_DIST", debug, 1.).value; //in AU
        
        if(problem_number == 2)
            monitor_output_index = num_cells/2; //TODO: Replace with index of sonic radius for dominant species?
        else
            monitor_output_index = 1;
        
        if(problem_number==2) {
            
            planet_mass     = read_parameter_from_file<double>(filename,"PARI_PLANET_MASS", debug, 1).value; //in Earth masses
            planet_mass     *= mearth;
            
            planet_position = read_parameter_from_file<double>(filename,"PARI_PLANET_POS", debug, 0.).value;  //inside the simulation domain
            rs              = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_LENGTH", debug, 0.).value; //Gravitational smoothing length in hill 
            rs_time         = read_parameter_from_file<double>(filename,"PARI_SMOOTHING_TIME", debug, 0.).value; //Time until we reach rs starting at rs_at_moment
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
        total_pressure    = np_zeros(num_cells+2);
        
        if(num_bands_in < 1) {
            cout<<" In INIT RADIATION, invalid num_bands_in = "<<num_bands_in<<" changing to num_bands_in = 1."<<endl;
            num_bands_in = 1;
        }
        if(num_bands_out < 1) {
            cout<<" In INIT RADIATION, invalid num_bands_out = "<<num_bands_out<<" changing to num_bands_out = 1."<<endl;
            num_bands_out = 1;
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
        photon_energies = np_somevalue(num_he_bands, 20. * ev_to_K * kb);
        for(int b=0; b<num_he_bands; b++) {
            photon_energies[b] = 1.24/( l_i_in[b + 1] ) * ev_to_K * kb ;
            //photon_energies[b] = 20. * ev_to_K * kb ;
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
            /*for(int i=1; i<= num_cells; i++) {
                
                if(x_i[i] > grid2_transition) {
                    grid2_transition_i = i;
                    cout<<" found grid2_transition_i = "<<i<<endl;
                    break;
                }
            }*/
            cout<<" grid_transition_i = "<<grid2_transition_i<<endl;
            
            for(int i=1; i<= num_cells; i++) {
                
                if(i <= grid2_transition_i) {
                    
                    templogx = dlogx;
                } else if (i >= (grid2_transition_i + 5)) {
                    templogx = dlogx2;
                } else {
                    templogx = dlogx + double(i-grid2_transition_i)/5. * (dlogx2 - dlogx);
                    //cout<<" in grid2generation, i/dlogx = "<<i<<"/"<<templogx<<" dlogx/dlogx2 = "<<dlogx<<"/"<<dlogx2<<endl;
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
        //TODO: Expand this with a few more sensible tests, NaN checks etc.
        
        int all_grid_ok = 1, broken_index;
        double broken_value;
        for(int i = 0; i <= num_cells+1; i++) {
            if( dx[i] < 0 ) {
                all_grid_ok = 0;
                broken_index = i;
                broken_value = dx[i];
            }
         
            //cout<<i<<" dx = "<<dx[i]/x_i12[i]<<endl;
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

        init_T_temp   = read_parameter_from_file<double>(filename,"INIT_T_TEMP", debug, 1.).value;

        if(debug > 0) cout<<"Init: About to setup species with num_species = "<<num_species<<endl;
        
        species.reserve(num_species);
        
        for(int s = 0; s < num_species; s++) {
            species.push_back(c_Species(this, filename, speciesfile, s, debug)); // Here we call the c_Species constructor
        }
        cout<<endl;
        //
        //
        // After successful init, compute all primitive quantities and determine first timestep
        //
        //
        
        for(int s = 0; s < num_species; s++)
            species[s].compute_pressure(species[s].u);
        
        dt = t_max ; // Old value needed by get_cfl_timestep()
        dt = get_cfl_timestep();
        
        if(debug > 0) cout<<"Init: Finished Init. Got initial dt = "<<dt<<" This is only temporary dt_init, not used for the first timestep."<<endl;
        //
        // Matrix init via Eigen
        //

        if(debug > 0) cout<<"Init: Setting up friction workspace."<<endl;


        alphas_sample  = Eigen::VectorXd::Zero(num_cells+2);
        friction_sample  = Eigen::VectorXd::Zero(num_cells+2);
        
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
        // Ion species feel each other strongly: Multiply ionic friction strength by the 'canonical' value of 100.
        //
        //
        for(int si=0; si < num_species; si++)
            for(int sj=0; sj < num_species; sj++) {
                if(species[si].is_dust_like == 1 && species[sj].is_dust_like == 1) {
                    friction_coeff_mask(si,sj) = 0.;
                    friction_coeff_mask(sj,si) = 0.;
                }
                if(std::abs(species[si].static_charge) > 0 && std::abs(species[sj].static_charge) > 0 && (si != sj)) {
                    cout<<" Setting friction mask to ionic value for species si/sj = "<<si<<"/"<<sj<<endl;
                    friction_coeff_mask(si,sj) = 10000.;
                    friction_coeff_mask(sj,si) = 10000.;
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
    S_band        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    dS_band       = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
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
        
        //if(debug > 2)
        //   cout<<" in planck_matrix, lT = "<<planck_matrix(p,0)<<" percentile = "<<planck_matrix(p,1)<<endl;
    }
    
//     debug = 4;
//     char inpit;
//     cout<<" Attempting to debug extreme temperatures:"<<endl;
//     cout<<compute_planck_function_integral3(l_i[0], l_i[1], 5777);
//     cin>>inpit;
//     double temp2 = 1e-40;
//     for(int i = 0; i<80; i++) {
//             cout<<"Energy fraction in global band at T = "<<temp2<<": "<<compute_planck_function_integral3(l_i[0], l_i[1], temp2)<<" planck2 "<<compute_planck_function_integral2(l_i[0], l_i[1], temp2)/(sigma_rad * pow(temp2,4.) / pi)<<endl;
//             temp2 *= 10.;
//     }
    
    //cout<<"Energy fraction in global band at T = 1e-40: "<<compute_planck_function_integral3(l_i[0], l_i[1], 1e-40)<<endl;
    //cin>>inpit;
    
    
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 1."<<endl;
    for(int s=0; s<num_species; s++) {
        species[s].dS = Eigen::VectorXd::Zero(num_cells+2,  1);
        species[s].dG = Eigen::VectorXd::Zero(num_cells+2,  1);
    }
    
    if(debug > 0) cout<<"Init: Assigning stellar luminosities 2."<<endl;
    
    double templumi = 0;
    for(int b=0; b<num_bands_in; b++) {
        if(debug > 0) cout<<"Init: Assigning stellar luminosities. b ="<<b<<endl;
        
        if(num_bands_in == 1) {
            solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.);
            solar_heating(b) += UV_star/(4.*pi*pow(planet_semimajor*au,2.));
            solar_heating(b) += X_star/(4.*pi*pow(planet_semimajor*au,2.));
            templumi += solar_heating(b);
            
            cout<<" Solar heating is "<<solar_heating(b)<<" T_star = "<<T_star<<" pow(R_star*rsolar,2.) "<<pow(R_star*rsolar,2.)<<" pow(planet_semimajor*au,2.) "<<planet_semimajor<<endl;
        }
        else{
            cout<<"SOLAR HEATING in bin "<<b;
            cout<<" from/to lmin/lmax"<<l_i_in[b];
            cout<<"/"<<l_i_in[b+1]<<" with T_star = "<<T_star;
            
            solar_heating(b)  = sigma_rad * pow(T_star,4.) * pow(R_star*rsolar,2.)/pow(planet_semimajor*au,2.) * compute_planck_function_integral4(l_i_in[b], l_i_in[b+1], T_star);
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
        }
        
        cout<<" is F = "<<solar_heating(b)<<endl;
            
    }
    
    cout<<"TOTAL SOLAR HEATING / Lumi = "<<templumi<<" lumi = "<<(templumi*4.*pi*rsolar*rsolar*pi)<<endl;
    
    double totallumi = 0;
    double totallumiincr = 0;
    for(int b=0; b < num_bands_in; b++) {
        
        totallumiincr = compute_planck_function_integral3(l_i_in[b], l_i_in[b+1], T_star);
        totallumi     += totallumiincr;
        if(debug > 0 ) cout<<" INIT BANDS, b = "<<b<<" l_i_in[b] = "<<l_i_in[b]<<" l_i[b+1] = "<<l_i_in[b+1]<<" l_i12[b] = "<<l_i12_in[b]<<" fraction = "<<totallumiincr<<" fraction for T=1430 K = "<<compute_planck_function_integral3(l_i_in[b], l_i_in[b+1], 1430.)<<endl;
    }
    cout<<" Total lumi / sigma T^4/pi is = "<<totallumi<<endl;
    
    total_opacity        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    cell_optical_depth   = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    radial_optical_depth = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);

    total_opacity_twotemp        = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    cell_optical_depth_twotemp   = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    radial_optical_depth_twotemp = Eigen::MatrixXd::Zero(num_cells+2, num_bands_in);
    
    Etot_corrected = Eigen::MatrixXd::Zero(num_cells+2, 1);
    
    Jrad_FLD       = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    Jrad_init      = Eigen::MatrixXd::Zero(num_cells+2, num_bands_out);
    tridiag        = BlockTriDiagSolver<Eigen::Dynamic>(num_cells+2, num_bands_out + num_species) ;
    
    double rade_l = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERL", debug, 1.).value;
    double rade_r = read_parameter_from_file<double>(filename,"PARI_RADSHOCK_ERR", debug, 1.).value;
    init_J_factor = read_parameter_from_file<double>(filename,"INIT_J_FACTOR", debug, -1.).value;
    
    
    for(int j = 0; j < num_cells+1; j++) {
        
        if(num_bands_out == 1)  {
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
        
    }
    
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

c_Species::c_Species(c_Sim *base_simulation, string filename, string species_filename, int species_index, int debug) {
    
        base               = base_simulation;
        this_species_index = species_index;
        num_cells          = base->num_cells;
        num_bands_in       = base->num_bands_in;
        num_bands_out      = base->num_bands_out;
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
        TEMPERATURE_BUMP_STRENGTH    = read_parameter_from_file<double>(filename,"TEMPERATURE_BUMP_STRENGTH", debug, 0.).value; 
        pressure_broadening_factor   = read_parameter_from_file<double>(filename,"PRESSURE_BROADENING", debug, 0.).value; 
        pressure_broadening_exponent = read_parameter_from_file<double>(filename,"BROADENING_EXP", debug, 1.).value; 
        
        if(debug > 0) cout<<"        Species["<<species_index<<"] Init: Finished reading boundaries."<<endl;
        if(debug > 0) cout<<"         Boundaries used in species["<<speciesname<<"]: "<<boundary_left<<" / "<<boundary_right<<endl;

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
        
        timesteps     = np_zeros(num_cells+2);		    
        timesteps_cs  = np_zeros(num_cells+2);	
        timesteps_rad = np_zeros(num_cells+2);
        finalstep     = np_zeros(num_cells+2);
        timesteps_de  = np_zeros(num_cells+2);
        
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
                const_rho_scale = read_parameter_from_file<double>(filename,"CONST_RHO_SCALE", debug, 550.e5).value;
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
        
        
        USE_WAVE = read_parameter_from_file<int>(filename,"WAVE_USE", debug, 0).value; 
        if(USE_WAVE==1) {
            WAVE_AMPLITUDE = read_parameter_from_file<double>(filename,"WAVE_AMPLITUDE", debug).value; 
            WAVE_PERIOD    = read_parameter_from_file<double>(filename,"WAVE_PERIOD", debug).value; 
        }
        
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
        
    if(base->photochemistry_level > 0 && std::abs(static_charge) > 0 )
        base->density_floor = read_parameter_from_file<double>(filename,"ION_FLOOR", debug, 1.e-20).value;
    else
        base->density_floor = read_parameter_from_file<double>(filename,"DENSITY_FLOOR", debug, 1.e-20).value;
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
            
            eos->get_p_over_rho_analytic(&T_outer, &factor_outer);
            eos->get_p_over_rho_analytic(&T_inner, &factor_inner);
            
            metric_outer = (base->phi[i+1] - base->phi[i]) * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
            metric_inner = (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
            
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
        
    } else {
            int negdens = 0; 
            
            for(int i=1; i<=iter_start; i++)  {
            
            //
            // Construct next density as to fulfil the hydrostatic condition
            //
            T_outer = prim[i+1].temperature ;
            T_inner = prim[i].temperature ;
            
            eos->get_p_over_rho_analytic(&T_outer, &factor_outer);
            eos->get_p_over_rho_analytic(&T_inner, &factor_inner);
            
            metric_outer = (base->phi[i+1] - base->phi[i]) * base->omegaplus[i+1] * base->dx[i+1] / (base->dx[i+1] + base->dx[i]);
            metric_inner = (base->phi[i+1] - base->phi[i]) * base->omegaminus[i]  * base->dx[i]   / (base->dx[i+1] + base->dx[i]);
            
            temp_rhofinal = u[i].u1 *  (factor_inner - metric_inner)/(factor_outer + metric_outer);
            
            //////////////////////////////////////////////////////////////////////////////////////
            double floor = base->density_floor * mass_amu ;
            if(temp_rhofinal < floor) {
                
                if(temp_rhofinal < 0.)
                    negdens = 1;
                
                temp_rhofinal = floor;
            }
            //////////////////////////////////////////////////////////////////////////////////////
            
            //if(this_species_index == 2) {
            //    cout<<" electrons, i/rho = "<<i<<"/"<<temp_rhofinal<<endl;
            //}
            
            u[i+1] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner);
            
            if(negdens == 2) {
                
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
        
        
    }
    
    if(base->type_of_grid == -2) {
        
        for(int i = iter_start+1; i<num_cells+1; i++) {
            double rr = pow(base->x_i[i]/base->x_i[iter_start],5.);
            u[i] = AOS(u[iter_start].u1 / rr, 0., cv * u[iter_start].u1 / rr * prim[i].temperature) ;
        }
    }
    
    
    //makeshift density fix
    //if(this_species_index == 0)
        //cout<<"FIXING SMALL DENSITIES to density floor ="<<base->density_floor*mass_amu<<endl;
    
    for(int i = 0; i<num_cells+1; i++) {
            double floor = base->density_floor * mass_amu ;
            if(u[i].u1 < floor) {
                //cout<<"FIXING SMALL DENSITIES in i ="<<i<<" for species = "<<this_species_index<<endl;
                u[i] = AOS(floor, 0., cv * floor * prim[i].temperature) ;
            }
    }
    
    if(base->init_T_temp > 2.7) {
        
            if(this_species_index == 0)
                cout<<" In INIT: Overwriting T with user-defined temperature:"<<base->init_T_temp;
                
            for(int i=num_cells; i>=0; i--)  {
                
                //if(x_i12[i]<2.e10) {
                    
                    double temprho = u[i].u1;
                
                    u[i] = AOS(temprho, 0., cv * temprho * base->init_T_temp) ;
                //}
                
                //prim[i].temperature = base->init_T_temp;
            }
        }
        else
            cout<<" In INIT:NO TOverwriting T with user-defined temperature:"<<base->init_T_temp;
    
    
    if(u[2].u1 > 1e40) {
        cout<<"    ---WARNING---"<<endl;
        cout<<"    IN CONSTRUCT HYDROSTATIC:"<<endl;
        cout<<"    species["<<speciesname<<"] has reached a very high initial density  of > 1e40 at the inner boundary. Density = "<<u[2].u1<<endl;
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
        
        cout<<"@@@ Ended hydrostatic construction for species "<<speciesname<<", last smoothed factor was "<<smoothed_factor<<" while  density_excess="<<density_excess<<" r_b="<<bondi_radius<<endl<<endl;
    }
    else
        cout<<"@@@ Ended hydrostatic construction for species "<<speciesname<<endl<<endl;
    
    //TODO:Give some measure if hydrostatic construction was also successful, and give warning if not.

}

void c_Species::initialize_exponential_atmosphere() {
    
    //cout<<" temperatures : ";
    
    for(int i=num_cells+1; i>=0; i--) {
            prim[i].temperature = const_T_space;
            //cout<<prim[i].temperature<<" ";
    }
    
    //cout<<endl<<" scale_H : "<<const_rho_scale<<endl;
    //cout<<" densities : ";
    
    for(int i=num_cells; i>=0; i--)  {
        double temp_rhofinal = u[i].u1*std::exp(-1./const_rho_scale*(base->x_i[i] - base->x_i[num_cells]));
        
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * prim[i].temperature) ;
        
        //cout<<temp_rhofinal<<" ";
    }
    cout<<endl<<endl;
    cout<<"            Ended expoenential density construction for species "<<speciesname<<endl;
    cout<<endl<<endl;
}



void c_Species::initialize_hydrostatic_atmosphere_iter(string filename) {
    
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
            for (int i=0; i < num_ghosts; i++) {
                
                double dens_wall;  
                
                if(base->problem_number == 1)
                    dens_wall = SHOCK_TUBE_UL.u1 *  mass_amu;
                else {
                    //cout<<" in fixed boundaries for s = "<<this_species_index<<", u1 = "<<BACKGROUND_U.u1<<" dens = "<<BACKGROUND_U.u1 *  mass_amu<<endl;
                    dens_wall = BACKGROUND_U.u1 *  mass_amu;
                }
                AOS cons_fixed = AOS(dens_wall, 0., cv * dens_wall * prim[i].temperature);
                eos->compute_primitive(&u[i],&prim[i], 1) ;
                //cout<<" found p = "<<prim[i].pres<<endl;
            }
            
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
                
                double r =base->x_i12[i];
                double mdot      = prim.density*prim.speed*r*r;
                double freefallv = -std::sqrt(2.*G*base->planet_mass/r);
                if(mdot < -1e18)
                    prim.density = -1e18/freefallv/r/r;
                    
                prim.pres = prim.pres;// -  prim.density * dphi ;
                //prim.pres = prim.pres -  prim.density * dphi ;
                prim.pres = std::max( prim.pres, 0.0) ;
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

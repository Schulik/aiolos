#include "main.h"

//extern AOS* init_AOS(int num);

/*
Initializer for a new simulation. Stores all parameters in the object, given to by the wrapper.
    
        par_control: List of control parameters
    
        par_numerics: List of numerical parameters
    
        par_physics: List of physics parameters
*/
hydro_run::hydro_run(string filename) {
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //  Numerical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double dx0;
        dx0         = read_parameter_from_file(filename,"PARI_DOMAIN_DX", TYPE_DOUBLE, debug).dvalue;
        domain_min  = read_parameter_from_file(filename,"PARI_DOMAIN_MIN", TYPE_DOUBLE, debug).dvalue;
        domain_max  = read_parameter_from_file(filename,"PARI_DOMAIN_MAX", TYPE_DOUBLE, debug).dvalue;
        num_cells = (int)((domain_max - domain_min)/dx0);
        //num_cells must always be read in pretty early, as all other memory allocations in this function depend on it.
        //Note-to-self: Don't change the order of read_parameter_from_file commands.
        
        dt          = read_parameter_from_file(filename,"PARI_TIME_DT", TYPE_DOUBLE, debug).dvalue;
        cflfactor   = read_parameter_from_file(filename,"PARI_CFLFACTOR", TYPE_DOUBLE, debug).dvalue;
        t_max       = read_parameter_from_file(filename,"PARI_TIME_TMAX", TYPE_DOUBLE, debug).dvalue;
        output_time = read_parameter_from_file(filename,"PARI_TIME_OUTPUT", TYPE_DOUBLE, debug).dvalue; 
        globalTime = 0.0;    
    
    
    
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //    Control parameters for users
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //In main, execution is about to start.

        simname = filename;
        //par_control  = par_control;
        //par_numerics = par_numerics;
        //par_physics  = par_physics;
        
        //debug             = read_parameter_from_file(filename,"PARI_DEBUGLEVEL", TYPE_INT, debug).ivalue;   
        
        //cout<<"Pos1"<<endl;
        boundaries_number = read_parameter_from_file(filename,"PARI_BOUND_TYPE", TYPE_INT, debug).ivalue;
        //cout<<"Pos2"<<endl;
        problem_number    = read_parameter_from_file(filename,"PARI_PROBLEM_NUMBER", TYPE_INT, debug).ivalue;
        //cout<<"Pos3"<<endl;
        
        use_self_gravity  = read_parameter_from_file(filename,"PARI_SELF_GRAV_SWITCH", TYPE_INT, debug).ivalue;
        use_linear_gravity= read_parameter_from_file(filename,"PARI_LINEAR_GRAV", TYPE_INT, debug).ivalue;
        use_rad_fluxes    = read_parameter_from_file(filename,"PARI_USE_RADIATION", TYPE_INT, debug).ivalue;
        
        double constopa;
        if(use_rad_fluxes==1)
            constopa   = read_parameter_from_file(filename,"PARI_CONST_OPAC", TYPE_DOUBLE, debug).dvalue;
        else
            constopa   = 1.;
        
        
        opacity        = np_somevalue(num_cells, constopa);
        opticaldepth   = np_zeros(num_cells);
        radiative_flux = np_zeros(num_cells);
        temperature    = np_somevalue(num_cells, 1.);
            
        
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
        // Physical
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        gamma_adiabat   = read_parameter_from_file(filename,"PARI_GAMMA", TYPE_DOUBLE, debug).dvalue; //ratio of specific heats
        ggminusone      = gamma_adiabat*(gamma_adiabat-1.);
        const_T         = read_parameter_from_file(filename,"PARI_CONST_TEMP", TYPE_DOUBLE, debug).dvalue;
        T_increment     = read_parameter_from_file(filename,"PARI_TEMP_INCR", TYPE_DOUBLE, debug).dvalue;
        
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Simulation data: Grid and variables
        //
        ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        x_i        = np_zeros(num_cells+1);		//The cell boundaries
        x_i12      = np_zeros(num_cells);		    //The cell mid positions
        timesteps  = np_zeros(num_cells);		    
        timesteps_cs= np_zeros(num_cells);	
        finalstep   = np_zeros(num_cells);
        dx         = np_zeros(num_cells);
        
        //Compute cell wall boundaries
        for(int i=0; i< num_cells+1; i++) {
            x_i[i] = (domain_min + dx0 * double(i));
        }
            

        //Compute cell mid positions
        for(int i=0; i<num_cells; i++) {
            x_i12[i] = (0.5 * (x_i[i+1] + x_i[i]) );
            dx[i]  = dx0;
        }
            

        u      = init_AOS(num_cells); //{ np_ones(num_cells*3);   //Conserved hyperbolic variables: density, mass flux, energy density
        oldu   = init_AOS(num_cells);   //Conserved hyperbolic variables: density, mass flux, energy density
        phi    = np_zeros(num_cells);  //Parabolic Variables: gravitational potential
        enclosed_mass = np_zeros(num_cells);
        source = init_AOS(num_cells);  //Parabolic Variables: gravitational potential
        flux   = init_AOS(num_cells+1);
        left_ghost = AOS(0,0,0);
        right_ghost= AOS(0,0,0);
        
        //u_output = [];    //Array of arrays to store snapshots of u
        //phi_output = [];  //Array of arrays to store snapshots of phi

        timecount = 0;
        plotcounter = 0;
        
        if(debug > 0) {
                cout<<"DEBUGLEVEL 1: Samples of generated values"<<endl;
                //cout<<"rho[0], v[0], e[0] = "<<u[0]<<" "<<u[num_cells]<<" "<<u[2*num_cells]<<endl;
                cout<<"First cell coordinates: |<--"<<x_i[0]<<" // "<<x_i12[0]<<" // "<<x_i[1]<<"-->|"<<endl;
                cout<<"Last cell coordinates:  |<--"<<x_i[num_cells-1]<<" // "<<x_i12[num_cells-1]<<" // "<<x_i[num_cells]<<"-->|"<<endl;
                cout<<"Value of gamma: "<<gamma_adiabat<<endl;
        } 
        
        //
        // Problem 1: two states, left and right, separated at a mid-position
        //
        if(problem_number == 1) {
            
            double u1l = read_parameter_from_file(filename,"PARI_INIT_DATA_U1L", TYPE_DOUBLE, debug).dvalue;
            double u2l = read_parameter_from_file(filename,"PARI_INIT_DATA_U2L", TYPE_DOUBLE, debug).dvalue;
            double u3l = read_parameter_from_file(filename,"PARI_INIT_DATA_U3L", TYPE_DOUBLE, debug).dvalue;
            
            double u1r = read_parameter_from_file(filename,"PARI_INIT_DATA_U1R", TYPE_DOUBLE, debug).dvalue;
            double u2r = read_parameter_from_file(filename,"PARI_INIT_DATA_U2R", TYPE_DOUBLE, debug).dvalue;
            double u3r = read_parameter_from_file(filename,"PARI_INIT_DATA_U3R", TYPE_DOUBLE, debug).dvalue;
            
            SHOCK_TUBE_MID = read_parameter_from_file(filename,"PARI_INIT_SHOCK_MID", TYPE_DOUBLE, debug).dvalue;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            SHOCK_TUBE_UL = AOS(u1l, u1l*u2l, 0.5*u1l*u2l*u2l + u3l/(gamma_adiabat-1.) );
            SHOCK_TUBE_UR = AOS(u1r, u1r*u2r, 0.5*u1r*u2r*u2r + u3r/(gamma_adiabat-1.));
            
            initialize_shock_tube_test(SHOCK_TUBE_UL, SHOCK_TUBE_UR);
            
            if(debug > 0) 
                cout<<"Successfully initialized problem 1."<<endl;
        }
        //
        // Problem 2: A constant background state, with a planet embedded into it, free to evolve as physics dictates
        //
        else if(problem_number == 2) {
            
            cout<<"Problem 2 pos0."<<endl;
            
            double u1 = read_parameter_from_file(filename,"PARI_INIT_DATA_U1", TYPE_DOUBLE, debug).dvalue;
            double u2 = read_parameter_from_file(filename,"PARI_INIT_DATA_U2", TYPE_DOUBLE, debug).dvalue;
            double u3 = read_parameter_from_file(filename,"PARI_INIT_DATA_U3", TYPE_DOUBLE, debug).dvalue;
            
            cout<<"Problem 2 pos1."<<endl;
            
            planet_mass     = read_parameter_from_file(filename,"PARI_PLANET_MASS", TYPE_DOUBLE, debug).dvalue; //in Earth masses
            planet_position = read_parameter_from_file(filename,"PARI_PLANET_POS", TYPE_DOUBLE, debug).dvalue;  //inside the simulation domain
            rs              = read_parameter_from_file(filename,"PARI_SMOOTHING_LENGTH", TYPE_DOUBLE, debug).dvalue; //Gravitational smoothing length in hill 
            rs_time         = read_parameter_from_file(filename,"PARI_SMOOTHING_TIME", TYPE_DOUBLE, debug).dvalue; //Time until we reach rs starting at rs_at_moment
            rs_at_moment    = 0.2;
            init_static_atmosphere = read_parameter_from_file(filename,"PARI_INIT_STATIC", TYPE_INT, debug).ivalue; //Yesno, isothermal at the moment
            
            //mdot boundaries
            if(boundaries_number == 3) {
                mdot = read_parameter_from_file(filename,"ACCRETION_RATE", TYPE_DOUBLE, debug).dvalue; //Accretion rate in numerical units
            }
            
            cout<<"Problem 2 pos2."<<endl;
            
            //Conversion from shock tube parameters (given as dens, velocity, pressure) to conserved variables (dens, momentum, internal energy)
            BACKGROUND_U = AOS(u1, u1*u2, 0.5*u1*u2*u2 + u3/(gamma_adiabat-1.) );
            
            cout<<"Problem 2 pos3."<<endl;
            
            initialize_background(BACKGROUND_U);
            init_grav_pot();
            
            if(debug > 0) 
                cout<<"Successfully initialized problem 2."<<endl;
        }
        else 
            initialize_default_test();
        
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // Last thing to do: Initialize derived quantities
        //
        //////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        pressure        = np_zeros(num_cells+2); //Those helper quantities are also defined on the ghost cells, so they get +2
        internal_energy = np_zeros(num_cells+2);
        speed           = np_zeros(num_cells+2);
        cs              = np_zeros(num_cells+2);
        
        //
        // If we want to initialize a hydrostatic atmosphere, then we overwrite the so far given density/pressure data and set every velocity to 0
        //
        if(init_static_atmosphere == 1) {
            initialize_hydrostatic_atmosphere();
        }
        
        if(debug > 0)
            cout<<"Ended init."<<endl;
}

//
//
//
void hydro_run::initialize_hydrostatic_atmosphere() {
    
    cout<<"ATTENTION: Initializing hydrostatic atmosphere and overwriting prior initial values."<<endl;
    
    //const_T = 1.; //read_parameter_from_file(filename,"PARI_CONST_TEMP", TYPE_DOUBLE, debug).dvalue; //Time until we reach rs starting at;
    cv      = 1.;
    
    long double temp_rhofinal;
    long double factor_inner, factor_outer;
    long double factor_grav;
    long double delta_phi;
    long double T_inner;
    long double T_outer;
    
    // debugvariables
    long double debug_dens;
    long double debug_energy;
    
    //
    // Start with the outermost cell and build up a hydrostatic atmosphere
    // Fulfilling KKM16, Eqn. 17
    //
    long double dcount = 3.;
    long double temp_increase = T_increment;
    long double tempgrad  = 0.;
    long double tempgrad2 = 0.;
    long double residual  = 0.;
    
    //
    // First, initialize (adiabatic) temperature
    //
    for(int i=0; i<num_cells; i++) {
            //temperature[i] = planet_mass / x_i12[i] / (cv * gamma_adiabat) + 1.;
            temperature[i] = - phi[i] / (cv * gamma_adiabat) + 1.;
        }
        
    
    //Last normal cell has to be awkwardly initialized
    u[num_cells-1] = AOS(u[num_cells-2].u1, 0., cv * u[num_cells-2].u1 * temperature[num_cells-2]);
    
    for(int i=num_cells-2; i>=0; i--)  {
        
        //
        // Increase temperatures for arbitrary stratifications
        //
        //T_outer = 100./(x_i12[i]   + 1.); //const_T + temp_increase * dcount ;
        //T_inner = 100./(x_i12[i+1] + 1.); //T_outer + temp_increase;
        
        T_outer = temperature[i+1]; //const_T + temp_increase * (domain_max - x_i12[i+1] + 0.1); //const_T + temp_increase * dcount ;
        T_inner = temperature[i];   //const_T + temp_increase * (domain_max - x_i12[1] + 0.1); //T_outer + temp_increase;
        
        dcount += 1.;
        
        /*T_inner = const_T / sqrt(dcount + 1.);
        T_outer = const_T / sqrt(dcount);
        dcount += 1.;
        
        if(T_inner < 1.)
            T_inner = 1.;
        if(T_outer < 1.)
            T_outer = 1;
        */
        
        //
        // Construct next density as to fulfil the hydrostatic condition
        //
        factor_outer = (gamma_adiabat-1.) * cv * T_outer; //TODO: Replace with T_outer for non-isothermal EOS
        factor_inner = (gamma_adiabat-1.) * cv * T_inner; //TODO: Replace with T_inner for non-isothermal EOS
        delta_phi    = phi[i+1] - phi[i];
        factor_grav  = - 0.5 * (x_i12[i+1] - x_i12[i]) * delta_phi;
        
        // TODO: Update with delta-x-dependent prefactors
        //temp_rhofinal = u[i+1].u1 * (1. + 0.5 * delta_phi / factor_outer) / (1. - 0.5 * delta_phi / factor_outer);
        //temp_rhofinal = u[i+1].u1 * (2.* factor_outer / delta_phi + 1.) / (2. * factor_inner / delta_phi - 1.);
        //temp_rhofinal = u[i+1].u1 * (2.* factor_outer + delta_phi) / (2. * factor_inner - delta_phi);
        temp_rhofinal = u[i+1].u1 * (2.* factor_outer + delta_phi) / (2. * factor_inner - delta_phi);
        
        //temp_rhofinal = u[i+1].u1 * (2. * factor_inner - delta_phi) / (2.* factor_outer + delta_phi);
        
        
        //temp_rhofinal = u[i+1].u1 * 1. / (1. - 0.5 * delta_phi / factor_outer); // Sub-hydrostatic envelope: accretion commences, resulting density is
                                                                                  // approx hydrostatic * 1e-2
        
        //
        // 1st order assignment
        // 
        u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner) ;
        
        //
        // 2nd order start
        //
        //residual       = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
        //temp_rhofinal += residual / (cv * (gamma_adiabat-1.) * T_inner);
        //u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner) ;
        
        //residual      = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
        //temp_rhofinal += residual / (cv * (gamma_adiabat-1.) * T_inner);
        //u[i] = AOS(temp_rhofinal, 0., cv * temp_rhofinal * T_inner) ;
        
        if( (i==20 || i== num_cells-20) && debug >= 0) {
            cout.precision(16);
            cout<<" i="<<i<<endl;
            cout<<"In hydostatic init: factor_outer/delta_phi = "<< factor_outer / delta_phi << " factor_inner/delta_phi = "<< factor_inner / delta_phi <<endl;
            cout<<"In hydostatic init: factor_dens = "<< (2.* factor_outer / delta_phi + 1.) / (2. * factor_inner / delta_phi - 1.) <<endl;
            cout<<"Ratio of densities inner/outer = "<< temp_rhofinal/u[i+1].u1 <<endl;
            cout<<"Ratio of temperatures inner/outer = "<<T_inner/T_outer<<endl;
            cout<<"Ratio of pressures inner/outer = "<<cv * temp_rhofinal * T_inner /u[i+1].u3<<endl;
            tempgrad = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            tempgrad2 = ((gamma_adiabat-1.)*cv*T_outer + 0.5 * delta_phi) * u[i+1].u1 - ((gamma_adiabat-1.)*cv*T_inner  + 0.5 * delta_phi ) * u[i].u1 ;
            cout<<"pressure diff "<<(gamma_adiabat-1.)*u[i+1].u3 - (gamma_adiabat-1.)*(u[i].u3)<<endl;
            cout<<"density diff "<<(u[i].u1 - u[i+1].u1)<<endl;
            cout<<"density sum " <<(u[i].u1 + u[i+1].u1)<<endl;
            cout<<"density sum with potential "<<0.5 * (u[i].u1 + u[i+1].u1) * delta_phi<<endl;
            cout<<"residual = "<<residual<<endl;
            //residual = (gamma_adiabat-1.)*cv*u[i+1].u1*T_outer - (gamma_adiabat-1.)*cv*u[i].u1*T_inner  + 0.5 * (u[i].u1 + u[i+1].u1) * delta_phi;
            //cout<<"residual2 = "<<residual<<endl;
            //cout<<"sum of hydrostatic gradients = "<<tempgrad<<endl;
            //cout<<"sum2 of hydrostatic gradients = "<<tempgrad2<<endl;
        }
    }
}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void hydro_run::initialize_default_test() {
    
    cout<<"Initialized default simulation."<<endl;
    
    for(int i=0; i<num_cells; i++) 
        u[i] = AOS(1,1,1);
}


//
// Initial conditions for shock tubes, dividing the domain into a left constant value and another right constant value
//
void hydro_run::initialize_shock_tube_test(const AOS &leftval,const AOS &rightval) {
    
    cout<<"Initialized shock tube test."<<endl;
    
    for(int i=0; i<num_cells; i++) {
            
        if(x_i[i] < SHOCK_TUBE_MID) 
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
    
    for(int i=0; i<num_cells; i++) 
        u[i] = background;
    
    left_ghost  = background;
    right_ghost = background;
}

//
// Non-changing boundary conditions. Give values as leftval and rightval (for example the initial values from SHOCK_TUBE_UL and UR)
//                                   and they are assigned to the ghost zones
//
//                                  In the case a planet is involved, the planet is by default assumed to be sitting at x=0 and 
//                                  then the left boundary is WALL and the right is CONST.
void hydro_run::boundaries_const_both(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost ){
    
    if(debug >= 3)
        cout<<"Inside boundaries_const_both. Assigning values left = "<<leftval.u1<<" "<<leftval.u2<<" "<<leftval.u3<<" and right = "<<leftval.u1<<" "<<rightval.u2<<" "<<rightval.u3<<endl;
    
    //For the shock tube test we use the initial values as boundaries
    if(problem_number == 1) {
        left_ghost  = SHOCK_TUBE_UL;
        right_ghost = SHOCK_TUBE_UR;
    } 
    else if(problem_number == 2) {
        
        left_ghost = AOS( leftval.u1, -leftval.u2, leftval.u3); 
        right_ghost= BACKGROUND_U;
        
    } else {
        left_ghost  = leftval;
        right_ghost = rightval;
    }
}




//
// Open boundaries on both ends: Takes the last and second last values on each side and extrapolates them to generate the ghost values
//
//                                  In the case a planet is involved, the planet is by default assumed to be sitting at x=0 and 
//                                  then the left boundary is WALL and the right is OPEN.
void hydro_run::boundaries_open_both(AOS &left_ghost, const AOS &leftval, const AOS &leftval2, const AOS &rightval2, const AOS &rightval, AOS &right_ghost ){
    
    if(debug >= 3)
        cout<<"Inside boundaries_open_both. Assigning values left = "<<leftval.u1<<" "<<leftval.u2<<" "<<leftval.u3<<" and right = "<<leftval.u1<<" "<<rightval.u2<<" "<<rightval.u3<<endl;
    
    double dl1 = (leftval2.u1 - leftval.u1); /// (x_i[1] - x_i[0])
    double dl2 = (leftval2.u2 - leftval.u2); /// (x_i[1] - x_i[0])
    double dl3 = (leftval2.u3 - leftval.u3); /// (x_i[1] - x_i[0])
    double dr1 = (rightval.u1 - rightval2.u1); /// (x_i[cell_number] - x_i[cell_number-1])
    double dr2 = (rightval.u2 - rightval2.u2); // (x_i[cell_number] - x_i[cell_number-1])
    double dr3 = (rightval.u3 - rightval2.u3); // (x_i[cell_number] - x_i[cell_number-1])

    
    if(problem_number == 2) {
        
        //double small_momentum = 1e-10;
        //double small_dens     = 1e-10;
        double tmp_r_dens     = rightval.u1 + dr1;
        double tmp_r_momentum = rightval.u2 + dr2;
        
        //Wall left
        left_ghost = AOS( leftval.u1, -leftval.u2, leftval.u3); 
        
        //Inflow/outflow boundary conditions, that do not allow sign change across the last cell boundary, otherwise horrible things happen
        
        //Sign change detect
        if((rightval.u2 + dr2)/rightval.u2 < 0) {
            tmp_r_momentum = rightval.u2;
            //tmp_r_momentum = small_momentum;
        }
        if(tmp_r_dens < 0.) {
            tmp_r_dens = rightval.u1;
            
        }
        double press_boundary   = (gamma_adiabat-1.) * (rightval.u3 - 0.5 * rightval.u2 * rightval.u2 / rightval.u1 );
        double hydrostat_energy = (press_boundary - rightval.u1 * (phi[num_cells-1] - phi_right_ghost) ) / (gamma_adiabat-1.);
        
        
        //right_ghost= AOS(tmp_r_dens, tmp_r_momentum, rightval.u3 + dr3 );
        //right_ghost = AOS(rightval.u1, rightval.u2, rightval.u3);
        right_ghost = AOS(rightval.u1, rightval.u2, hydrostat_energy);
    }
    else {
        
        left_ghost  = AOS( leftval.u1 - dl1, leftval.u2  - dl2,  leftval.u3 - dl3 ); 
        right_ghost = AOS(rightval.u1 + dr1, rightval.u2 + dr2, rightval.u3 + dr3 );
        
    }
        
}

void hydro_run::boundaries_planet_mdot(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost ) {
    
    if(debug >= 3)
        cout<<"Inside boundaries_planet_mdot. Assigning values left = "<<leftval.u1<<" "<<leftval.u2<<" "<<leftval.u3<<" and right = "<<leftval.u1<<" "<<rightval.u2<<" "<<rightval.u3<<endl;
    
    left_ghost  = AOS( leftval.u1, -leftval.u2, leftval.u3); 
    
    //right_ghost = BACKGROUND_U; //This is identical to const background boundaries, but we will develop this further with time-dependent factors etc.

    double rho0 = BACKGROUND_U.u1;
    double e0   = BACKGROUND_U.u3 - 0.5 * BACKGROUND_U.u2 * BACKGROUND_U.u2 / rho0;
    
    right_ghost = AOS(rho0, -mdot, e0 + 0.5 * mdot * mdot / rho0 );
    //we compute back and forth between the two momenta and E's because the initial conditions could be 0 momentum, then we would have mdot=0    
    
}

//
// Wall boundary conditions on both. 
//                                   

void hydro_run::boundaries_wall_both(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost ){
    
    if(debug >= 3)
        cout<<"Inside boundaries_wall_both. Assigning values left = "<<leftval.u1<<" "<<leftval.u2<<" "<<leftval.u3<<" and right = "<<leftval.u1<<" "<<rightval.u2<<" "<<rightval.u3<<endl;
    
    left_ghost = AOS( leftval.u1, -leftval.u2, leftval.u3); 
    right_ghost= AOS( rightval.u1, -rightval.u2, rightval.u3);
}


hydro_run::~hydro_run() {
    
    delete x_i;
    delete x_i12;
    delete u;
    delete oldu;
    delete phi; 
    
    cout<<"The destructor is called, but its empty so far."<<endl;
    
}

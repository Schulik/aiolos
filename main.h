#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <exception>
#include <string>
#include <vector>
#include <math.h>
#include "advection.h"
#include "source.h"

using namespace std;

//
// Definition of helper constants, mainly used to initialize simulation parameter from file
//
const int TYPE_INT = 0;
const int TYPE_DOUBLE = 1;
const int TYPE_STRING = 2;

//Basic physics quantities
const double G        = 6.678e-8; //cgs units
const double pi       = 3.141592;
const double navo     = 6e23; // particles per mole
const double sigma    = 5.67e-5;   //erg cm-2 s-1 K-4
const double kb       = 1.38e-16;  //erg/K
const double km       = 1e5; //kilometers in cm
const double mearth   = 5.98e27;  //g
const double msolar   = 2.0e33;   //g
const double au       = 1.49e13;  //cm
const double year     = 365*24*3600;
const double myear    = 1e6*year;
const double gyear    = 1e9*year;
const double rearth   = 6370e5;
const double rjupiter = 74000*km;
const double rsolar   = 695510*km;
const double pc       = 3.08567758e18; //cm
const double kpc      = 1e3 * pc;
const double mpc      = 1e6 * pc;

// For entropy computation, to be set to sensible parameters e.g. at setting of initial conditions
const double P_ref   = 1e-20;
const double Rho_ref = 1e-20; 

const double cflfactor = 0.9;

//
// Functions mimicking certain numpy functionalities
//
double *np_zeros(int);
double *np_ones(int);
double *np_somevalue(int, double);

//
// Everything to define and read simulation parameters
//
struct simulation_parameter {
    string name;
    int type;
    double dvalue;
    int ivalue;
    string svalue;
};

//
// Array of strings structure, for comfortable use of 3-vectors
//
struct AOS {
 double u1;
 double u2;
 double u3;
 
    AOS(double u1=0, double u2=0, double u3=0) 
        : u1(u1), u2(u2), u3(u3)
    {
    }
    
    // addop. doesn't modify object. therefore const.
    AOS operator+(const AOS& a) const
    {
        return AOS(a.u1+u1, a.u2+u2, a.u3+u3);
    }
    
    AOS operator-(const AOS& a) const
    {
        return AOS(u1-a.u1, u2-a.u2, u3-a.u3);
    }
    
    AOS operator*(const double& a) const
    {
        return AOS(u1*a, u2*a, u3*a);
    }
    AOS operator/(const double& a) const
    {
        return AOS(u1/a, u2/a, u3/a);
    }
 
}; 



class hydro_run
{
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Control variables
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    string simname;
    int boundaries_number;
    int problem_number;
    int debug;
    int solver;
    int use_self_gravity;
    int use_linear_gravity;
    int use_rad_fluxes;
    int suppress_warnings;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Numerical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    double *dx;
    double domain_min, domain_max;
    int num_cells;
    double *x_i;    		//The cell boundaries
    double *x_i12;  		//The cell mid positions

    double dt;
    double cflfactor;
    double t_max;
    double globalTime;
    double *timesteps;
    double *timesteps_cs;
    double *finalstep;
    double output_time;
    double snd_crs_time;
    int plotskip;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    double planet_mass;     //in Earth masses
    double planet_position; //inside the simulation domain
    double rs;
    double rs_at_moment;
    double rs_time;
    int    init_static_atmosphere;
    int    static_atmosphere_tempprofile;
    double gamma_adiabat;           //ratio of specific heats
    double ggminusone;
    double mdot;
    double cv;
    double const_T;
    double T_increment;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Simulation data: Grid and variables
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    AOS *u;              //Conserved hyperbolic variables: density, mass flux, energy density
    AOS *oldu;           //Conserved hyperbolic variables: density, mass flux, energy density
    AOS *source;
    AOS *flux;
    double *phi;            //Parabolic Variables: gravitational potential
    double phi_left_ghost, phi_right_ghost;
    double ghost_xi_left, ghost_xi_right;
    AOS left_ghost, right_ghost; //Ghost cells for boundaries

    AOS *u_output;       //Array of arrays to store snapshots of u
    double *phi_output;     //Array of arrays to store snapshots of phi
    
    int timecount;
    int plotcounter;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // temp variables: all helper variables, like intercell fluxes that are stored in between funciton calls etc.
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    double *cs;         //Speed of sound array
    double *pressure;   
    double *internal_energy;
    double *speed;
    double *enclosed_mass;
    double *opticaldepth;
    double *opacity;
    double *temperature;
    double *radiative_flux;
    
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Variables for different scenarios
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //Shock tube tests
    double SHOCK_TUBE_MID = 0.3;
    AOS SHOCK_TUBE_UL, SHOCK_TUBE_UR; //For problem_number==1
    AOS BACKGROUND_U;                 //For problem_number==2
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Functions
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //
    // Init
    //
    
    void initialize_shock_tube_test(const AOS &left,const AOS &right); //For Problem 1
    void initialize_background(const AOS &state); //For problem 2
    void initialize_default_test();
    void initialize_space_test(AOS background);
    void initialize_custom_setup();
    void initialize_hydrostatic_atmosphere();
    void initialize_gravitational_potential();
    
    
    //
    // Utilities
    //
    double get_phi_grav(); //Returns the value of the gravitational potential at a position
    double get_max_soundspeed();
    double get_cfl_timestep();
    void compute_pressure();
    void print_AOS_component_tofile(double *x, AOS* data, AOS* fluxes , int timestepnumber);
    
    
    //
    // Boundaries
    //
    //AOS get_boundaries_1();
    //AOS get_boundaries_2();
    //AOS get_boundaries_3();
    //AOS get_boundaries_4();
    
    void boundaries_const_both(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost );
    void boundaries_open_both(AOS &left_ghost, const AOS &leftval, const AOS &leftval2, const AOS &rightval2, const AOS &rightval, AOS &right_ghost );
    void boundaries_planet_mdot(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost );
    void boundaries_wall_both(AOS &left_ghost, const AOS &leftval, const AOS &rightval, AOS &right_ghost );

    //
    // Source terms
    //
    
    //Gravity
    void init_grav_pot();
    double get_p_hydrostatic(AOS &u, double &phi_l, double &phi_r, const int &i);
    double get_phi_grav(double &r, double &mass);
    void update_mass_and_pot();
    AOS source_grav(AOS &u, double &phi_l, double &phi_r);
    
    //Radiation
    void update_radiation();
    AOS source_radflux(int i);
    
    // 
    // Riemann solver
    //
    double get_wave_speeds(AOS &input_vec, const double &sign, const int &j);
    AOS analytic_flux(AOS &input_vec, const int &j);
    
    AOS hll_flux (AOS &leftvalue, AOS &rightvalue,                               const int &jleft, const int &jright);
    AOS hllc_flux(AOS &leftvalue, AOS &rightvalue, double &phi_l, double &phi_r, const int &jleft, const int &jright);
    
public:
    
    hydro_run(string);
    ~hydro_run();
    
    void execute(); //Main loop
    void set_debug(int);
    void set_suppress_warnings(int j) {suppress_warnings = j;}
};

AOS* init_AOS(int num);
vector<string> stringsplit(const string& str, const string& delim);
simulation_parameter read_parameter_from_file(string, string, int, int);

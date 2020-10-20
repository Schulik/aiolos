#ifndef _AIOLOS_MAIN_H_
#define _AIOLOS_MAIN_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <new>
#include <exception>
#include <string>
#include <vector>
#include <math.h>
#include <type_traits>
#include "advection.h"
#include "source.h"

#include "enum.h"
#include "eos.h"

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


////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Functions
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//
// Functions mimicking certain numpy functionalities
//
std::vector<double> np_zeros(int);
std::vector<double> np_ones(int);
std::vector<double> np_somevalue(int, double);

//
// Everything to define and read simulation parameters
//
template<typename T>
struct simulation_parameter {
    string name;
    T value ;
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
    
    AOS& operator+=(const AOS& a) {
        u1 += a.u1 ;
        u2 += a.u2 ;
        u3 += a.u3 ;
        return *this ;
    }

    AOS& operator-=(const AOS& a) {
        u1 -= a.u1 ;
        u2 -= a.u2 ;
        u3 -= a.u3 ;
        return *this ;
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

struct AOS_prim {
    AOS_prim(double d=0, double s=0, double p=0, double cs=0, double u=0, double T=0)
        : density(d), speed(s), pres(p), 
          sound_speed(cs), internal_energy(u), temperature(T)
    { } ;
    double density, speed, pres ;
    double sound_speed, internal_energy, temperature;
} ;


std::vector<AOS> init_AOS(int num);
vector<string> stringsplit(const string& str, const string& delim);

template<typename T>
simulation_parameter<T> read_parameter_from_file(string, string, int);

template<typename T>
simulation_parameter<T> read_parameter_from_file(string, string, int, T);



class c_Species;


////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  CLASS SIMULATION
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



class c_Sim
{
public:   
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Control variables
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    string simname;
    
    int num_species;
    std::vector<c_Species> species;
    
    int problem_number;
    int debug;
    int use_self_gravity;
    int use_linear_gravity;
    int use_rad_fluxes;
    int suppress_warnings;
    int init_geometry;
    
    Geometry geometry ;
    

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Numerical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    std::vector<double>dx;
    double domain_min, domain_max;
    int num_cells;
    double cells_per_decade;
    int type_of_grid;       //0 = cartesian, 1 = log
    std::vector<double> x_i;    		//The cell boundaries
    std::vector<double> x_i12;  		//The cell mid positions
    std::vector<double> x_iVC;         // Volumetric centres
    std::vector<double> surf;
    std::vector<double> vol;
    std::vector<double> omegaplus;
    std::vector<double> omegaminus;
    std::vector<double> source_pressure_prefactor_left;
    std::vector<double> source_pressure_prefactor_right;

    double dt;
    double cflfactor;
    double t_max;
    double globalTime;
    double output_time;
    double max_snd_crs_time;
    int steps;
    int timecount;

    IntegrationType order ;
    int num_ghosts ;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physical, global
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    double planet_mass;     //in Earth masses
    double planet_position; //inside the simulation domain
    double rs;
    double rs_at_moment = 0 ;
    double rs_time;
    
    std::vector<double> phi;            //Parabolic Variables: gravitational potential
    std::vector<double> enclosed_mass;
    

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Variables for different scenarios
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    double mdot;
    double T_increment;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Functions
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //
    // Init
    //
    
    void initialize_gravitational_potential();
    
    //
    // Utilities
    //
    double get_max_soundspeed();
    double get_cfl_timestep();
    
    //
    // Source terms
    //
    
    //Gravity
    void init_grav_pot();
    void update_mass_and_pot();
    double get_phi_grav(double &r, double &mass);
    
    
public:
    
    c_Sim() {};
    c_Sim(string parameter_filename, string species_filename, int debug=0);
    ~c_Sim();
    
    void execute(); //Main loop
    void compute_friction_step(); 
    void set_debug(int);
    void set_suppress_warnings(int j) {suppress_warnings = j;}
    
    int get_num_cells() {return num_cells; }
};


////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  CLASS SPECIES
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class c_Species
{
public:
    
    c_Sim *base;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Numerical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    string name;
    int this_species_index;
    int boundaries_number;
    BoundaryType boundary_left;
    BoundaryType boundary_right;
    int num_cells;
    int debug;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physics variables
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    double mass_amu;
    double initial_fraction;
    
    std::vector<AOS> u;              // Conserved hyperbolic variables: density, mass flux, energy density
    std::vector<AOS> dudt[2];        // Time derivative of u at each stage ;      
    std::vector<AOS> source;         // Gravitational source term
    std::vector<AOS> source_pressure;// Geometric source term
    std::vector<AOS> flux;

    
    double const_T_space;
    double const_opacity;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Primitives and other helper quantities
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<double> timesteps;
    std::vector<double> timesteps_cs;
    std::vector<double> finalstep;
    double snd_crs_time;

    std::vector<AOS_prim> prim ;
    std::vector<AOS_prim> prim_l ; // Reconstructed left/ right edges
    std::vector<AOS_prim> prim_r ;
    std::vector<double> opticaldepth;
    std::vector<double> opacity;
    std::vector<double> radiative_flux;
    
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Init and scenarios
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int  init_static_atmosphere;
    
    //Wave tests
    int    USE_WAVE;
    double WAVE_AMPLITUDE;
    double WAVE_PERIOD;
    
    //Shock tube tests
    double SHOCK_TUBE_MID = 0.3;
    AOS SHOCK_TUBE_UL, SHOCK_TUBE_UR; //For problem_number==1
    AOS BACKGROUND_U;                 //For problem_number==2
    
    //Initial temperature profile
    double TEMPERATURE_BUMP_STRENGTH = 1.;
    
    void initialize_shock_tube_test(const AOS &left,const AOS &right); //For Problem 1
    void initialize_background(const AOS &state); //For problem 2
    void initialize_default_test();
    void initialize_space_test(AOS background);
    void initialize_custom_setup();
    void initialize_hydrostatic_atmosphere();
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physics computations
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // 
    // Riemann solver and source functions
    //
    double get_p_hydrostatic(AOS &u, double &phi_l, double &phi_r, const int &i);
    double get_dp_hydrostatic_nonuniform(const int &i, const int &plusminus);
    void reconstruct_edge_states() ;
    
    AOS hllc_flux(int);
    AOS source_grav(AOS &u, int &j);

    void compute_pressure(std::vector<AOS>& u) {
        eos->compute_primitive(&(u[0]), &(prim[0]), num_cells+2) ;    
        eos->compute_auxillary(&(prim[0]), num_cells+2);
    }
    
    //
    // Radiation
    //
    int use_rad_fluxes;
    void update_radiation(std::vector<AOS>& u);
    AOS source_radflux(int i);
    
    //
    // Equations of state
    //
    double cv, gamma_adiabat ;
    EOS_Base* eos ;
    
    //
    // Boundaries
    //
    void add_wave(std::vector<AOS>& u, double time);
    
    void apply_boundary_left(std::vector<AOS>& u);
    void apply_boundary_right(std::vector<AOS>& u);
    void user_boundary_left(std::vector<AOS>& u);
    void user_boundary_right(std::vector<AOS>& u);

    void user_initial_conditions();
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // File operations etc.
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int suppress_warnings;
    void print_AOS_component_tofile(int timestepnumber);
    
    void read_species_data(string filename, int species_index);
    void set_debug(int);
    void set_suppress_warnings(int j) {suppress_warnings = j;}    
    
    //c_Species();
    c_Species(c_Sim *base_simulation, string filename, string species_filename, int species_index, int debug=0);
    ~c_Species();
    void execute(std::vector<AOS>& u, std::vector<AOS>& dudt);
    
};

#endif//_AIOLOS_MAIN_H

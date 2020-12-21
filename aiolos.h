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
#include <gsl/gsl_sf_lambert.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include "block_tridiag_solve.h"

#include "enum.h"
#include "eos.h"

#ifndef NUM_SPECIES
#define NUM_SPECIES Eigen::Dynamic
#endif


using namespace std;

//Basic physics quantities
const double G        = 6.678e-8; //cgs units
const double pi       = 3.141592653589793;
const double c_light   = 2.99792458e10;
const double navo     = 6.02214e23; // particles per mole
const double amu      = 1.66054e-24; //g
const double h_planck = 6.62607015e-27;//  cm^2 g/s
const double Rgas     = 8.31446261815324e7;  //erg/K/g
//const double Rgas_fake = 1.;
const double kb       = 1.380649e-16;  //erg/K
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
const double angstroem= 1e-4; //cm
const double ergcm2_to_wattperm2 = 1e-3;
const double sigma_rad = 5.670374419e-5;   //erg cm-2 s-1 K-4
const double sigma_rad2= 2*h_planck*c_light*c_light/pow(angstroem,4.);

// For entropy computation, to be set to sensible parameters e.g. at setting of initial conditions
const double P_ref   = 1e-20;
const double Rho_ref = 1e-20; 

const double cflfactor = 0.9;
const int debug2 = 0;

////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  Functions
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double delta_ij(int i, int j);

//
// Functions mimicking certain numpy functionalities
//
std::vector<double> np_zeros(int);
std::vector<double> np_ones(int);
std::vector<double> np_somevalue(int, double);

double compute_planck_function_integral(double lmin, double lmax, double temperature);
double compute_planck_function_integral2(double lmin, double lmax, double temperature);

//
// Everything to define and read simulation parameters
//
template<typename T>
struct simulation_parameter {
    string name;
    T value ;
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
    AOS_prim(double d=0, double s=0, double p=0, double cs=0, double u=0, double T=0, double n=0)
        : density(d), speed(s), pres(p), 
          sound_speed(cs), internal_energy(u), temperature(T), number_density(n)
    { } ;
    double density, speed, pres ;
    double sound_speed, internal_energy, temperature, number_density;
} ;

struct Monitored_Quantities {
    
    Monitored_Quantities(double u1=0, double u2=0, double u3=0)
        : u1_tot(u1), u2_tot(u2), u3_tot(u3)
    { } ;
    double u1_tot, u2_tot, u3_tot;
} ;

std::vector<AOS> init_AOS(int num);
vector<string> stringsplit(const string& str, const string& delim);

template<typename T>
simulation_parameter<T> read_parameter_from_file(string, string, int);

template<typename T>
simulation_parameter<T> read_parameter_from_file(string, string, int, T);

void update_cons_prim_after_friction(AOS* cons, AOS_prim* prim, double dEkin, double newv, double mass, double gamma, double cv);


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
    string workingdir;
    
    int problem_number;
    int debug;
    int use_self_gravity;
    int use_linear_gravity;
    int use_rad_fluxes;
    int suppress_warnings;
    int do_hydrodynamics;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Numerical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    int num_species;
    std::vector<c_Species> species;
    
    int num_bands;
    int num_plancks = 512;
    int i_wien;
    int i_rayleighjeans;
    double lT_spacing;
    double T_crit_min;
    double T_crit_max; //Critical min/max temperatures beyond which one should extend the wavelength grid
    
    int num_cells;
    std::vector<double>dx;
    double domain_min, domain_max;
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

    double lminglobal = 1e-6, lmaxglobal = 1e14; //in mum; corresponds to Wien-temperature maxima from 3e9 to 1e-7 K
    double lambda_min, lambda_max, lambda_per_decade;
    std::vector<double> l_i;
    std::vector<double> l_i12;
    
    int init_geometry;
    char collision_model;
    char opacity_model;
    
    Geometry geometry ;
    
    double dt;
    double cflfactor;
    double t_max;
    double timestep_rad2;
    double globalTime;
    double output_time;
    double monitor_time;
    double max_snd_crs_time;
    double rad_energy_multiplier;
    int steps;
    int timecount;
    int monitor_output_index;

    IntegrationType order ;
    int num_ghosts ;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physical, global
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //
    // Conserved quantities to be monitored etc.
    //
    
    Monitored_Quantities monitored_quantities;
    Monitored_Quantities initial_monitored_quantities;
    
    //
    // Gravitation
    //
    
    double star_mass;
    
    double planet_mass;     //in Earth masses
    double planet_position; //inside the simulation domain
    double planet_semimajor;
    double rs;
    double rs_at_moment = 0 ;
    double rs_time;
    
    double scale_rb;
    double scale_rh;
    double scale_cs;
    double scale_vk;
    double scale_time;
    
    std::vector<double> phi;
    std::vector<double> enclosed_mass;
    
    //
    // Friction
    //
    
    int friction_solver;
    
    using Matrix_t = Eigen::Matrix<double, NUM_SPECIES,NUM_SPECIES, Eigen::RowMajor>;
    using Vector_t = Eigen::Matrix<double, NUM_SPECIES, 1>;
    Matrix_t friction_matrix_T;
    Matrix_t identity_matrix;
    Matrix_t friction_coefficients;
    Matrix_t friction_coeff_mask;
    Vector_t friction_vec_input;
    Vector_t friction_vec_output;
    Vector_t friction_dEkin;
    Vector_t unity_vector;
    Eigen::PartialPivLU<Matrix_t> LU;
    Vector_t dens_vector;
    Vector_t numdens_vector;
    Vector_t mass_vector;

    Matrix_t radiation_matrix_T;
    Matrix_t radiation_matrix_M;
    Vector_t radiation_vec_input;
    Vector_t radiation_vec_output;
    Vector_t radiation_cv_vector;
    Vector_t radiation_T3_vector;
    
    Eigen::VectorXd alphas_sample;
    
    //
    // Radiation
    //
    double T_star;
    double UV_star;
    double R_star;
    
    int radiation_matter_equilibrium_test; //If set to 1, sets J = J_init in update_radiation()
    int radiation_diffusion_test_linear;
    int radiation_diffusion_test_nonlinear;
    double CFL_break_time; //Numerical time after which cflfactor=0.9. Used in get_cfl_timestep()
    
    std::vector<double> previous_monitor_J;
    std::vector<double> previous_monitor_T;
    
    //Eigen::MatrixXd F_up;      //num_cells * num_bands each
    //Eigen::MatrixXd F_down;
    //Eigen::MatrixXd F_plus;
    //Eigen::MatrixXd F_minus;
    //Eigen::MatrixXd dJrad;
    //Eigen::MatrixXd S_total;
    Eigen::MatrixXd S_band;
    Eigen::MatrixXd dS_band;
    Eigen::MatrixXd solar_heating;
    
    Eigen::MatrixXd total_opacity;
    Eigen::MatrixXd cell_optical_depth;
    Eigen::MatrixXd radial_optical_depth;

    Eigen::MatrixXd T_FLD ;
    Eigen::MatrixXd T_FLD2 ;
    Eigen::MatrixXd T_FLD3 ;
    
    Eigen::MatrixXd planck_matrix;
    
    Eigen::MatrixXd Etot_corrected ;
    
    Eigen::MatrixXd Jrad_FLD ;

    BlockTriDiagSolver<Eigen::Dynamic> tridiag ;
    Eigen::MatrixXd Jrad_FLD2 ;
    Eigen::MatrixXd Jrad_FLD3 ;
    Eigen::MatrixXd Jrad_init;
    Eigen::MatrixXd Jrad_FLD_total ;
    Eigen::MatrixXd Jrad_FLD_total2 ;
    
    int rad_solver_max_iter = 1;
    double epsilon_rad_min = 1e-1;  // Some convergence measure for the radiation solver, if needed
    
    void fill_alpha_basis_arrays(int j);
    void fill_rad_basis_arrays(int, double, Eigen::MatrixXd &, Eigen::MatrixXd &);
    void compute_alpha_matrix(int j, int actually_compute_beta);
    
    void transport_radiation();     //  Called if use_rad_fluxes == 1
    void update_opacities();
    void update_fluxes(double timestep);           //  Called from transport_radiation#   
    void update_fluxes_FLD();           //  Called from transport_radiation#
    void update_fluxes_FLD2(double, Eigen::MatrixXd &,Eigen::MatrixXd &,Eigen::MatrixXd &);           //  Called from transport_radiation#
    void update_temperatures(double, Eigen::MatrixXd &,Eigen::MatrixXd &,Eigen::MatrixXd &);
    void correct_totalenergy(double, Eigen::MatrixXd &);
    //AOS source_radflux(int i);
    double compute_planck_function_integral3(double lmin, double lmax, double temperature);
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Variables for different scenarios
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int init_wind;
    double mdot;
    double T_increment;
    double alpha_collision;
    
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
    void compute_friction_analytical(); 
    void compute_friction_numerical(); 
    
    //Gravity
    void init_grav_pot();
    void update_mass_and_pot();
    double get_phi_grav(double &r, double &mass);
    
    void print_monitor(int i);
    
 
    
public:
    
    c_Sim() {};
    c_Sim(string parameter_filename, string species_filename, string workingdir, int debug=0);
    ~c_Sim();
    
    void execute(); //Main loop
    
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

    string speciesname;
    string workingdir;
    
    int this_species_index;
    int boundaries_number;
    BoundaryType boundary_left;
    BoundaryType boundary_right;
    int num_cells;
    int num_bands;
    int debug;
    
    string opacity_data_string;
    int num_opacity_datas;

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
    std::vector<double> u_analytic;
    
    Eigen::MatrixXd opacity_data;   //
    Eigen::MatrixXd opacity;        //num_cells * num_bands
    Eigen::MatrixXd opacity_planck; //num_cells
    Eigen::MatrixXd fraction_total_opacity; //num_cells * num_bands
    Eigen::MatrixXd dS;
    
    //std::vector<double> opticaldepth;
    //std::vector<double> opacity;
    //std::vector<double> radiative_flux;
    
    
    double density_excess;
    double bondi_radius;
    
    double const_T_space;
    double const_opacity;
    int    is_dust_like;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Primitives and other helper quantities
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<double> timesteps;
    std::vector<double> timesteps_cs;
    std::vector<double> timesteps_rad;
    std::vector<double> finalstep;
    double snd_crs_time;

    std::vector<AOS_prim> prim ;
    std::vector<AOS_prim> primlast ;
    std::vector<double> de_e;
    std::vector<AOS_prim> prim_l ; // Reconstructed left/ right edges
    std::vector<AOS_prim> prim_r ;
    std::vector<double> temp_temperature;
    
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
    void initialize_sound_wave();
    void compute_analytic_solution();
    
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physics computations
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // 
    // Riemann solver and source functions
    //
    double get_dp_hydrostatic(const int &i, const int &plusminus) const {
        
        double dp_final;
    
        if(plusminus == -1) {
                dp_final = - base->dx[i] * base->omegaplus[i] * (base->phi[i-1] - base->phi[i]) / (base->dx[i-1] + base->dx[i]);
        } else {
                dp_final = - base->dx[i] * base->omegaminus[i] * (base->phi[i+1] - base->phi[i]) / (base->dx[i+1] + base->dx[i]);
        }
            
        return dp_final;
    }

    void reconstruct_edge_states() ;
    
    AOS hllc_flux(int);
    AOS dust_flux(int);
    AOS source_grav(AOS &u, int &j);

    void compute_pressure(std::vector<AOS>& u) {
        eos->compute_primitive(&(u[0]), &(prim[0]), num_cells+2) ;    
        eos->compute_auxillary(&(prim[0]), num_cells+2);
    }
    
    void update_opacities();
    
    //
    // Equations of state
    //
    double cv;            //Specific heat capacity. Has to be computed from molar heat capacity (from degrees of freedom) and molecular mass
    double gamma_adiabat; 
    double degrees_of_freedom;
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

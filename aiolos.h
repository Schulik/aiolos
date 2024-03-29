#ifndef _AIOLOS_MAIN_H_
#define _AIOLOS_MAIN_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <new>
#include <exception>
#include <functional>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>
#include <ctime>
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

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif

using namespace std;

using Matrix_t = Eigen::Matrix<double, NUM_SPECIES,NUM_SPECIES, Eigen::RowMajor>;
using Vector_t = Eigen::Matrix<double, NUM_SPECIES, 1>;


//Basic physics quantities
const double G        = 6.678e-8; //cgs units
const double pi       = 3.141592653589793;
const double c_light   = 2.99792458e10;
const double navo     = 6.02214e23; // particles per mole
const double amu      = 1.66054e-24; //g
const double h_planck = 6.62607015e-27;//  cm^2 g/s
const double Rgas     = 8.31446261815324e7;  //erg/K/g    or //erg/K/mole?
//const double Rgas_fake = 1.;
const double kb       = 1.380649e-16;  //erg/K
const double km       = 1e5; //kilometers in cm
const double mearth   = 5.98e27;  //g
const double msolar   = 2.0e33;   //g
const double au       = 1.495978707e13;  //cm
const double year     = 365*24*3600;
const double myear    = 1e6*year;
const double gyear    = 1e9*year;
const double rearth   = 6370e5;
const double rjupiter = 69911*km;
const double rsolar   = 695510*km;
const double pc       = 3.08567758e18; //cm
const double kpc      = 1e3 * pc;
const double mpc      = 1e6 * pc;
const double angstroem= 1e-4; //cm
const double ergcm2_to_wattperm2 = 1e-3;
const double sigma_rad = 5.670374419e-5;   //erg cm-2 s-1 K-4
const double sigma_rad2= 2*h_planck*c_light*c_light/pow(angstroem,4.);
const double K_to_eV    = 8.621738e-5;
const double ev_to_K    = 1./K_to_eV;
const double de_broglie_e = 2.*pi* amu/2000.*kb/h_planck/h_planck;
const double de_broglie_p = 2.*pi* amu * kb/h_planck/h_planck;
const double elm_charge   = 4.80320425e-10; //statcoulomb

// For entropy computation, to be set to sensible parameters e.g. at setting of initial conditions

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
double lint(double xa, int N, double* X, double* RI);
double logint(double xa, int N, double* X, double* RI);
const float log10ff = std::log(10.);

inline double dfdx(const function<double(double)>& f, double x0, double dx) {
    
    return (f(x0+dx)-f(x0-dx))/(dx+dx);
}

inline float __int_as_float (int32_t a) { float r; memcpy (&r, &a, sizeof r); return r;} 
inline int32_t __float_as_int (float a) { int32_t r; memcpy (&r, &a, sizeof r); return r;}

//from user njuffa https://stackoverflow.com/questions/39821367/very-fast-approximate-logarithm-natural-log-function-in-c

/* natural log on [0x1.f7a5ecp-127, 0x1.fffffep127]. Maximum relative error 9.4529e-5 */
inline float njuffas_logf (float a) 
{
    float m, r, s, t, i, f;
    int32_t e;

    e = (__float_as_int (a) - 0x3f2aaaab) & 0xff800000;
    m = __int_as_float (__float_as_int (a) - e);
    i = (float)e * 1.19209290e-7f; // 0x1.0p-23
    /* m in [2/3, 4/3] */
    f = m - 1.0f;
    s = f * f;
    /* Compute log1p(f) for f in [-1/3, 1/3] */
    r = fmaf (0.230836749f, f, -0.279208571f); // 0x1.d8c0f0p-3, -0x1.1de8dap-2
    t = fmaf (0.331826031f, f, -0.498910338f); // 0x1.53ca34p-2, -0x1.fee25ap-2
    r = fmaf (r, s, t);
    r = fmaf (r, s, f);
    r = fmaf (i, 0.693147182f, r); // 0x1.62e430p-1 // log(2) 
    return r;
}
inline float njuffas_log10f(float a) {    
    return njuffas_logf(a)/log10ff;
} 

inline float fastpow1(float base, int exp) {

    if( exp == 0)
       return 1;
    float temp = fastpow1(base, exp/2);       
    if (exp%2 == 0)
        return temp*temp;
    else {
        if(exp > 0)
            return base*temp*temp;
        else
            return (temp*temp)/base; //negative exponent computation 
    }

} 

inline
double fastexp2(double x) {
  x = 1.0 + x / 1024;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}

inline double fastexpm1(double x) {
  return x + 0.5*x*x + 0.33333333333*x*x*x;
}

inline double fastexpm1_2(double x) {
  return x + 0.5*x*x;
}


//
// Functions mimicking certain numpy functionalities
//
std::vector<double> np_zeros(int);
std::vector<double> np_ones(int);
std::vector<double> np_somevalue(int, double);

std::vector<int> inp_zeros(int);
std::vector<int> inp_ones(int);
std::vector<int> inp_somevalue(int, int);


//double compute_planck_function_integral(double lmin, double lmax, double temperature);
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
class c_Sim;

/**
 * Thermochemistry reactions
 *
 * Reactions of the form A + B + ... ->  A' + B' + ... (only unidirectional!), implementing 
 * terms of the form dn_A/dt = - k * n_A * n_B ..., dn_B/dt = - k * n_A * n_B ..., dn_A'/dt = + k * n_A * n_B ...
 * If one wants chemical equilibrium to occur, the reverse reaction has to be added manually
 */
class c_reaction
{    
public:
    c_Sim *base;
    
    int reaction_number = -1;
    std::vector<int> educts;   //A list of indices
    std::vector<int> products; //A list of indices
    std::vector<int> e_stoch_i;
    std::vector<int> p_stoch_i;
    std::vector<double> e_stoch;
    std::vector<double> p_stoch;
    double delta_stoch;             //Stochiometric exponent difference between reactant and product side
    double products_total_mass;     //Total mass of products
    
    int e_num = 0.;      //Number reactants/educts (educt = reactant in German)
    int p_num = 0.;      //Number products
    
    double r = 0.;       //Reaction coefficient (a constant for now)
    double reac_a = 0.;  //Reaction coefficients of the form a * T^b * exp(-c/T)
    double reac_b = 0.;
    double reac_c = 0.;
    bool is_reverse_reac;    //Will try to automatically find the reverse reaction by computing r * exp(-dG/T), where the Gibbs energy difference dG has to be in the chemistry lookup tables
    bool mtype_reaction = 0; //Reaction has a catalyst in it?
    double current_dG = 0.;
    double dndt_old;           //Total momentum correction term, split on reaction products according to their mass
    std::vector<double> dndts; //Momentum correction term per reactant, needs to be known for all reactants for all reactions
    
public:
    
    c_reaction(bool is_reverse, bool mtype_reaction, int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double a, double b, double c);
    c_reaction(bool is_reverse, bool mtype_reaction, int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double reaction_rate);
    void c_reaction_real(int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double reaction_rate);
    void set_base_pointer(c_Sim *base_simulation);
    
    void update_reaction_rate(double T);
    double get_reaction_rate(double T);
    void set_reac_number(int num) {reaction_number = num;}
};

/**
 * Photoochemistry reactions
 *
 * Reactions of the form A ->  A' + B' + ... , implementing 
 * terms of the form dn_A/dt = - G * n_A, dn_A'/dt = + G * n_A, 
 * where G will be typically along the lines of flux/photon energy * opacity.
 */
class c_photochem_reaction 
{
public:
    c_Sim *base;
    
    int reaction_number = -1;
    std::vector<int> educts;   //A list of indices
    std::vector<int> products; //A list of indices
    std::vector<double> e_stoch;
    std::vector<double> p_stoch;
    
    int e_num = 0.;
    int p_num = 0.;
    int count_p;
    double dndt_old;   //Total momentum correction term. No reactant-term is needed as we assume incoming photons carry negligible momentum
    
    //Photochem specific
    int band;                    //Minimum band to start this reaction being relevant (user needs to make sure this is compatible with threshold_energy!!!)
    double branching_ratio;
    double threshold_energy;
    double products_total_mass;
    double energy_split_factor;  //Ensure most of the heating goes to light products when mass ratios between products are extreme and an even split between reactants of equal mass
    
    c_photochem_reaction(int num_species, int num_bands, int band, std::vector<int> e_indices, std::vector<int> p_indices, std::vector<double> e_stoch, std::vector<double> p_stoch, double branching, double threshold);
    void set_reac_number(int num) {reaction_number = num;}
    void set_base_pointer(c_Sim *base_simulation);
};


/**
 * CLASS SIMULATION
 * 
 * Container for all simulation data.
 * Class Methods execute all physics modules, read input and write output data.
 * 
 */
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
    int cdebug = 0;
    int debug_cell = 40;
    int debug_steps = 8965;
    int use_self_gravity;
    int use_linear_gravity;
    int use_rad_fluxes;
    int suppress_warnings;
    int do_hydrodynamics;
    int photochemistry_level;
    int use_collisional_heating;
    int use_drag_predictor_step;
    int use_convective_fluxes;
    double convect_boundary_strength;
    double start_hydro_time;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  Numerical
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    int num_species;
    std::vector<c_Species> species;
    
    string speciesfile;
    string speciesfile_solo;
    string parfile;
    string fluxfile;
    string reactionfile;
    
    int num_bands_in;
    int num_bands_out;
    int num_he_bands;
    int num_plancks = 512;
    int i_wien;
    int i_rayleighjeans;
    double lT_spacing;
    char temperature_model;
    double T_crit_min;
    double T_crit_max; //Critical min/max temperatures beyond which one should extend the wavelength grid
    double fluxmultiplier; //In case one forgot a unit conversion...
    
    int num_cells;
    std::vector<double>dx;
    double domain_min, domain_max;
    double cells_per_decade;
    double grid2_cells_per_decade;
    double grid2_transition;
    int grid2_transition_i;
    int reverse_hydrostat_constrution;
    
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
    double lambda_min_in, lambda_max_in, lambda_per_decade_in;
    double lambda_min_out, lambda_max_out, lambda_per_decade_out;
    std::vector<double> l_i_in;
    std::vector<double> l_i_out;
    std::vector<double> l_i12_in;
    std::vector<double> l_i12_out;
    
    int init_geometry;
    char collision_model;
    char opacity_model;
    double ion_precision;
    double ion_heating_precision;
    double ion_floor;
    int ion_maxiter;
    int ion_heating_maxiter;
    double temperature_floor;
    double max_temperature;
    
    Geometry geometry ;
    
    double dt;
    double cflfactor;
    double t_max;
    double max_timestep_change;
    double dt_min_init;
    double dt_max;               //When wanting to cap the timestep in static simulations
    double timestep_rad2;
    double cfl_step ;
    double energy_epsilon;
    double globalTime;
    double output_time;
    double output_time_offset;
    double monitor_time;
    double max_snd_crs_time;
    double rad_energy_multiplier;
    signed long long int steps;
    int timecount;
    int monitor_output_index;

    IntegrationType order ;
    int num_ghosts ;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Boundaries
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    double right_extrap_press_multiplier;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Initial conditions
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int use_init_discont_smoothing;
    
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
    int    use_tides;
    double rhill;
    double init_star_mass;
    double ramp_star_mass_t0;
    double ramp_star_mass_t1;

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
    std::vector<double> enclosed_mass_tmp;
    //std::vector<double> total_pressure;
    std::vector<double> total_press_l ; // Reconstructed left/ right edges
    std::vector<double> total_press_r ;
    std::vector<double> total_adiabatic_index;
    int use_total_pressure;
    
    //
    // Friction
    //
    
    int friction_solver;
    int update_coll_frequently;
    double init_coll_factor;
    double coll_rampup_time;
    int use_avg_temperature;
    double avg_temperature_t0;
    double avg_temperature_t1;

    //using Matrix_t = Eigen::Matrix<double, NUM_SPECIES,NUM_SPECIES, Eigen::RowMajor>;
    //using Vector_t = Eigen::Matrix<double, NUM_SPECIES, 1>;
    Matrix_t friction_matrix_T;
    Matrix_t identity_matrix;
    Matrix_t friction_coefficients;
    Matrix_t friction_coeff_mask;
    Matrix_t inv_totmasses;
    Matrix_t resonant_pair_matrix;
    Vector_t friction_vec_input;
    Vector_t friction_vec_output;
    Vector_t friction_dEkin;
    Vector_t unity_vector;
    Eigen::PartialPivLU<Matrix_t> LU;
    Vector_t dens_vector;
    Vector_t numdens_vector;
    Vector_t mass_vector;
    Vector_t temperature_vector;
    Vector_t temperature_vector_augment;

    Matrix_t radiation_matrix_T;
    Matrix_t radiation_matrix_M;
    Vector_t radiation_vec_input;
    Vector_t radiation_vec_output;
    Vector_t radiation_cv_vector;
    Vector_t radiation_T3_vector;
    
    Eigen::VectorXd alphas_sample;
    Eigen::VectorXd friction_sample;
    
    double K_zz_init;
    std::vector<double> K_zz;
    
    //
    // Radiation
    //
    double T_star;
    double UV_star;
    double X_star;
    double Lyalpha_star;
    double R_star;
    
    double T_other;
    double R_other;
    double d_other;
    
    double T_int;  //Internal heat flux
    double T_planet; //Radiation flux making it through the atmosphere
    double R_core;
    double core_cv; //Heat capacity of the surface layer
    double bond_albedo;
    
    double init_radiation_factor;
    double radiation_rampup_time;
    
    //Indices for highenergy cooling
    int hnull_idx; 
    int hplus_idx;
    int e_idx;
    int C_idx;
    int Cp_idx;
    int Cpp_idx;
    int O_idx;
    int Op_idx;
    int Opp_idx;
    int h3plus_idx;
    
    //int radiation_solver;
    int use_planetary_temperature;
    int closed_radiative_boundaries ;
    int radiation_matter_equilibrium_test; //If set to 1, sets J = J_init in update_radiation()
    int radiation_diffusion_test_linear;
    int radiation_diffusion_test_nonlinear;
    int couple_J_into_T;
    double no_rad_trans;      // Multiplier for the div F radiation transport in the radiation solver to compare to models which don't cool thermally
    double CFL_break_time; //Numerical time after which cflfactor=0.9. Used in get_cfl_timestep()
    double photocooling_multiplier;
    double photocooling_expansion;
    
    std::vector<double> previous_monitor_J;
    std::vector<double> previous_monitor_T;

    Eigen::MatrixXd S_band;
    Eigen::MatrixXd dS_band;
    Eigen::MatrixXd dS_band_special; //Only used to document high-energy flux for C2ray solver case
    Eigen::MatrixXd dS_band_zero;
    Eigen::MatrixXd solar_heating;
    Eigen::MatrixXd solar_heating_final;
    std::vector<int> BAND_IS_HIGHENERGY;
    std::vector<double> photon_energies;
    const int HIGHENERGY_BAND_TRUE = 1;
    const int HIGHENERGY_BAND_FALSE= 0;
    
    Eigen::MatrixXd total_opacity;           //num_bands_out
    Eigen::MatrixXd cell_optical_depth;
    Eigen::MatrixXd radial_optical_depth;
    
    Eigen::MatrixXd total_opacity_twotemp;   //num_bands_in
    Eigen::MatrixXd cell_optical_depth_twotemp;
    Eigen::MatrixXd radial_optical_depth_twotemp;
    
    Eigen::MatrixXd cell_optical_depth_highenergy;     // Book-keeping for already used photons
    Eigen::MatrixXd highenergy_switch;                 // Switch off thermal heating for species in high-energy bands
    
    double *data_opacity[4], *opa_gas_tscale, *opa_gas_pscale, *opa_gas_ross, *opa_gas_planck;
    int opacity_gas_rows = 126, opacity_gas_cols = 94;
    double const_opacity_solar_factor;
    double const_opacity_rosseland_factor;
    double const_opacity_planck_factor;
    double init_J_factor;
    double init_T_temp;
    
    Eigen::MatrixXd planck_matrix;
    Eigen::MatrixXd Etot_corrected;
    Eigen::MatrixXd Jrad_FLD;

    BlockTriDiagSolver<Eigen::Dynamic> tridiag;
    Eigen::MatrixXd Jrad_init;
    
    int rad_solver_max_iter = 1;
    double epsilon_rad_min = 1e-1;  // Some convergence measure for the radiation solver, if needed
    double density_floor;
    double xi_rad;
    int solve_for_j;
    //
    // (Photo)Chemistry
    //
    int use_chemistry;
    int num_reactions;
    int num_photoreactions;
    int read_reactions_from_species_file;
    std::vector<c_reaction> reactions;
    std::vector<c_photochem_reaction> photoreactions;
    
    Vector_t reaction_b;
    Vector_t momentum_b;
    //Vector_t n_olds;
    //Vector_t n_news;
    std::vector<double> n_init;
    std::vector<double> n_tmp;
    Matrix_t reaction_matrix;
    Eigen::MatrixXd *reaction_matrix_ptr;
    Eigen::VectorXd *reaction_b_ptr;
    //Vector_t n_news;
    Eigen::PartialPivLU<Matrix_t> LUchem;
    
    Eigen::PartialPivLU<Matrix_t> *LUchem_ptr;
    
    Matrix_t chem_momentum_matrix;
    Eigen::PartialPivLU<Matrix_t> LUchem_mom;
    
    double chemistry_precision;
    int chemistry_maxiter;
    int chemistry_miniter;
    double chemistry_numberdens_floor;
    
    int intermediate_chemfloor_check;
    int chem_momentum_correction;
    int chem_ekin_correction;
    void init_reactions(int cdebug);
    void interpret_chem_reaction_list(string dir, string filename);
    void find_resonant_pairs(string dir, string filename);
    void do_chemistry(double timestep);
    
       int dt_skip_ichem;
    double dt_skip_dchem;
    
    Vector_t solver_cchem_implicit_general(double dt, int num_spec, int cdebug, const Vector_t& n_normalized, double ntot);
    int solver_cchem_implicit_specialized_cochem(double dt, int num_spec, int cdebug);
    
    std::vector<double> thermo_poly(double T,double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9);
    std::vector<double> get_thermo_variables(double T,string species_string);
    
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Variables for different scenarios
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int init_wind;
    double max_mdot;
    double init_sonic_radius;
    double T_increment;
    double alpha_collision;
    double alpha_collision_ions;
    double dust_to_gas_ratio;
    
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
    
    void print_monitor(int i);
    void print_diagnostic_file(int i);
    void write_into_execution_log(string dir, string par, string spcfile);
    
    void compute_total_pressure();
    int get_species_index(const string name, const int verbose);
    int find_closest_band(double energy_threshold);
    
    //
    // Friction
    //
    void compute_friction_analytical(); 
    void compute_friction_numerical(); 
    void compute_drag_update();
    
    void fill_alpha_basis_arrays(int j);
    void fill_rad_basis_arrays(int, double, Eigen::MatrixXd &, Eigen::MatrixXd &);
    void compute_alpha_matrix(int j);
    void compute_collisional_heat_exchange_matrix(int j) ;
    void compute_collisional_heat_exchange(); 
    
    //Gravity
    void init_grav_pot();
    void update_mass_and_pot();
    double get_phi_grav(double &r, double &mass);
    
    //Opacities 
    void init_malygin_opacities();
    double opacity_semenov_malygin(int rosseland, double temperature, double rho, double pressure, int caller_is_dust) ;
    double get_gas_malygin(int rosseland, double T_gas, double pressure);
    double freedman_opacity(double P, double T, double _Z);
    void kappa_landscape();
    
    //Radiation transport and chemistry
    
    void transport_radiation();     //  Called if use_rad_fluxes == 1
    void reset_dS();
    void update_dS();
    void update_dS_jb(int j, int b);
    void update_dS_jb_photochem(int j);
    void do_highenergy_cooling(int j);
    void update_tau_s_jb(int j, int b);
    void update_opacities();
    
    void do_photochemistry();
    void init_highenergy_cooling_indices();
    
    void update_fluxes(double timestep);           //  Called from transport_radiation#   
    void update_fluxes_FLD();           //  Called from transport_radiation#
    void update_fluxes_FLD_simple(double ddt);           //  Called from transport_radiation#
    void update_temperatures_simple();           //  Called from transport_radiation#
    void update_fluxes_FLD2(double, Eigen::MatrixXd &,Eigen::MatrixXd &,Eigen::MatrixXd &);           //  Called from transport_radiation#
    void update_temperatures(double, Eigen::MatrixXd &,Eigen::MatrixXd &,Eigen::MatrixXd &);
    double compute_planck_function_integral3(double lmin, double lmax, double temperature);
    double compute_planck_function_integral4(double lmin, double lmax, double temperature);
    
    
public:
    
    c_Sim() {};
    c_Sim(string parameter_filename, string species_filename, string workingdir, int debug=0, int debug_cell=1, int debug_steps=99999999);
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
    int num_bands_in;
    int num_bands_out;
    int debug;
    
    string opacity_data_string;
    string opacity_corrk_string;
    int num_opacity_datas;

    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Physics variables
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    double mass_amu;
    double inv_mass;
    int static_charge;
    double initial_fraction;
    
    std::vector<AOS> u, u0;          // Conserved hyperbolic variables: density, mass flux, energy density
    std::vector<AOS> dudt[2];        // Time derivative of u at each stage ;      
    std::vector<AOS> source;         // Gravitational source term
    std::vector<AOS> source_pressure;// Geometric source term
    std::vector<AOS> flux;
    std::vector<double> lconvect;
    std::vector<double> u_analytic;
    std::vector<double> K_zzf;     //(1+K_zz*n/b*mu/m_s) / (1+Kzz*n/b), which can be computed globally, so that every species s in cell i just needs to replace phi_i -> phi_i (1 + k_zzf / m_s)
    
    Eigen::MatrixXd opacity_data;  //superhigh-res lamdba grid, resolution dependent on minimum wl-distance in opacity data
    Eigen::VectorXd opacity_avg_solar;   //num_bands_in, this is the non-weighted average of opacity_data over the bands.
    Eigen::VectorXd opacity_avg_planck;   //num_bands_in, this is the non-weighted average of opacity_data over the bands.
    Eigen::VectorXd opacity_avg_rosseland;   //num_bands_in, this is the non-weighted average of opacity_data over the bands.
    
    Eigen::MatrixXd opacity;        //num_cells * num_bands_out, stands in for the Rosseland mean
    Eigen::MatrixXd opacity_planck; //num_cells * num_bands_out
    Eigen::MatrixXd opacity_twotemp;//num_cells * num_bands_in
    
    Eigen::MatrixXd fraction_total_solar_opacity; //num_cells * num_bands_in, what fraction of the total opacity is represented by this species in this cell. gives heating of the species
    Eigen::MatrixXd dS;  //num_cells * num_bands_out, Heating, low-energy and high-energy
    Eigen::MatrixXd dG;  //num_cells * num_bands_out, Cooling, high-energy
    Eigen::MatrixXd dGdT; //num_cells * num_bands_out, Cooling, high-energy
    double dlogOpa;
    
    //////////////
    //Custom opacity table stuff
    //////////////
    int opa_pgrid_size;
    int opa_tgrid_size;
    Eigen::VectorXd opa_pgrid;
    Eigen::VectorXd opa_tgrid;
    Eigen::VectorXd opa_grid_rosseland;
    Eigen::VectorXd opa_grid_planck;
    Eigen::VectorXd opa_grid_solar;
    
    Eigen::VectorXd opa_tgrid_log;
    Eigen::VectorXd opa_pgrid_log;
    Eigen::VectorXd opa_grid_rosseland_log;
    Eigen::VectorXd opa_grid_planck_log;
    Eigen::VectorXd opa_grid_solar_log;
    
    double density_excess;
    double bondi_radius;
    
    double const_T_space;
    double const_rho_scale;
    double const_opacity;
    int    is_dust_like;
    double pressure_broadening_factor;
    double pressure_broadening_exponent;
    
    //
    // Boundary condensation mass reservoir variables
    //
    double mass_reservoir;
    double p_sat;
    double t_evap;
    double latent_heat;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // Primitives and other helper quantities
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    std::vector<double> timesteps;
    std::vector<double> timesteps_cs;
    std::vector<double> timesteps_rad;
    std::vector<double> finalstep;
    std::vector<double> timesteps_de;
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
    void initialize_hydrostatic_atmosphere(string);
    void initialize_hydrostatic_atmosphere_iter(string);
    void initialize_exponential_atmosphere();
    void initialize_sound_wave();
    void compute_analytic_solution();
    
    int  init_static_atmosphere;
    void init_analytic_wind_solution();
    
    
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
                dp_final = - base->dx[i] * base->omegaplus[i] * (phi_s[i-1] - phi_s[i]) / (base->dx[i-1] + base->dx[i]);
                //dp_final = - base->dx[i] * base->omegaplus[i] * (base->phi[i-1]*K_zzf[i-1] - base->phi[i]*K_zzf[i]) / (base->dx[i-1] + base->dx[i]);
        } else {
                dp_final = - base->dx[i] * base->omegaminus[i] * (phi_s[i+1] - phi_s[i]) / (base->dx[i+1] + base->dx[i]);
                //dp_final = - base->dx[i] * base->omegaminus[i] * (base->phi[i+1]*K_zzf[i+1] - base->phi[i]*K_zzf[i]) / (base->dx[i+1] + base->dx[i]);
        }
            
        return dp_final;
    }

    void reconstruct_edge_states() ;
    
    AOS hllc_flux(int);
    AOS dust_flux(int);
    AOS source_grav(AOS &u, int &j);
    std::vector<double> phi_s;

    void compute_pressure(std::vector<AOS>& u) {
        eos->compute_primitive(&(u[0]), &(prim[0]), num_cells+2) ;    
        eos->compute_auxillary(&(prim[0]), num_cells+2);
    }
    
    void update_kzz_and_gravpot(int argument);
    
    void update_opacities();
    double interpol_tabulated_opacity(const Eigen::VectorXd& array, int band, double T_gas, double pressure);
    
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
    void user_opacity() ;

    void user_initial_conditions();
    void user_species_loop_function() ;
    
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    // File operations etc.
    //
    ////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    int suppress_warnings;
    void print_AOS_component_tofile(int timestepnumber);
    
    int read_species_data(string filename, int species_index);
    void read_opacity_table(string filename);
    void set_debug(int);
    void set_suppress_warnings(int j) {suppress_warnings = j;}    
    
    //c_Species();
    c_Species(c_Sim *base_simulation, string filename, string species_filename, int species_index, int debug=0);
    ~c_Species();
    void execute(std::vector<AOS>& u, std::vector<AOS>& dudt);
    
};

#endif//_AIOLOS_MAIN_H

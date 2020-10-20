#include "aiolos.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Pressure computation and matters of (Equation of) State
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//
// Help for variable definitions:
//
// E_total = u[i].u3 = E + 0.5 \rho u^2 = E + 0.5 * u[i].u2 * u[i].u2 / u[i].u1
// then we have E = \rho eint 
//
// The adiabatic equation of state usually looks P = (gamma_adiabat - 1.) E = (gamma_adiabat - 1.) \rho eint
// and the ideal gas EOS can be interpreted as   eint = c_v T
// from which the more familiar                  P = (gamma_adiabat-1.) cv \rho T = \rho k_B T / mu  follows
// when observing that \gamma_adiabat = cp/cv and cp-cv = kb / mu

//
// Compute pressure: The most important function for all kinds of applications inside the code.
//                   Is also the location of the equation of state (EOS). This should be the only instance in the code where P and the conserved variables are connected.
//
void c_Species::compute_pressure(std::vector<AOS>& u) {
    
    //Pressure now defined also on ghost cells, so that the analytic fluxes can be computed there
    for(int i=0;i<=num_cells+1;i++) {
        speed[i]           = u[i].u2 / u[i].u1;
        internal_energy[i] = u[i].u3 / u[i].u1 - 0.5 * speed[i] * speed[i];
        
        pressure[i]        = eos_pressure(u[i].u1, internal_energy[i], i); //(gamma_adiabat[i]-1.)*(u[i].u3 - 0.5* u[i].u2 * speed[i] );
        
        //TODO: Eliminate all other references in the code where the pressure is not computed via the EOS!
        
        if(i > 0)
            pressure_l[i]      = get_p_hydrostatic_nonuniform(i, u[i],  -1);
        if(i < num_cells+1)
            pressure_r[i]      = get_p_hydrostatic_nonuniform(i, u[i],  +1);
        
        internal_energy[i] = u[i].u3 / u[i].u1 - 0.5 * speed[i] * speed[i];
        temperature[i]     = internal_energy[i] / cv[i];
        cs[i]              = std::sqrt(gamma_adiabat[i] * pressure[i] / u[i].u1);
    }
    
}


double c_Species::eos_pressure(double dens, double eint, int i){
    double p;
    
    switch(eos_pressure_type) {
        case EOS_pressure_type::adiabatic:
            
            p = (gamma_adiabat[i] - 1.) * dens * eint;
            break;
        case EOS_pressure_type::polytropic:
            
            p = 1. * pow(dens, (gamma_adiabat[i] - 1.)/gamma_adiabat[i]);
            break;
        
        case EOS_pressure_type::supernova:
        {   
            double k1 = 4.897e14;
            double rhonuc = 1e14; 
            
            if(dens < rhonuc)
                p = k1 * pow(dens, (gamma_adiabat[i] - 1.)/gamma_adiabat[i]);
            else
                p = k1 / 1.5 * pow(rhonuc, gamma_adiabat[i] - 2.5) / (gamma_adiabat[i] - 1.) * pow(dens, 2.5);
            
            break;
        }
        case EOS_pressure_type::tabulated:
            
            p = interpolated_values_on_eos_p_grid(dens, eint);
            break;
        case EOS_pressure_type::user:
            
            p = eos_p_user(dens, eint);
            break;
        default:
            p = -1; //If we get this value something went wrong/
    }
    
    return p;
}

double c_Species::eos_eint(double dens, double temperature, int i) {
    double eint;
    
    switch(eos_internal_energy_type) {
        case EOS_internal_energy_type::thermal:
            
            eint = cv[i] * temperature;
            break;
        
        case EOS_internal_energy_type::constant:
            
            eint = 1.;
            break;
            
        case EOS_internal_energy_type::supernova:
        {
            double k1 = 4.897e14;
            double rhonuc = 1e14; 
            
            if(dens < 1e14)
                eint = k1 / (gamma_adiabat[i] - 1.) * pow(dens, gamma_adiabat[i]);
            else
                eint= k1 / 1.5 * pow(rhonuc, gamma_adiabat[i] - 2.5) * pow(dens, 2.5) + (2.5 - gamma_adiabat[i])/1.5*k1/(gamma_adiabat[i] - 1.) * dens * pow(rhonuc, gamma_adiabat[i] - 1);
            
            break;
        }
        case EOS_internal_energy_type::tabulated:
            
            eint = interpolated_values_on_eos_eint_grid(dens, temperature);
            break;
        case EOS_internal_energy_type::user:
            
            eint = eos_e_user(dens, temperature);
            break;
        default:
            eint = -1; //If we get this value something went wrong/
    }
    
    return eint;
}


//
// Tabulation functions
//
//

double c_Species::interpolated_values_on_eos_p_grid(double dens, double eint) {
    
    double p;
    
    throw std::runtime_error(
        "hydro::interpolated_values_on_eos_p_grid not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
    
    cout<<dens<<eint;    //Just use the variable so that the compiler shuts up about unused variables.
        
    // Step 1: Find dens in tab_p_x;
    // Step 2: Find eth in  tab_p_y;
    // Step 3: Define nearest support points on grid and interpolate p on tab_p
    
    return p;
}

double c_Species::interpolated_values_on_eos_eint_grid(double dens, double temperature) {
    
    double eint;
    
    throw std::runtime_error(
        "hydro::interpolated_values_on_eos_eint_grid not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
    
    cout<<dens<<temperature;    
        
    // Step 1: Find dens in tab_eint_x;
    // Step 2: Find eth in  tab_eint_y;
    // Step 3: Define nearest support points on grid and interpolate eint on tab_eint
    
    return eint;
}

//
// Function to read in pressure as function of density and internal energy from a file
//
//

void c_Species::read_tabulated_eos_data_pressure(string filename) {
    
    cout<<"Reading file "<<filename<<" in order to find me some good EOS data..."<<endl;
    //TODO:
    //tab_p = ...read in stuff...
    
    throw std::runtime_error(
        "hydro::read_tabulated_eos_data_pressure not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
}

//
// Function to read in internal energy as function of denisty and temperature from a file
//
//

void c_Species::read_tabulated_eos_data_eint(string filename) {
    
    cout<<"Reading file "<<filename<<" in order to find me some good EOS data..."<<endl;
    //TODO
    //tab_e = ...read in stuff...
    
    throw std::runtime_error(
        "hydro::read_tabulated_eos_data_eint not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
}


 void c_Sim::set_debug(int debug) {
     this->debug = debug;
}

 void c_Species::set_debug(int debug) {
     this->debug = debug;
}

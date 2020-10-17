
#include <stdexcept>

#include "aiolos.h"

void c_Species::user_boundary_left() {
    throw std::runtime_error(
        "hydro::user_boundary_left() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}


void c_Species::user_boundary_right() {
    throw std::runtime_error(
        "hydro::user_boundary_right() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}

void c_Species::user_initial_conditions(){
    throw std::runtime_error(
        "hydro::user_initial_conditions() not implemented. You must provide this "
        "method if you want to use user defined initial conditions.") ;
}

double c_Species::eos_p_user(double density, double eint) {
        throw std::runtime_error(
        "hydro::eos_p_user() not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
        
        cout<<density<<eint;
}

double c_Species::eos_e_user(double density, double temperature) {
        throw std::runtime_error(
        "hydro::eos_e_user() not implemented. You must provide this "
        "method if you want to use a user-defined equation of state.") ;
        
        cout<<density<<temperature;
}

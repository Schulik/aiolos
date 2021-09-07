
#include <stdexcept>

#include "aiolos.h"

void c_Species::user_boundary_left(std::vector<AOS>&) {
    throw std::runtime_error(
        "hydro::user_boundary_left() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}


void c_Species::user_boundary_right(std::vector<AOS>&) {
    throw std::runtime_error(
        "hydro::user_boundary_right() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}

void c_Species::user_initial_conditions(){
    throw std::runtime_error(
        "hydro::user_initial_conditions() not implemented. You must provide this "
        "method if you want to use user defined initial conditions.") ;
}

void c_Species::user_opacity() {
    throw std::runtime_error(
        "hydro::user_opacity() not implemented. You must provide this "
        "method if you want to use user defined initial conditions.") ;
}

void c_Species::user_species_loop_function() {} ;
void c_Sim::user_heating_function() {} ;
void c_Sim::user_output_function(int) {} ;
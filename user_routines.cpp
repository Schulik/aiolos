
#include <stdexcept>

#include "aiolos.h"

void hydro_run::user_boundary_left(std::vector<AOS>&) {
    throw std::runtime_error(
        "hydro::user_boundary_left() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}


void hydro_run::user_boundary_right(std::vector<AOS>&) {
    throw std::runtime_error(
        "hydro::user_boundary_right() not implemented. You must provide this "
        "method if you want to use user defined boundaries.") ;
}

void hydro_run::user_initial_conditions(){
    throw std::runtime_error(
        "hydro::user_initial_conditions() not implemented. You must provide this "
        "method if you want to use user defined initial conditions.") ;
}
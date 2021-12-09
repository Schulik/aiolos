
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

void c_Sim::user_heating_function() {} ;
void c_Sim::user_loop_function() {

  double gmm  = 5/3. ;
  double beta = 25 ; 

  double A=0, B=0, k=2*pi;
  for (int j=2; j < num_cells; j++) {
    double d = (species[0].prim[j].density - 1)*dx[j];
    A += 2*d*std::sin(k*x_i[j]) ;
    B += 2*d*std::cos(k*x_i[j]) ;
  }
    
  double x = dt / beta ;
  for (int j=1; j <num_cells+1; j++) {
    species[0].prim[j].pres = (species[0].prim[j].pres + x/gmm) / (1 + x) ;

    double g = B*std::sin(k*x_i[j]) - A*std::cos(k*x_i[j]) ;
    species[0].prim[j].speed -= dt * g ;
  }
  
  species[0].eos->compute_auxillary(&species[0].prim[0], num_cells+2);
  species[0].eos->compute_conserved(&species[0].prim[0], &species[0].u[0], num_cells+2);
  
} ;
void c_Sim::user_output_function(int) {} ;

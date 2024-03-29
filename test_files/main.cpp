#include "aiolos.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Begin Main
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//
// Important: Debug levels are ordered according to amount of information printed and depth into the algorithms
//
//  debug = 0: Only barebone info is printed, like "program started", "program ended."
//  debug = 1: General info about initialization is provided, nothing that is output every step
//  debug = 2: General info about every step is provided, like step number and determined dt. Also individual steps inside one timestep
//  debug = 3: General info is given about results of functions that are visited every timestep, like hllc_flux
//  debug = 4: Detailed info is given about every result of functions taht are visited every timestep. I.e. in hllc_flux all base data would now be provided how result was computed
//
//  

bool run_test(std::string test_name, std::string spc = "hydrogen.spc",
	      int suppress_warn=1, int debug=0) {
      
    try {
        c_Sim test(test_name, spc, "test_files/", debug) ;
        test.set_suppress_warnings(suppress_warn);

        test.execute(); 
    } catch (...) {
        std::cerr << "Test " << test_name << " crashed" << std::endl ;
        return false ;
    }
    std::cerr << "Test " << test_name << " completed" << std::endl ;

    return true ;
}

int main()
{
    int count = 0 ;
    count += !run_test("shock_tube1.par") ;
    count += !run_test("shock_tube2.par") ;
    count += !run_test("shock_tube3.par") ;
    count += !run_test("shock_tube4.par") ;
    count += !run_test("shock_tube5.par") ;
    count += !run_test("shock_tube6.par") ;
    count += !run_test("shock_tube7.par") ;

    count += !run_test("soundwave_32.par") ;
    count += !run_test("soundwave_64.par") ;
    count += !run_test("soundwave_128.par") ;
    count += !run_test("soundwave_256.par") ;
    count += !run_test("soundwave_512.par") ;
    
    count += !run_test("friction_2spc.par", "friction_2spc.spc") ;
    count += !run_test("friction_2spc_an.par", "friction_2spc.spc") ;
    count += !run_test("friction_2spc_phys.par", "friction_2spc.spc") ;
    count += !run_test("friction_2spc_phys_an.par", "friction_2spc.spc") ;
    count += !run_test("friction_6spc.par", "friction_6spc.spc") ;
    
    count += !run_test("collheat_2spc.par", "collheat_2spc.spc") ;
    count += !run_test("collheat_2spc_rad.par", "collheat_2spc.spc") ;

    count += !run_test("dustywave_nonstiff.par") ;
    count += !run_test("dustywave_stiff.par") ;

    count += !run_test("dusty_shock.par") ;
    
    count += !run_test("planet_cartesian.par") ;
    count += !run_test("planet_spherical.par") ;

    count += !run_test("irradiation.par", "irradiation.spc") ;

    
    return count ;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// End Main
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


************************************
************************************
*** AIOLOS - A 1-D radiation hydrodynamics multi-species code
************************************
************************************

In order to install and run Aiolos, we recommend the following steps:

1.) Clone this git repository (This presumably has just been done) via 'git clone https://github.com/Schulik/aiolos/'
    The default version of the code, with which the tests in our documentation paper were performed (Schulik&Booth2023) is 'v0.2', 
    as opposed to the main branch, on which active code development might be ongoing.
    We recommend switching to this branch.
2.) Install the EIGEN library via https://eigen.tuxfamily.org
3.) Compile, details below
4.) Run, details below


************************************
*** Compilation
************************************

aiolos can be compiled with make and gcc, just type 'make'. If you wish to provide
your own initial conditions or boundary conditions, add them to your own problem
file, e.g. "problems/my_problem.cpp". These can then be compiled into aiolos
via:

    make PROBLEM=my_problem


* Fixed number of species

By default, aiolos is compiled in a flexible mode that can run problesm with any
number of input species. However, a significant (10-20%) speed up can be
realised if the number of species is set at compile time. This can be
done via:

     make PROBLEM=my_problem NUM_SPECIES=4

where the "4" refers to a target of 4 species total.
A few warnings might occur concerning unused temporary variables. Those can be ignored.

************************************
*** In order to execute Aiolos, type:
************************************

   ./aiolos -par planet_spherical.par -spc mix3.spc

with the *.par and the *.spc file present into the command line.
This is a simple hydrostatic test problem, i.e. the initial density profile that is constructed should be kept perfectly.

A simple wind solution:

   ./aiolos -dir test_files/ -par planet_wind.par -spc mix_wind.spc
   
   
************************************
*** Brief description
************************************



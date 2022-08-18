In order to execute Aiolos in multispecies mode, type:

   ./aiolos -par planet_spherical.par -spc mix1.spc

with the *.par and the *.spc file present into the line command.


Compilation

aiolos can be compiled with make and gcc, just type make. If you wish to provide
your own initial conditions or boundary conditions, add them to your own problem
file, e.g. "problems/my_problem.cpp". These can then be compiled into aiolos
via:

    make PROBLEM=my_problem


Fixed number of species

By default, aiolos is compiled in a flexible mode that can run problesm with any
number of input species. However, a significant (10-20%) speed up can be
realised if the number of species is set at compile time. This can be
done via:

     make PROBLEM=my_problem NUM_SPECIES=4

where the "4" refers to a target of 4 species total.


Parallelization

aiolos uses a parallelised chemistry module to speed up complex networks. When compiled with the -fopenmp flag, -n <int> will run aiolos from the main command line 
with <int> OMP threads solving the chemistry loop in parallel.

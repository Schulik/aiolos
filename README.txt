In order to execute Aiolos in multispecies mode, type:

   ./aiolos -par planet_spherical.par -spc mix1.spc

with the *.par and the *.spc file present into the command line.


Compilation

aiolos can be compiled with make and gcc, just type make. If you wish to provide
your own initial conditions or boundary conditions, add them to your own problem
file, e.g. "problems/my_problem.cpp". These can then be compiled into aiolos
via:

    make PROBLEM=my_problem

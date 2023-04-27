************************************
************************************
*** AIOLOS - A 1-D radiation hydrodynamics multi-species code
************************************
************************************

************************************
*** 0. Download and branches
************************************

In order to install and run Aiolos, we recommend the following steps:

1.) Clone this git repository (This presumably has just been done) via 'git clone https://github.com/Schulik/aiolos/'
    The default version of the code, with which the tests in our documentation paper were performed (Schulik&Booth2023, hereafter SB23) is 'v0.2', 
    as opposed to the main branch, on which active code development might be ongoing.
    We recommend switching to this branch.
2.) Install the EIGEN library via https://eigen.tuxfamily.org
3.) Compile, details below
4.) Run, details below


************************************
*** 1. Compilation
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
*** 2. In order to execute Aiolos, type:
************************************

   ./aiolos -par planet_spherical.par -spc mix3.spc

with the *.par (parameter) and the *.spc (species) file present into the command line.
This is a simple hydrostatic test problem, i.e. the initial density profile that is constructed should be kept perfectly.
They command line does not require the species file to be added to it, it can be added to the parameter file instead.
Add into the parameter file
PARAMETER_FILE mix3.spc
then aiolos can be executed via 

   ./aiolos -par planet_spherical.par

A simple wind solution:

   ./aiolos -dir test_files/ -par planet_wind.par -spc mix_wind.spc
   
   
************************************
## 3. Execution parameters
************************************

As outlined in SB23, the code solves various differential equations, corresponding to code modules.
Code modules are hydrodynamics, friction, (thermal) radiation transport and (photo)chemistry.
Add the following keywords into the parameter file (**no tabs, only whitespaces**) to activate the modules.
Some parameters are prefixed with "PARI_", this is a leftover from code development and does not signify anything special now.
Parameters can be added in any line.

_________________________________
# 3.0 Basic execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation starts at t=0 and runs until tmax. 

_________________________________
3.1 Hydrodynamics module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DO_HYDRO 1
//options: 0: off, 1: on
//Off will allow to ignore the hydrodynamical CFL condition.

USE_TIDES 1
//options: 0: off, 1: use tidal field based on PARI_MSTAR and PARI_PLANET_DIST
_________________________________
3.2 Friction module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FRICTION_SOLVER 2
//options: 0: off, 1: analytic implicit, 2: numeric implicit

PARI_COLL_MODEL P
//options: C: constant, P: physical

PARI_ALPHA_COLL         1e0
//If coll_modell==C, then this is the constant value

_________________________________
3.3 Radiation transport module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PARI_USE_RADIATION  1
//options: 0: Off, 1: Full solver for coupled J_b and T_s variables. No limitation on the number of in or outgoing bands. 
//         2: 'Fast' solver, limited to one outgoing radiation band, any number of instellation bands. Radiation and temperatures are solved in a decoupled manner.

NO_RAD_TRANS 1.
//options: multiplier value for thermal radiative losses. If 0, only nonthermal cooling is activated. If no nonthermal cooling functions are specified, then absorption will be balanced by adiabatic cooling, i.e. energy-limited escape.

PARI_OPACITY_MODEL      P
//options: 
/* U: User defined opacity. Go into user_opacity() and specify your custom algorithm.
 * C: Constant opacities. Modifiers for solar, planck and rosseland opas exist individually.
 * P: 'Physical': Constant with simple pressure-broadening parameterization powerlaw above 0.1 bars
 * F: Freedman model. Only the fit for Rosseland opacities is currently included.
 * M: Malygin2014 (Planck Rosseland Gas opa)/Semenov2003 (Planck Rosseland Dust opa) combined opacities. Solar opas taken from *opa files.
 * D: Dust only model, given a certain dust size.
 * T: Tabulated opacities. Takes p-T dependent solar, planck and rosseland data from *aiopa files.
 * K: 'Kombined' opacities. Planck and Rosseland is p-T dependent from *aiopa data and solar are p-T-constant, but spectrally resolved from *opa files.
 */

_________________________________
3.4 (Photo) chemistry module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PHOTOCHEM_LEVEL 1
//options: 0: off, 1: C2Ray scheme as described in SB23, 2: General photo and thermochemistry solver

_________________________________
3.5 Species files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SPECIES_FILE  mix1.spc
//options: Species filename

In *spc files, column descriptions:
# 1:Number 2:Name 3:mass in amu 4:dof              5:electrostatic charge    6:relative amount 7:initial density excess 8: is_dust_like 9:opacity.opa
@ 0        S0     1.0         3.                   0                         1.0                 -0.0                      0              highenergy_atomichydrogen.opa

//options: number:       running number from 0 to s-1
// name:                 string, used to designate output files for species
// mass in amu:          self-explanatory
// dof:                  degrees of freedom for this species, defining adiabatic index as gamma_ad = (dof+2)/dof
// electrostatic charge: in elementary charges, can be positive or negative
// relative amount:      relative density in boundary cell
// initial density excess: density jump at discontinuity, in case PARI_INIT_WIND 1  e.g. if value is -0.9 density jump of ratio 10, for -0.99 we get a jump of ratio 100 etc.
// is_dust_like:          0 or 1
// opacity.opa            opacity files, depending on the preset opacity model, *opa, *op2, *op3 or *aiopa files are expected here.

_________________________________
3.6 Opacity files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See example files in inpufiles/







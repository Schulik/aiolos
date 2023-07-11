#include "aiolos.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Begin Main
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * Main routine.
 * 
 * Intent for debug levels (might not be consistent in the entire program)
 *   debug = 0(default): Only barebone info is printed, like "program started", "program ended."
 *   debug = 1: General info about initialization is provided, nothing that is output every step
 *   debug = 2: General info about every step is provided, like step number and determined dt. Also individual steps inside one timestep
 *   debug = 3: General info is given about results of functions that are visited every timestep, like hllc_flux
 *   debug = 4: Detailed info is given about every result of functions taht are visited every timestep. I.e. in hllc_flux all base data would now be provided how result was computed
 * 
 * @param[in] argc Standard C++ argument counter
 * @param[in] argv Standard C++ command line string list. Used to communicate key parameters like the *par and *spc files, as well as the debug level to the program.
 */
int main(int argc, char** argv)
{
    string simulationname;
    string speciesfile;
    string workingdir = "./";
    string tempintent = "---";
    int debug;
    int debug_cell = 40;
    int debug_steps = 8965;
    int suppress_warnings_global = 0;
    int external_thread_num = 1;
    cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
          cout<<"~~~ Welcome to AIOLOS! May a gentle breeze lead your way through the bugs."<<endl;
          cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
    
    //
    // Loop through all args, in order to find the name of the parameter file (KEY), debug options etc.
    //
    cout<<std::scientific;
    cout<<"Use -par name.par to let the program know the name of the required parameter file."<<endl;
    cout<<"Use -spc name.spc to let the program know the name of the required species file."<<endl;
    cout<<"Optional command line parametrs: -debug <int>, for debug level; -war <int>, for suppressing warnings."<<endl;
    
    int parameterfile_found = 0;
    int speciesfile_found   = 0;
    debug = 0;
    
    for(int i=0; i<argc; i++) {
        
        string tmpstring = argv[i];
        
        if(tmpstring.compare("-par") == 0) {
            simulationname      = argv[i+1];
            parameterfile_found = 1;
            cout<<"Attempting at assingning a simulationname to simulationname: "<<simulationname<<endl;
            i++;
        }
        if(tmpstring.compare("-spc") == 0) {
            speciesfile         = argv[i+1];
            speciesfile_found   = 1;
            cout<<"Attempting at assingning the following string to speciesfile: "<<speciesfile<<endl;
            i++;
        }
        if(tmpstring.compare("-dir") == 0) {
            workingdir         = argv[i+1];
            cout<<"Attempting at assingning the following string to workingdir: "<<workingdir<<endl;
            i++;
        }
        if(tmpstring.compare("-debug") == 0) {
            string yet_another_string = argv[i+1];
            debug      = std::stoi(yet_another_string);
            cout<<"Accepted debug modus "<<debug<<endl;
            i++;
        }
        if(tmpstring.compare("-war") == 0) {
            string yet_another_string = argv[i+1];
            suppress_warnings_global =   std::stoi(yet_another_string);
            
            cout<<"Suppressing warnings yes/no = 1/0 : "<<suppress_warnings_global<<endl;
            i++;
        }
        if(tmpstring.compare("-dcell") == 0) {
            string yet_another_string = argv[i+1];
            debug_cell    = std::stoi(yet_another_string);
            cout<<"Debugging cell "<<debug<<endl;
            i++;
        }
        if(tmpstring.compare("-dsteps") == 0) {
            string yet_another_string = argv[i+1];
            debug_steps    = std::stoi(yet_another_string);
            cout<<"Debugging info output after steps number "<<debug<<endl;
            i++;
        }
        if(tmpstring.compare("-n") == 0) {
            string yet_another_string = argv[i+1];
            external_thread_num      = std::stoi(yet_another_string);
            i++;
        }
        if(tmpstring.compare("-int") == 0) {
            tempintent      = argv[i+1];
            cout<<"Intent comment found as "<<tempintent<<endl;
            i++;
        }
    }
    
    if(!parameterfile_found) {
            cout<<"No simulationname found, chosing default simulationname: simulation.par"<<endl;
            simulationname = "simulation.par";
    }
    if(!speciesfile_found) {
            cout<<"No speciesfile found on command line, chosing default speciesfile until parameterfile tells us otherwise: default.spc"<<endl;
            speciesfile = "default.spc";
    }
        
    #if defined(_OPENMP)
        omp_set_num_threads(external_thread_num);    
        cout<<"Running with OMP, omp_thread_num = "<<omp_get_max_threads()<<endl;
    #endif
    
    try {
        
        cout<<endl<<"In main, construction of simulation is about to start."<<endl;
       
        //Main simulation class object, is initialized with the simulation parameters from a file
        c_Sim simulation1(simulationname, speciesfile, workingdir, tempintent, debug, debug_cell, debug_steps);
        
        simulation1.set_suppress_warnings(suppress_warnings_global);

        cout<<"In main, execution is about to start."<<endl;
        
        simulation1.execute();
    }
    catch (int err){
        
        switch(err) {
            case 0: cout<<"np_zeros allocation failed!"<<endl; break ;
            case 1: cout<<"np_ones allocation failed!"<<endl; break ;
            case 2: cout<<"np_zeros allocation failed!"<<endl; break ;
            default:cout<<"Unknown integer error occured!"<<endl; break ;
        } 
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc caught in main simulation block: " << ba.what() <<endl;
    }
    catch(std::exception &e) {
        cout << e.what() << endl;
    }
    catch(...) {
        cout<<"Unknown error in initialization of main variables!"<<endl;
    }
    
    //Running class desctructors manually apparently not needed in modern C++ anymore
    
    return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// End Main
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


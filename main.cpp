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

int main(int argc, char** argv)
{
    string simulationname;
    int debug;
    int suppress_warnings_global = 0;
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    cout<<"Welcome to AIOLOS! May a gentle breeze lead your way home. argc="<<argc<<endl;
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
    
    //
    // Loop through all args, in order to find the name of the parameter file (KEY), debug options etc.
    //
    cout<<std::scientific;
    cout<<"Use -par name.par to give key parameters to the executable."<<endl;
    
    int parameterfile_found = 0;
    debug = 0;
    
    for(int i=0; i<argc; i++) {
        
        string tmpstring = argv[i];
        
        if(tmpstring.compare("-par") == 0) {
            simulationname      = argv[i+1];
            parameterfile_found = 1;
            cout<<"Attempting at assingning a simulationname to simulationname: "<<simulationname<<endl;
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
        
    }
    
    if(!parameterfile_found) {
            cout<<"No simulationname found, chosing default simulationname: simulation.par"<<endl;
            simulationname = "simulation.par";
    }
        
    
    try {

        cout<<"In main, construction of simulation is about to start."<<endl;
       
        //Main simulation class object, is initialized with the simulation parameters from a file
        hydro_run simulation1(simulationname, debug);
        
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
        std::cerr << "bad_alloc caught in hydro init block: " << ba.what() <<endl;
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


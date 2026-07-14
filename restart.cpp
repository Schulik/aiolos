/*
 * Restart.cpp 
 * 
 * Contains routines to restart aiolos from a previous output
 */

#define EIGEN_RUNTIME_NO_MALLOC

#include <array>
#include <cassert>

#include "aiolos.h"

/*
 * Restartmain
 * 
 * Searches for restart data, reads it in, and overwrites initial conditions with that data, so that the simulation can continue from the savenumber.
 * 
 */
void c_Sim::restart_from_outputnumber(int restartnumber, double restarttime_cmdline, int restartmode) {
    
    
    cout<<" In restart, restartnumber = "<<restartnumber<<" enter a char to continue. retime/mode = "<<restarttime_cmdline<<" / "<<restartmode<<endl;
    //char a;
    //cin>>a;
    
    string fappendix = "";
    if(restartmode > 0.5)
        fappendix = "RE"; //Restart from already restarted files
    
    //
    // Locate output files to read in with right number
    //
    int allfilesexist = 1;
    for(int s=0; s<num_species; s++) {
        string filename ;
        {
            stringstream filenamedummy;
            string truncated_name = stringsplit(simname,".")[0];
            filenamedummy<<workingdir<<"output_"<<truncated_name<<"_"<<species[s].speciesname<<"_t"<<restartnumber<<fappendix<<".dat";
            filename = filenamedummy.str() ;
        }
        ifstream infile(filename, ios::in);
        if (!infile.is_open()) {
            cout<<filename<<" does not exist!"<<endl;
            allfilesexist = 0;
        }
            
    }
    
    //
    // If all files exist, continue and read data
    //
    if(allfilesexist) {
        
        cout<<" All files exist, assigning data.."<<endl;
        
        for(int s=0; s<num_species; s++) {
            string filename ;
            {
                stringstream filenamedummy;
                string truncated_name = stringsplit(simname,".")[0];
                filenamedummy<<workingdir<<"output_"<<truncated_name<<"_"<<species[s].speciesname<<"_t"<<restartnumber<<fappendix<<".dat";
                filename = filenamedummy.str() ;
            }
            ifstream infile(filename, ios::in);
            if (!infile.is_open())
                cout<<filename<<" not open although it exists!"<<endl;
        
            //Read data
            std::vector<AOS> tmp = init_AOS(num_cells+2);
            string line;
            int num_lines = 0;
            
            while(std::getline( infile, line )) //The allocated array should not have changed after the restart
            {
                std::vector<string> datalist = stringsplit(line,"\t");
                //tmp[num_lines].u1 = std::stod(datalist[1]);
                //tmp[num_lines].u2 = std::stod(datalist[2]);
                //tmp[num_lines].u3 = std::stod(datalist[3]);
                
                species[s].u[num_lines].u1 = std::stod(datalist[1]);
                species[s].u[num_lines].u2 = std::stod(datalist[2]);
                species[s].u[num_lines].u3 = std::stod(datalist[3]);
                
                num_lines++;
                
                if(num_lines==1) {
                    //cout<<"reading "<<tmp[num_lines].u1
                    cout<<"Restarting species = "<<species[s].speciesname<<" reading boundary dens = "<<species[s].u[num_lines].u1<<" from "<<filename<<endl;
                    //char cc;
                    //cin>>cc;
                }
            }
            
            //for(int i=1; i<= num_cells; i++) {}
            
            //No need to recompute primitives and auxiliaries, because directly after restart_from_outputnumber(), compute_pressure() is called, which does both
            
        }
        
    } else {
        cout<<" ERROR IN RESTART: Not all files exist for chosen input! "<<endl;
        char b;
        cin>>b;
    }
    
    restarttime = output_time * ((int)restartnumber);   
    if (restarttime_cmdline > 0) {
        restarttime = restarttime_cmdline;
        globalTime  = restarttime_cmdline;
        cout<<" found nonzero restarttime "<<restarttime_cmdline<<endl;
        //char a;
        //cin>>a;
    }
    
    
    
    //
    //
    //
    for(int s = 0; s < num_species; s++)
        species[s].compute_pressure(species[s].u);
    
}

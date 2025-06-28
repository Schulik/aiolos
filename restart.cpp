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
void c_Sim::restart_from_outputnumber(int restartnumber) {
    
    char a;
    cout<<" In restart, restartnumber = "<<restartnumber<<" enter a char to continue."<<endl;
    cin>>a;
    
    
    
    //
    // Locate output files to read in with right number
    //
    int allfilesexist = 1;
    for(int s=0; s<num_species; s++) {
        string filename ;
        {
            stringstream filenamedummy;
            string truncated_name = stringsplit(simname,".")[0];
            filenamedummy<<workingdir<<"output_"<<truncated_name<<"_"<<species[s].speciesname<<"_t"<<restartnumber<<".dat";
            filename = filenamedummy.str() ;
        }
        ifstream infile(filename, ios::in);
        if (!infile.is_open())
            allfilesexist = 0;
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
                filenamedummy<<workingdir<<"output_"<<truncated_name<<"_"<<species[s].speciesname<<"_t"<<restartnumber<<".dat";
                filename = filenamedummy.str() ;
            }
            ifstream infile(filename, ios::in);
            if (!infile.is_open())
                cout<<filename<<" not open although it exists!"<<endl;
        
            //Read data
            std::vector<AOS> tmp = init_AOS(num_cells+2);
            string line;
            int num_lines = 0;
            
            
            while(std::getline( infile, line )) 
            {
                std::vector<string> datalist = stringsplit(line,"\t");
                //tmp[num_lines].u1 = std::stod(datalist[1]);
                //tmp[num_lines].u2 = std::stod(datalist[2]);
                //tmp[num_lines].u3 = std::stod(datalist[3]);
                
                species[s].u[num_lines].u1 = std::stod(datalist[1]);
                species[s].u[num_lines].u2 = std::stod(datalist[2]);
                species[s].u[num_lines].u3 = std::stod(datalist[3]);
                
                num_lines++;
                
                if(num_lines==100)
                    //cout<<"reading "<<tmp[num_lines].u1
                    cout<<"Restarting species = "<<species[s].speciesname<<" "<<species[s].u[num_lines].u1<<endl;
            }
            
            //for(int i=1; i<= num_cells; i++) {}
            
            //No need to recompute primitives and auxiliaries, because directly after restart_from_outputnumber(), compute_pressure() is called, which does both
            
        }
        
    } else {
        cout<<" ERROR IN RESTART: Not all files exist for chosen input! "<<endl;
    }
    
    restarttime = output_time * ((int)restartnumber);   
    
    
    
    
    //
    //
    //
    for(int s = 0; s < num_species; s++)
                species[s].compute_pressure(species[s].u);
    
}

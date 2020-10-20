#include <sstream>
#include <stdexcept>
#include "aiolos.h"




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// CFL Timestep
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


double c_Sim::get_cfl_timestep() {

    //
    // Compute individual max wave crossing timesteps per cell
    //  t = delta x / v = delta x / momentum / density
    //
    double minstep = 0.; 
    //TODO: Compute sound crossing time for entire domain
    //double finalstep = 0.;
    
    max_snd_crs_time=0;
    for(int s=0; s < num_species; s++) {
        
        species[s].snd_crs_time = 0;
        for(int i=1; i<=num_cells; i++) {
            
            //Computing the inverse timesteps first
            species[s].timesteps[i]    = std::abs(species[s].prim[i].speed / dx[i]); 
            species[s].timesteps_cs[i] = species[s].prim[i].sound_speed / dx[i];
            species[s].finalstep[i]    = species[s].timesteps[i] + species[s].timesteps_cs[i] ;
            
            species[s].snd_crs_time += 2.* dx[i] / species[s].prim[i].sound_speed ;
            
            minstep = std::max(minstep, species[s].finalstep[i]) ;
        }
        max_snd_crs_time = std::max(max_snd_crs_time, species[s].snd_crs_time) ;
    }
    //Invert and apply CFL secutiry factor
    minstep = cflfactor / minstep;
    
    //Check if the simulation step is too large (e.g. when v=0 in most of the domain), then make the step artificially smaller
    if ( minstep > t_max)
        return t_max * 1e-2;
    
    return minstep;
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Helper functions
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


std::vector<AOS> init_AOS(int num) {   
    return std::vector<AOS>(num);
}

//

//
// Return a 1-D array of zeroes, identical to the numpy function
//
std::vector<double> np_zeros(int size) {

    return std::vector<double>(size, 0.0) ;
}

//
// Return a 1-D array of one, identical to the numpy function
//
std::vector<double> nnp_ones(int size) {

    return std::vector<double>(size, 1.0) ;
}
//
// Return a 1-D array of a constant value, identical to the numpy function
//
std::vector<double> np_somevalue(int size, double set_value) {

     return std::vector<double>(size, set_value) ;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// String split function, returns vector of split strings delimited by delim of initial string str
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vector<string> stringsplit(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Read parameters: Looks for the value of one parameter in a file
// 
// Takes:
//      string filename: The filename that should contain all parameters and their values
//      string variablename: the name of one individual parameter
//      int desired_vartype: what kind of value are we looking for, int, double or string?
//
// Returns:
//      simulation_parameter structure containing the name, type and value of a parameter found in a file
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void c_Species::read_species_data(string filename, int species_index) {
    
    
    ifstream file(filename);
    string line;
    if(!file) {
        cout<<"Couldnt open species file "<<filename<<"!!!!!!!!!!1111"<<endl;
        
    }
    //simulation_parameter tmp_parameter = {"NaN",0,0.,0,"NaN"};
    int found = 0;
    
    if(debug > 0) cout<<"          In read species Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        std::vector<string> stringlist = stringsplit(line," ");

        //if(line.find(variablename) != string::npos) {
        if(stringlist[0].find("@") != string::npos) {
            
            if(std::stoi(stringlist[1]) == species_index) {
                
                found = 1;
                
                this->name             = stringlist[2];
                this->mass_amu         = std::stod(stringlist[3]);
                this->cv               = std::stod(stringlist[4]);
                this->gamma_adiabat    = std::stod(stringlist[5]);
                this->initial_fraction = std::stod(stringlist[6]);
                
                if(debug > 0)
                    cout<<"Found species called "<<name<<" with a mass of "<<mass_amu<<" cv "<<cv<<" gamma "<<gamma_adiabat<<" and zeta_0 = "<<initial_fraction<<endl;
                
            }

        }
    
    }
    
    //cout<<"    In read_parameter Pos2"<<endl;
    
    if(found == 0) {
            
            if(debug > 0)
                cout<<"WARNING: Species number "<<species_index<<" not found in parameterfile!"<<endl;
    }
    if(found > 1) {
            cout<<"WARNING: Species number "<<species_index<<" defined more than once in parameterfile!"<<endl;
    }
    
    file.close();
    
    if(debug > 0) cout<<"         Leaving species readin now. Bye!"<<endl;
}


template<typename T> 
T parse_item(string item) {
    std::stringstream ss(item);
    T val;
    ss >> val;
    return val;
}

class ParameterError : public std::invalid_argument {
  public:
    ParameterError(int count_, string msg)
      : std::invalid_argument(msg), count(count_)
    { } ;
  
  int count ;
} ;

template<typename T>
simulation_parameter<T> read_parameter_from_file(string filename, string variablename, int debug) {
    
    ifstream file(filename);
    string line;
    if(!file) {
        std::stringstream error ;
        error << "Could not open parameter file:" << filename << endl ;
        throw std::runtime_error(error.str()) ;
    }

    simulation_parameter<T> tmp_parameter = {"NaN",T(),"NaN"};
    int found = 0;
    
    //cout<<"    In read_parameter Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        //cout<<"    In read_parameter Pos1.5: "<<line<<endl;
        
        if(line.find(variablename) != string::npos) {
            tmp_parameter.name = variablename;
            
            tmp_parameter.value = parse_item<T>(stringsplit(line," ")[1]) ;

            if(debug > 0)
                cout<<"Found variable called "<<variablename<<" and set it to "<<tmp_parameter.value<<endl;
            found++;
        }
    }
    
    //cout<<"    In read_parameter Pos2"<<endl;
    
    if(found == 0) {
        std:: stringstream error ;
        error << "ERROR: Variable "<<variablename<<" not found in parameterfile!"<<endl;
        throw ParameterError(0, error.str()) ;
    }
    if(found > 1) {
        std:: stringstream error ;
        error <<"ERROR: Variable "<<variablename<<" defined more than once in parameterfile!"<<endl;
        throw ParameterError(found, error.str()) ;
    }
    
    file.close();
    return tmp_parameter;
}


template<typename T>
simulation_parameter<T> read_parameter_from_file(string filename, string variablename, int debug, T default_) {
    simulation_parameter<T> parameter ;

    try {
       parameter = read_parameter_from_file<T>(filename, variablename, debug) ;
    } catch (ParameterError& e) {
        if (e.count == 0) {
            parameter.name = variablename ;
            parameter.value = default_ ;
            cout<<"WARNING: Parameter "<<variablename<<" not found in "<<filename<<". Assigning default value of "<<default_<<endl;
        }
        else
            throw ; // re-throw because we can't help with duplicated parameters
    }

    return parameter ;
}


template simulation_parameter<bool>  read_parameter_from_file(string, string, int, bool);
template simulation_parameter<int>   read_parameter_from_file(string, string, int, int);
template simulation_parameter<double> read_parameter_from_file(string, string, int, double);
template simulation_parameter<string> read_parameter_from_file(string, string, int, string);

template simulation_parameter<Geometry> read_parameter_from_file(string, string, int, Geometry);
template simulation_parameter<BoundaryType> read_parameter_from_file(string, string, int, BoundaryType);
template simulation_parameter<IntegrationType> read_parameter_from_file(string, string, int, IntegrationType);

//Print 2 
void c_Species::print_AOS_component_tofile(int timestepnumber) {
                                              
    // Choose a sensible default output name
    string filename ;
    {
        stringstream filenamedummy;
        string truncated_name = stringsplit(base->simname,".")[0];
        filenamedummy<<"output_"<<"_"<<truncated_name;
        filename = filenamedummy.str() ;
    }

    filename = read_parameter_from_file<string>(base->simname, "OUTPUT_FILENAME", debug, filename).value ;

    // Add the snap number
    {
        stringstream filenamedummy;
        filenamedummy << filename<<"_"<< name << "_t" << timestepnumber << ".dat";
        filename = filenamedummy.str() ;
    }

    if(debug > 1)
        cout<<"Trying to open file "<<filename<<endl;
    
    ofstream outfile(filename, ios::out);
    
    if (outfile.is_open())
    {
        //outfile.precision(16);
        
        double hydrostat = 0., hydrostat2 = 0., hydrostat3 = 0.;
        
        //Print left ghost stuff
        outfile<<base->x_i12[0]<<'\t'<<u[0].u1<<'\t'<<u[0].u2<<'\t'<<u[0].u3<<'\t'<<flux[0].u1<<'\t'<<flux[0].u2<<'\t'<<flux[0].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[0].pres<<'\t'<<u[0].u2/u[0].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[0]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
        
        //Print the domain
        for(int i=1; i<=num_cells; i++) {
            
            hydrostat = flux[i-1].u2/base->dx[i] ; //hydrostat2 + hydrostat3 ; 
            hydrostat2 = flux[i].u2/base->dx[i];//pressure[i+1] - pressure[i];
            hydrostat3 = source[i].u2;//0.5 * (u[i].u1 + u[i+1].u1) * (phi[i+1] - phi[i]);
            
            //outfile<<x[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<flux[i+1].u1<<'\t'<<flux[i+1].u2<<'\t'<<flux[i+1].u3<<'\t'<<pressure[i]<<'\t'<<u[i].u2/u[i].u1<<'\t'<<internal_energy[i] <<'\t'<<timesteps[i]<<'\t'<<phi[i]<<'\t'<<timesteps_cs[i]<<'\t'<<opticaldepth[i]<<'\t'<<radiative_flux[i] <<endl;
            
            //flux[i] - flux[i+1] + source[i]
            
            outfile<<base->x_i12[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<((flux[i-1].u1 - flux[i].u1)/base->dx[i] + source[i].u1)<<'\t'<<((flux[i-1].u2 - flux[i].u2)/base->dx[i] + source[i].u2)<<'\t'<<((flux[i-1].u3 - flux[i].u3)/base->dx[i] + source[i].u3)<<'\t'<<prim[i].pres<<'\t'<<u[i].u2/u[i].u1<<'\t'<<prim[i].internal_energy <<'\t'<<timesteps[i]<<'\t'<<base->phi[i]<<'\t'<<timesteps_cs[i]<<'\t'<<opticaldepth[i]<<'\t'<<hydrostat<<'\t'<<hydrostat2<<'\t'<<hydrostat3<<'\t'<<base->enclosed_mass[i]<<endl;
        }
        
        //Print right ghost stuff
        outfile<<base->x_i12[num_cells+1]<<'\t'<<u[num_cells+1].u1<<'\t'<<u[num_cells+1].u2<<'\t'<<u[num_cells+1].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[num_cells+1].pres<<'\t'<<u[num_cells+1].u2/u[num_cells+1].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[num_cells+1]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
   
        cout<<"    Sucessfully written file "<<filename<<" for species = "<<name<<" at time "<<base->globalTime<<" and dt="<<base->dt<<endl;
    }
    else cout << "Unable to open file" << filename << endl; 
    outfile.close();
    
 
}

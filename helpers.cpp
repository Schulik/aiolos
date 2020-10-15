#include <sstream>
#include <stdexcept>
#include "aiolos.h"

std::vector<AOS> init_AOS(int num) {   
    return std::vector<AOS>(num);
}

//
//
//
void hydro_run::compute_pressure(const std::vector<AOS>& u) {
    
    //Pressure now defined also on ghost cells, so that the analytic fluxes can be computed there
    for(int i=0;i<=num_cells+1;i++) {
        speed[i]           = u[i].u2 / u[i].u1;
        pressure[i]        = (gamma_adiabat-1.)*(u[i].u3 - 0.5* u[i].u2 * speed[i] );
        
        if(i > 0)
            pressure_l[i]      = get_p_hydrostatic_nonuniform(i,  -1);
        if(i < num_cells+1)
            pressure_r[i]      = get_p_hydrostatic_nonuniform(i,  +1);
        
        internal_energy[i] = u[i].u3 / u[i].u1 - 0.5 * speed[i] * speed[i];
        temperature[i]     = internal_energy[i] / cv;
        cs[i]              = std::sqrt(gamma_adiabat * pressure[i] / u[i].u1);
    }
    
}

 void hydro_run::set_debug(int debug) {
     this->debug = debug;
}

AOS hydro_run::analytic_flux(AOS &input_vec, const int &j) {
     
    return AOS (input_vec.u2, input_vec.u2*input_vec.u2/input_vec.u1 + pressure[j], input_vec.u2/input_vec.u1 * (input_vec.u3 + pressure[j]) );
}

double hydro_run::get_cfl_timestep() {

    //
    // Compute individual max wave crossing timesteps per cell
    //  t = delta x / v = delta x / momentum / density
    //
    double minstep = 0.; 
    //TODO: Compute sound crossing time for entire domain
    //double finalstep = 0.;
    
    snd_crs_time=0;
    for(int i=1; i<=num_cells; i++) {
        
        //Computing the inverse timesteps first
        timesteps[i]    = std::abs(u[i].u2 / u[i].u1 / dx[i]); 
        timesteps_cs[i] = std::abs(cs[i] / dx[i]);
        finalstep[i]    = std::sqrt(timesteps[i]*timesteps[i] + timesteps_cs[i] * timesteps_cs[i]);
        
        //cout<<" timestep["<<i<<"]="<<timesteps[i]<<endl;
        snd_crs_time += 2.* dx[i] / cs[i];
        
        if (finalstep[i] > minstep)
            minstep = finalstep[i];  
    }
    
    //Invert and apply CFL secutiry factor
    minstep = cflfactor * 0.4 / minstep;
    
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
        }
        else
            throw ; // re-throw because we can't help with duplicated parameters
    }

    return parameter ;
}


template simulation_parameter<bool>  read_parameter_from_file(string, string, int, bool);
template simulation_parameter<int> read_parameter_from_file(string, string, int, int);
template simulation_parameter<double> read_parameter_from_file(string, string, int, double);
template simulation_parameter<string> read_parameter_from_file(string, string, int, string);

template simulation_parameter<Geometry> read_parameter_from_file(string, string, int, Geometry);
template simulation_parameter<BoundaryType> read_parameter_from_file(string, string, int, BoundaryType);


//Print 2 
void hydro_run::print_AOS_component_tofile(const std::vector<double>& x, 
                                          const std::vector<AOS>& data,
                                          const std::vector<AOS>& fluxes,
                                          int timestepnumber) {
                                              
    // Choose a sensible default output name
    string filename ;
    {
        stringstream filenamedummy;
        string truncated_name = stringsplit(simname,".")[0];
        filenamedummy<<"output_"<<truncated_name;
        filename = filenamedummy.str() ;
    }

    filename = read_parameter_from_file<string>(simname, "OUTPUT_FILENAME", debug, filename).value ;

    // Add the snap number
    {
        stringstream filenamedummy;
        filenamedummy << filename << "_t" << timestepnumber << ".dat";
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
        outfile<<x[0]<<'\t'<<data[0].u1<<'\t'<<data[0].u2<<'\t'<<data[0].u3<<'\t'<<fluxes[0].u1<<'\t'<<fluxes[0].u2<<'\t'<<fluxes[0].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<pressure[0]<<'\t'<<data[0].u2/data[0].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<phi[0]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
        
        //Print the domain
        for(int i=1; i<=num_cells; i++) {
            
            hydrostat = flux[i].u2/dx[i] ; //hydrostat2 + hydrostat3 ; 
            hydrostat2 = flux[i+1].u2/dx[i];//pressure[i+1] - pressure[i];
            hydrostat3 = source[i].u2;//0.5 * (data[i].u1 + data[i+1].u1) * (phi[i+1] - phi[i]);
            
            
            //outfile<<x[i]<<'\t'<<data[i].u1<<'\t'<<data[i].u2<<'\t'<<data[i].u3<<'\t'<<fluxes[i].u1<<'\t'<<fluxes[i].u2<<'\t'<<fluxes[i].u3<<'\t'<<fluxes[i+1].u1<<'\t'<<fluxes[i+1].u2<<'\t'<<fluxes[i+1].u3<<'\t'<<pressure[i]<<'\t'<<data[i].u2/data[i].u1<<'\t'<<internal_energy[i] <<'\t'<<timesteps[i]<<'\t'<<phi[i]<<'\t'<<timesteps_cs[i]<<'\t'<<opticaldepth[i]<<'\t'<<radiative_flux[i] <<endl;
            
            //flux[i] - flux[i+1] + source[i]
            
            outfile<<x[i]<<'\t'<<data[i].u1<<'\t'<<data[i].u2<<'\t'<<data[i].u3<<'\t'<<fluxes[i].u1<<'\t'<<fluxes[i].u2<<'\t'<<fluxes[i].u3<<'\t'<<((flux[i-1].u1 - flux[i].u1)/dx[i] + source[i].u1)<<'\t'<<((flux[i-1].u2 - flux[i].u2)/dx[i] + source[i].u2)<<'\t'<<((flux[i-1].u3 - flux[i].u3)/dx[i] + source[i].u3)<<'\t'<<pressure[i]<<'\t'<<data[i].u2/data[i].u1<<'\t'<<internal_energy[i] <<'\t'<<timesteps[i]<<'\t'<<phi[i]<<'\t'<<timesteps_cs[i]<<'\t'<<opticaldepth[i]<<'\t'<<hydrostat<<'\t'<<hydrostat2<<'\t'<<hydrostat3<<'\t'<<enclosed_mass[i]<<endl;
        }
        
        //Print right ghost stuff
        outfile<<x[num_cells+1]<<'\t'<<data[num_cells+1].u1<<'\t'<<data[num_cells+1].u2<<'\t'<<data[num_cells+1].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<pressure[num_cells+1]<<'\t'<<data[num_cells+1].u2/data[num_cells+1].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<phi[num_cells+1]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
   
        cout<<"Sucessfully written file "<<filename<<" with smoothed gravity = "<<rs_at_moment<<" at time "<<globalTime<<" and dt="<<dt<<endl;
    }
    else cout << "Unable to open file" << filename << endl; 
    outfile.close();
    
 
}

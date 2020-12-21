///////////////////////////////////////////////////////////
//
//
//  io.cpp
//
// This file contains input/output routines.
//
//
//
///////////////////////////////////////////////////////////

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"


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
    this->num_opacity_datas = -1;
    
    if(debug > 0) cout<<"          In read species Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        std::vector<string> stringlist = stringsplit(line," ");

        //if(line.find(variablename) != string::npos) {
        if(stringlist[0].find("@") != string::npos) {
            
            if(std::stoi(stringlist[1]) == species_index) {
                
                found = 1;
                
                this->speciesname             = stringlist[2];
                this->mass_amu         = std::stod(stringlist[3]);
                this->degrees_of_freedom = std::stod(stringlist[4]);
                this->gamma_adiabat    = std::stod(stringlist[5]);
                this->initial_fraction = std::stod(stringlist[6]);
                this->density_excess   = std::stod(stringlist[7]);
                this->is_dust_like     = std::stod(stringlist[8]);
                
                this->opacity_data_string = stringlist[9];
                this->num_opacity_datas= std::stod(stringlist[10]);
                
                if(debug > 0)
                    cout<<"Found species called "<<speciesname<<" with a mass of "<<mass_amu<<" dof "<<degrees_of_freedom<<" gamma "<<gamma_adiabat<<" and zeta_0 = "<<initial_fraction<<endl;
                
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
    
    if(num_opacity_datas > -1) {
        
        cout<<"Starting reading opacity data. string = "<<opacity_data_string<<" num_opacity_datas = "<<num_opacity_datas<<endl;
        
        opacity_data = Eigen::MatrixXd::Zero(num_opacity_datas, 2);   //
        
        ifstream file2("inputdata/" + opacity_data_string);
        string line2;
        
        int data_count=0;
        
        while(std::getline( file2, line2 )) {
            
            std::vector<string> stringlist = stringsplit(line2," ");
            
            opacity_data(data_count, 0) = std::stod(stringlist[0]);
            opacity_data(data_count, 1) = std::stod(stringlist[1]);
            
            data_count++;
        }
        
        file2.close();
    }
    
    cout<<" IN INIT SPECIES, OPACITY_DATA for file "<<opacity_data_string<<" is "<<endl<<opacity_data<<endl;
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

    simulation_parameter<T> tmp_parameter = {"NaN",T()};
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
template simulation_parameter<bool>  read_parameter_from_file(string, string, int);
template simulation_parameter<int>   read_parameter_from_file(string, string, int, int);
template simulation_parameter<int>   read_parameter_from_file(string, string, int);
template simulation_parameter<double> read_parameter_from_file(string, string, int, double);
template simulation_parameter<double> read_parameter_from_file(string, string, int);
template simulation_parameter<char> read_parameter_from_file(string, string, int, char);
template simulation_parameter<char> read_parameter_from_file(string, string, int);
template simulation_parameter<string> read_parameter_from_file(string, string, int, string);
template simulation_parameter<string> read_parameter_from_file(string, string, int);

template simulation_parameter<Geometry> read_parameter_from_file(string, string, int, Geometry);
template simulation_parameter<Geometry> read_parameter_from_file(string, string, int);
template simulation_parameter<BoundaryType> read_parameter_from_file(string, string, int, BoundaryType);
template simulation_parameter<BoundaryType> read_parameter_from_file(string, string, int);
template simulation_parameter<IntegrationType> read_parameter_from_file(string, string, int, IntegrationType);
template simulation_parameter<IntegrationType> read_parameter_from_file(string, string, int);

//Print 2 
void c_Species::print_AOS_component_tofile(int timestepnumber) {
                           
    //
    // If we are computing winds, compute the analytic solution for outputting it
    //
    if(base->init_wind==1 && base->problem_number == 2) 
        compute_analytic_solution();
    
    // Choose a sensible default output name
    string filename ;
    {
        stringstream filenamedummy;
        string truncated_name = stringsplit(base->simname,".")[0];
        filenamedummy<<base->workingdir<<"output_"<<truncated_name<<"_"<<speciesname<<"_t"<<timestepnumber<<".dat";
        filename = filenamedummy.str() ;
    }
    /*
    filename = read_parameter_from_file<string>(base->workingdir+base->simname, "OUTPUT_FILENAME", debug-1, filename).value ;

    // Add the snap number
    {
        stringstream filenamedummy;
        filenamedummy << filename<<"_"<< speciesname << "_t" << timestepnumber << ".dat";
        filename = filenamedummy.str() ;
    }*/

    if(debug > 1)
        cout<<"Trying to open file "<<filename<<endl;
    
    ofstream outfile(filename, ios::out);
    outfile << std::setprecision(15) ;

    if (outfile.is_open())
    {
        //outfile.precision(16);
        
        double hydrostat2 = 0., hydrostat3 = 0.;
        
        //Print left ghost stuff
        outfile<<base->x_i12[0]<<'\t'<<u[0].u1<<'\t'<<u[0].u2<<'\t'<<u[0].u3<<'\t'<<flux[0].u1<<'\t'<<flux[0].u2<<'\t'<<flux[0].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[0].pres<<'\t'<<u[0].u2/u[0].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[0]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
        
        //Print the domain
        
        for(int i=1; i<=num_cells; i++) {
            
            double balance1 = ((flux[i-1].u1 * base->surf[i-1] - flux[i].u1 * base->surf[i]) / base->vol[i] + (source[i].u1 +source_pressure[i].u1));
            double balance2 = ((flux[i-1].u2 * base->surf[i-1] - flux[i].u2 * base->surf[i]) / base->vol[i] + (source[i].u2 +source_pressure[i].u2));
            double balance3 = ((flux[i-1].u3 * base->surf[i-1] - flux[i].u3 * base->surf[i]) / base->vol[i] + (source[i].u3 +source_pressure[i].u3));
            
            //hydrostat = flux[i-1].u2/base->dx[i] ; //hydrostat2 + hydrostat3 ; 
            hydrostat2 = flux[i].u2/base->dx[i];//pressure[i+1] - pressure[i];
            hydrostat3 = source[i].u2;//0.5 * (u[i].u1 + u[i+1].u1) * (phi[i+1] - phi[i]);
            
            outfile<<base->x_i12[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<balance1<<'\t'<<balance2<<'\t'<<balance3<<'\t'<<prim[i].pres<<'\t'<<u[i].u2/u[i].u1<<'\t'<<prim[i].temperature <<'\t'<<timesteps[i]<<'\t'<<base->phi[i]<<'\t'<<prim[i].sound_speed<<'\t'<<1.<<'\t'<<u_analytic[i]<<'\t'<<base->alphas_sample(i)<<'\t'<<hydrostat3<<'\t'<<base->enclosed_mass[i]-base->planet_mass<<'\t'<<(base->Jrad_FLD(i,0)/c_light*4.*pi)<<'\t'<<base->S_band(i,0)<<endl;
        }
        
        //Print right ghost stuff
        outfile<<base->x_i12[num_cells+1]<<'\t'<<u[num_cells+1].u1<<'\t'<<u[num_cells+1].u2<<'\t'<<u[num_cells+1].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[num_cells+1].pres<<'\t'<<u[num_cells+1].u2/u[num_cells+1].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[num_cells+1]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
   
        cout<<"    Sucessfully written file "<<filename<<" for species = "<<speciesname<<" at time "<<base->globalTime<<" and dt="<<base->dt<<", cflfactor = "<<base->cflfactor<<endl;
    }
    else cout << "Unable to open file" << filename << endl; 
    outfile.close();

}

void c_Sim::print_monitor(int num_steps) {

    
    
    //
    //
    // First compute the volumetric sum of all initial conserved quantities
    //
    //
    if(num_steps==0) {
        initial_monitored_quantities.u1_tot = 0.;
        initial_monitored_quantities.u2_tot = 0.;
        initial_monitored_quantities.u3_tot = 0.;
        
        for(int i=2; i<num_cells-1; i++)
            for(int s=0; s<num_species; s++) {
                initial_monitored_quantities.u1_tot += species[s].u[i].u1 * vol[i];
                initial_monitored_quantities.u2_tot += species[s].u[i].u2 * vol[i];
                initial_monitored_quantities.u3_tot += species[s].u[i].u3 * vol[i];
            }

    }
    
    //
    //
    // Open and write into the monitor file - the file monitoring quantitiess as function of time, not space
    //
    //
    
    string filename2 ;
    {
        stringstream filenamedummy;
        string truncated_name = stringsplit(simname,".")[0];
        filenamedummy<<workingdir<<"monitor_"<<truncated_name<<".dat";
        filename2 = filenamedummy.str() ;
    }
    
    //
    // Delete old file if we start a new one
    //
    if(steps==0) {
        ofstream monitor(filename2);
        monitor.close();
    }
        
    
    ofstream monitor(filename2, std::ios_base::app);
    monitor << std::setprecision(24);
    
    if(monitor.is_open()){
        
        //
        // Compute conserved quantities on entire domain over all species 
        //
        monitored_quantities.u1_tot = 0.;
        monitored_quantities.u2_tot = 0.;
        monitored_quantities.u3_tot = 0.;
        
        for(int i=2; i<num_cells-1; i++)
            for(int s=0; s<num_species; s++) {
                monitored_quantities.u1_tot += species[s].u[i].u1 * vol[i];
                monitored_quantities.u2_tot += species[s].u[i].u2 * vol[i];
                monitored_quantities.u3_tot += species[s].u[i].u3 * vol[i];
            }
        
        //
        // Print the stuff
        //
        
        //Cols 1-5
        monitor<<globalTime<<'\t'<<dt<<'\t';//<<monitored_quantities.u1_tot/initial_monitored_quantities.u1_tot<<'\t'<<monitored_quantities.u2_tot/initial_monitored_quantities.u2_tot<<'\t'<<monitored_quantities.u3_tot/initial_monitored_quantities.u3_tot<<'\t';
        
        double etot = 0;
        for(int b=0; b<num_bands; b++) 
            etot += Jrad_FLD(5,b)*4*pi/c_light; 
        monitor<<etot<<'\t';    //col 3  total radiative energy densities of all bands added up,  unit = erg/cm^3
        
        etot=0;
        for(int s=0; s<num_species; s++) 
            etot += species[s].prim[5].density*species[s].cv*species[s].prim[5].temperature;
        monitor<<etot<<'\t';     //col 4:      total internal energy densities of all species added up,  unit = erg/cm^3
        
        etot = 0;
        for(int b=0; b<num_bands; b++) 
            etot += Jrad_FLD(5,b)*4*pi/c_light; 
        monitor<<pow(c_light * etot/sigma_rad/4., 0.25)<<'\t';      //col 5: The radiative temperature 
        
        for(int s=0; s<num_species; s++)
            monitor<<species[s].prim[5].temperature<<'\t'; //col 6-...: all species temperatures 
        
        //for(int s=0; s<num_species; s++)
        //    monitor<<compute_planck_function_integral3(l_i[0],l_i[1],species[s].prim[5].temperature)<<'\t';
            
        ///(sigma_rad/pi*pow(species[s].prim[5].temperature,4.))<<'\t'; // col..? species 4 sigma_T4/c, unit = erg/cm^3
        
        
        for(int s=0; s<num_species; s++)
            monitor<<species[s].prim[5].sound_speed / dx[5]<<'\t';
            
        
        
        //Starting with col 6
        for(int b=0; b<num_bands; b++) 
            monitor<<Jrad_FLD(5,b)*4*pi/c_light<<'\t'; 
        
        
        
        //Convergence measures for radiative quantities
        for(int b=0; b<num_bands; b++)  {
            monitor<<((Jrad_FLD(5,b)-previous_monitor_J[b])/previous_monitor_J[b])<<'\t'; 
            previous_monitor_J[b] = Jrad_FLD(5,b);
        }
        
        for(int s=0; s<num_species; s++)  {
            monitor<<((sigma_rad/pi*pow(species[s].prim[5].temperature,4.)-previous_monitor_T[s])/previous_monitor_T[s])<<'\t'; 
            previous_monitor_T[s] = sigma_rad/pi*pow(species[s].prim[5].temperature,4.);
        }
        
        //Starting with col 6
        //for(int s=0; s<num_species; s++)
        //    monitor<<species[s].prim[1].speed<<'\t';
        
        for(int s=0; s<num_species; s++)
            monitor<<species[s].cv*species[s].prim[1].temperature<<'\t';
        
        monitor<<endl;
    }
    else if (num_steps == 0) cout << "WARNING: Unable to open file " << filename2 << endl;
    monitor.close();
    
    if(num_steps == -1) 
        cout<<"Finished writing monitor_file "<<filename2<<endl;
    
}

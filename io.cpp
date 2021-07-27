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
                
                //if(opacity_model == 'P') {}
                this->opacity_data_string = stringlist[9];
                //this->num_opacity_datas= std::stod(stringlist[10]);
                
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
    
    //
    // If the opacity model is 'P', physical, then we read in opacity data from files
    //
    opacity_avg_solar            = Eigen::VectorXd::Zero(num_bands_in); //num_cells * num_bands
    opacity_avg_planck           = Eigen::VectorXd::Zero(num_bands_out); //num_cells * num_bands
    opacity_avg_rosseland        = Eigen::VectorXd::Zero(num_bands_out); //num_cells * num_bands
    
    //if(num_opacity_datas > -1) 
    if(base->opacity_model == 'P' || base->opacity_model == 'M' || base->opacity_model == 'C') {
        
        cout<<"P or M or C opacity chosen & enough data to read in files. Reading file = "<<("inputdata/" +opacity_data_string)<<endl;
        //
        // Start reading opacity data
        //
        
        string opacityinputfile = "inputdata/" + opacity_data_string;
        ifstream file2( opacityinputfile);
        string line2;
        
        int num_readin_columns = 1;
        std::vector<string> stringending = stringsplit(opacityinputfile,".");
        
        if(debug > 1)
            cout<<" DEBUG LV2, namepart1 = "<<stringending[0]<<" namepart2 = "<<stringending[1]<<endl;
        
        if(stringending[1].compare("op2") == 0)
            num_readin_columns = 2;
        if(stringending[1].compare("op3") == 0)
            num_readin_columns = 3;
        
        int data_count=0;
        
        while(std::getline( file2, line2 ))
            data_count++;
        
        num_opacity_datas = data_count;
        
        file2.clear();
        file2.seekg(0, ios::beg);
        
        if(debug > 1)
            cout<<" num_readin_columns = "<<num_readin_columns<<endl;
        
        Eigen::MatrixXd file_opacity_data = Eigen::MatrixXd::Zero(num_opacity_datas, num_readin_columns + 1);   //Saves the data in format wl,opa,opa,opa as it is in the file
        
        data_count = 0;
        
        while(std::getline( file2, line2 )) {
            std::vector<string> stringlist = stringsplit(line2," ");
            //std::vector<string> stringlist = split(line2," ");
            
            file_opacity_data(data_count, 0) = std::stod(stringlist[0]); // Wavelength
            file_opacity_data(data_count, 1) = std::stod(stringlist[1]); // Opacity solar
            if(num_readin_columns >= 2) file_opacity_data(data_count, 2) = std::stod(stringlist[2]); // Opacity planck
            if(num_readin_columns >= 3) file_opacity_data(data_count, 3) = std::stod(stringlist[3]); // Opacity planck
            
            data_count++;
        }
        
        file2.close();
        //
        // Interpolate the possibly non-uniformely binned opacity data on a uniform grid for fast interpolation later
        //
        double minDeltaL = 9999999999.;
        double deltaL;
        // Determine smallest distance
        for(int j = 0; j < num_opacity_datas-1; j++) {
            deltaL    = file_opacity_data(j+1, 0) - file_opacity_data(j, 0); // Wavelength
            minDeltaL = (deltaL < minDeltaL)?deltaL:minDeltaL;
            //cout<<" found deltaL = "<<deltaL<<" at j = "<<j<<endl;
        }
        
        // Create grid with resolution as smallest distance (artificial high-res grid, that is later used for averaging)
        int num_tmp_lambdas = (file_opacity_data(num_opacity_datas-1,0) - file_opacity_data(0,0))/minDeltaL + 2;    
        opacity_data = Eigen::MatrixXd::Zero(num_tmp_lambdas, num_readin_columns + 1);   
    
        // Find high-res opacity in read-in opacity 
        for (int col = 1; col < num_readin_columns+1; col++) {
            
            int j = 0;
            for(int i = 0; i< num_tmp_lambdas; i++) {
                opacity_data(i,0) = file_opacity_data(0,0) + ((double)i) * minDeltaL;
                
                //for(int j = 0; j < num_opacity_datas-1; j++) {
                double wl = opacity_data(i,0);
                double lmin = file_opacity_data(j,   0);
                double lmax = file_opacity_data(j+1, 0);
                
                if(wl > lmax) {
                    j++;
                    lmin = file_opacity_data(j,   0);
                    lmax = file_opacity_data(j+1, 0);
                }
                
                //Interpolation on highres grid
                double m  =  std::log10(file_opacity_data(j+1,col) / file_opacity_data(j,col) ) / std::log10(lmax/lmin);
                opacity_data(i,col) = pow(10., std::log10(file_opacity_data(j,col)) + m * std::log10(wl/lmin) );
                
                //Boundaries of wl grid
                if(i==0)
                    opacity_data(i,col) = file_opacity_data(0,col);
                if(i==num_tmp_lambdas-1)
                    opacity_data(i,col) = file_opacity_data(num_opacity_datas-1,col);
                
                    //else if (wl < file_opacity_data(0, 0))
                    //    opacity_data(i,1) = file_opacity_data(0,1);
                    //else if (wl > file_opacity_data(num_tmp_lambdas-1, 0))
                    //    opacity_data(i,1) = file_opacity_data(num_tmp_lambdas-1,1);
                //}
                
            }
            
            //cout<<" ~~~~~~~~~~~~~~~ Read in opacities for col_num = "<<num_readin_columns<<endl;
        }
        
        
        //
        // Compute *the* average opacity per band (planck/rosseland/flux means have to be done on the fly)
        //
        if(debug > 1) {
            cout<<" minDeltaL = "<<minDeltaL<<" num_tmp_lambdas = "<<num_tmp_lambdas<<endl;
            cout<<"tmp_opacity_data = "<<endl<<file_opacity_data<<endl;
            //cout<<"opacity_data = "<<endl<<opacity_data<<endl;
        }
        
        for(int col = 1; col< num_readin_columns+1; col++) {
            
            int num_bands_readin;
            std::vector<double>* lgrid;
            Eigen::VectorXd*     opacity_avg;
            
            if(col==1) {
                num_bands_readin = num_bands_in;
                lgrid            = &base->l_i_in;
                opacity_avg      = &opacity_avg_solar;
            }
            else {
                num_bands_readin = num_bands_out;
                lgrid            = &base->l_i_out;
                if(col==2) opacity_avg      = &opacity_avg_planck;
                if(col==3) opacity_avg      = &opacity_avg_rosseland;
            }
        
            for(int b = 0; b< num_bands_readin; b++) {
                
                int    wlcount = 0;
                double lmin = (*lgrid)[b]; //base->l_i[b];
                double lmax = (*lgrid)[b+1]; //base->l_i[b+1];
                (*opacity_avg)[b] = 1e-10;
                
                //cout<<" DEBUG band b = "<<b<<" lmin/lmax = "<<lmin<<"/"<<lmax<<endl;
                
                for(int i = 0; i < num_tmp_lambdas; i++) {
                    double wl = opacity_data(i,0);
                    //cout<<" i = "<<i<<" wl = "<<wl<<" opa = "<<opacity_data(i,1)<<endl;
                    
                    if( wl < lmax && wl > lmin) {
                        if(opacity_data(i,col) > 1e10 || opacity_data(i,col) < 1.1e-20) {
                            //cout<<"Band opacity_data("<<i<<",1) = "<<opacity_data(i,1)<<"is faulty! wlcount ="<<wlcount<<endl;
                            //cout<<" Band b, lmax "<<lmax<<" opa_data(0,0) = "<<opacity_data(0,0)<<endl;
                            //cout<<" Band b, lmin "<<lmin<<" opa_data(-1,0) = "<<opacity_data(num_tmp_lambdas,0)<<endl;
                            
                        } else {
                            (*opacity_avg)[b]+= opacity_data(i,col); 
                            //opacity_avg(b) += opacity_data(i,col); 
                            wlcount++;
                        }
                    }
                    if(i==0) {
                        
                        //TODO: When data boundary and band boundary do not coincide
                    }
                    else if(i==num_tmp_lambdas-2) {
                        //TODO: When data boundary and band boundary do not coincide
                    }
                }
                
                if(wlcount > 0) {
                    (*opacity_avg)[b] /= (double)wlcount;
                }
                else 
                    cout<<" Band "<<b<<" has wlcount==0!"<<endl;
                    
                //
                // For bands which have no opacity data given / first and last band, we assume the nearest datapoint
                //
                if(lmax < opacity_data(0,0) ) {
                    (*opacity_avg)[b] = opacity_data(0,col);
                    cout<<" Band b, lmax "<<lmax<<"< opa_data(0,0)"<<opacity_data(0,0)<<endl;
                }
                if(lmin > opacity_data(num_tmp_lambdas-1,0) ) {
                    (*opacity_avg)[b] = opacity_data(num_tmp_lambdas-1,col);
                    cout<<" Band b, lmin "<<lmin<<"< opa_data(-1,0)"<<opacity_data(num_tmp_lambdas,0)<<endl;
                }
                if((*opacity_avg)[b] > 1e10 || (*opacity_avg)[b] < 0.) {
                    
                    cout<<"Band b = "<<b<<"is faulty! opacity = "<<(*opacity_avg)[b]<<" wlcount ="<<wlcount<<endl;
                    cout<<" Band b, lmax "<<lmax<<" opa_data(0,0) = "<<opacity_data(0,0)<<endl;
                    cout<<" Band b, lmin "<<lmin<<" opa_data(-1,0) = "<<opacity_data(num_tmp_lambdas,0)<<endl;
                }
            }
        }
        
        if(num_readin_columns < 3)
            for(int b = 0; b < num_bands_out; b++) opacity_avg_rosseland(b) = const_opacity;
        if(num_readin_columns < 2)
            for(int b = 0; b < num_bands_out; b++) opacity_avg_planck(b) = const_opacity;
            
        //cout<<"pos4"<<endl;
        //Done! Now debug plot stuff
        cout<<"Done reading in data for species = "<<species_index<<". Wavelength grids in/out = ";
        for(int b = 0; b <= num_bands_in; b++) cout<<base->l_i_in[b]<<"/";
        cout<<" ||| ";
        for(int b = 0; b <= num_bands_out; b++) cout<<base->l_i_out[b]<<"/";
        cout<<endl;
        
        cout<<"avg opacities solar = ";
        for(int b = 0; b < num_bands_in; b++) cout<<opacity_avg_solar(b)<<"/";
        cout<<endl;
        
        cout<<"avg opacities plnck = ";
        for(int b = 0; b < num_bands_out; b++) cout<<opacity_avg_planck(b)<<"/";
        cout<<endl;
        
        cout<<"avg opacities rossl = ";
        for(int b = 0; b < num_bands_out; b++) cout<<opacity_avg_rosseland(b)<<"/";
        cout<<endl;
        
    }
    else {
        for(int b = 0; b < num_bands_in; b++) opacity_avg_solar(b)      = const_opacity;
        for(int b = 0; b < num_bands_out; b++) opacity_avg_planck(b)    = const_opacity;
        for(int b = 0; b < num_bands_out; b++) opacity_avg_rosseland(b) = const_opacity;
    }
    
    
    //cout<<" IN INIT SPECIES, OPACITY_DATA for file "<<opacity_data_string<<" is "<<endl<<opacity_data<<endl;
    if(debug > 1) {
        
        char stopchar;
        cin>>stopchar;
    }
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
        if(debug > 1)
            error << "ERROR: Variable "<<variablename<<" not found in parameterfile!"<<endl;
        else
            error<<"";
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
            if(debug > 1)
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


void c_Species::print_AOS_component_tofile(int timestepnumber) {
                           
    //
    // If we are computing winds, compute the analytic solution for outputting it
    //
    if(base->init_wind==1 && base->problem_number == 2) 
        compute_analytic_solution();
    
    // 
    //
    // Step 1: Print general information file with fixed column numbers for all setups
    //
    //
    string filename ;
    {
        stringstream filenamedummy;
        string truncated_name = stringsplit(base->simname,".")[0];
        filenamedummy<<base->workingdir<<"output_"<<truncated_name<<"_"<<speciesname<<"_t"<<timestepnumber<<".dat";
        filename = filenamedummy.str() ;
    }

    if(debug > 1)
        cout<<"Trying to open file "<<filename<<endl;
    
    ofstream outfile(filename, ios::out);
    outfile << std::setprecision(15) ;

    if (outfile.is_open())
    {
        //outfile.precision(16);
        
        double hydrostat2 = 0., hydrostat3 = 0.;
        
        //Print left ghost stuff
        outfile<<base->x_i12[0]<<'\t'<<u[0].u1<<'\t'<<u[0].u2<<'\t'<<u[0].u3<<'\t'<<flux[0].u1<<'\t'<<flux[0].u2<<'\t'<<flux[0].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[0].pres<<'\t'<<u[0].u2/u[0].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[0]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
        
        //Print the domain
        base->enclosed_mass_tmp[0] = 0.;
        for(int i = 1; i <= num_cells+1; i++) {
            base->enclosed_mass_tmp[i] = base->enclosed_mass_tmp[i-1] +  4. * 3.141592 * (pow(base->x_i[i],3.)-pow(base->x_i[i-1],3.) )/3. * u[i].u1;
        }
        
        for(int i=1; i<=num_cells; i++) {
            
            double balance1 = ((flux[i-1].u1 * base->surf[i-1] - flux[i].u1 * base->surf[i]) / base->vol[i] + (source[i].u1 +source_pressure[i].u1));
            double balance2 = ((flux[i-1].u2 * base->surf[i-1] - flux[i].u2 * base->surf[i]) / base->vol[i] + (source[i].u2 +source_pressure[i].u2));
            double balance3 = ((flux[i-1].u3 * base->surf[i-1] - flux[i].u3 * base->surf[i]) / base->vol[i] + (source[i].u3 +source_pressure[i].u3));
            
            //hydrostat = flux[i-1].u2/base->dx[i] ; //hydrostat2 + hydrostat3 ; 
            hydrostat2 = flux[i].u2/base->dx[i];//pressure[i+1] - pressure[i];
            hydrostat3 = source[i].u2;//0.5 * (u[i].u1 + u[i+1].u1) * (phi[i+1] - phi[i]);
            
            double Jtot = 0.;
            double Stot = 0.;
            for(int b=0; b<num_bands_in; b++) {
                Stot += base->S_band(i,b);
            }
            for(int b=0; b<num_bands_out; b++) {
                Jtot += base->Jrad_FLD(i,b)/c_light*4.*pi;
            }
            
            outfile<<base->x_i12[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<balance1<<'\t'<<balance2<<'\t'<<balance3<<'\t'<<prim[i].pres<<'\t'<<u[i].u2/u[i].u1<<'\t'<<prim[i].temperature <<'\t'<<timesteps_cs[i]<<'\t'<<base->cflfactor/timesteps[i]<<'\t'<<prim[i].sound_speed<<'\t'<<timesteps_de[i]<<'\t'<<u_analytic[i]<<'\t'<<base->alphas_sample(i)<<'\t'<<base->phi[i]<<'\t'<<base->enclosed_mass_tmp[i]<<'\t'<<Jtot<<'\t'<<Stot<<'\t'<<base->friction_sample(i)<<endl;
        } //prim[i].sound_speed
        
        //Print right ghost stuff
        outfile<<base->x_i12[num_cells+1]<<'\t'<<u[num_cells+1].u1<<'\t'<<u[num_cells+1].u2<<'\t'<<u[num_cells+1].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[num_cells+1].pres<<'\t'<<u[num_cells+1].u2/u[num_cells+1].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[num_cells+1]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
   
        cout<<"    Sucessfully written file "<<filename<<" for species = "<<speciesname<<" t = "<<base->globalTime<<" dt = "<<base->dt<<", cfl = "<<base->cflfactor<<" steps = "<<base->steps<<endl;
    }
    else cout << "Unable to open file" << filename << endl; 
    outfile.close();
    
    

}

void c_Sim::print_diagnostic_file(int outputnumber) {
    
    //
    //
    // Step 2: Print radiation + species diagnostics file. Column numbers in this file will vary depending on setup
    //
    //
    
    //Only print this file once, as it contains all the info 
    
        
        string filenameDiagnostic ;
        {
            stringstream filenamedummy;
            string truncated_name = stringsplit(simname,".")[0];
            filenamedummy<<workingdir<<"diagnostic_"<<truncated_name<<"_t"<<outputnumber<<".dat";
            filenameDiagnostic = filenamedummy.str() ;
        }

        if(debug > 1)
            cout<<"Trying to open file "<<filenameDiagnostic<<endl;
        
        ofstream outfileDiagnostic(filenameDiagnostic, ios::out);
        outfileDiagnostic << std::setprecision(15) ;

        if (outfileDiagnostic.is_open())
        {
            
            //outfileDiagnostic<<""<<endl;
            
            for(int i=1; i<=num_cells; i++) {
                
                //
                //Band-by-band info
                // 
                double Jtot = 0.;
                double Stot = 0.;
                for(int b=0; b<num_bands_in; b++) {
                    Stot += S_band(i,b);
                }
                for(int b=0; b<num_bands_out; b++) {
                    Jtot += Jrad_FLD(i,b)/c_light*4.*pi;
                }
                outfileDiagnostic<<x_i12[i]<<'\t'<<Jtot<<'\t'<<Stot; //Rows 1 2 3
                
                //Radiation field
                for(int b=0; b<num_bands_out; b++) {
                    outfileDiagnostic<<'\t'<<Jrad_FLD(i,b);
                }
                
                //Solar radiation field
                for(int b=0; b<num_bands_in; b++) {
                    outfileDiagnostic<<'\t'<<S_band(i,b);
                }
                
                //Solar heating profile
                for(int b=0; b<num_bands_in; b++) {
                    outfileDiagnostic<<'\t'<<dS_band(i,b);
                }
                
                //Total opacity from all species
                for(int b=0; b<num_bands_out; b++) {
                    outfileDiagnostic<<'\t'<<total_opacity(i,b);
                }
                
                //Total opacity from all species
                for(int b=0; b<num_bands_out; b++) {
                    outfileDiagnostic<<'\t'<<cell_optical_depth(i,b);
                }
                
                //Optical depths for solar radiation
                for(int b=0; b<num_bands_in; b++) {
                    //outfileDiagnostic<<'\t'<<radial_optical_depth(i,b);
                    outfileDiagnostic<<'\t'<<radial_optical_depth_twotemp(i,b);
                    
                }
                
                //
                //Species-by-species and band opacities
                //
                for(int b=0; b<num_bands_in; b++) {
                    for(int s=0; s<num_species; s++) {
                        outfileDiagnostic<<'\t'<<const_opacity_solar_factor*species[s].opacity_twotemp(i,b);
                    }
                }
                
                for(int b=0; b<num_bands_out; b++) {
                    for(int s=0; s<num_species; s++) {
                        outfileDiagnostic<<'\t'<<species[s].opacity(i,b);
                        outfileDiagnostic<<'\t'<<species[s].opacity_planck(i,b);
                    }
                }
                
                //Radiative and species temperatures
                outfileDiagnostic<<'\t'<<pow(Jtot*pi/sigma_rad*c_light/4./3.141, 0.25);      //col 5: The radiative temperature 
                
                for(int s=0; s<num_species; s++)
                    outfileDiagnostic<<'\t'<<species[s].prim[i].temperature; //col 6-...: all species temperatures 
                
                //Pressure
                double tempP = 0.;
                for(int s=0; s<num_species; s++)
                        tempP += species[s].prim[i].pres;
                outfileDiagnostic<<'\t'<<tempP;
                
                //Fluxes
                for(int b=0; b<num_bands_out; b++) {
                    
                    double dx      = (x_i12[i+1]-x_i12[i]) ;
                    
                    double rhokr   = max(2.*(total_opacity(i,b)*total_opacity(i+1,b))/(total_opacity(i,b) + total_opacity(i+1,b)), 4./3./dx );
                           rhokr   = min( 0.5*( total_opacity(i,b) + total_opacity(i+1,b)) , rhokr);
                    double tau_inv = 0.5 / (dx * rhokr) ;
                
                    //double tau_inv = 0.5 / (dx * (total_opacity(i,b) + total_opacity(i+1,b))) ;
                    double R       = 2 * tau_inv * std::abs(Jrad_FLD(i+1,b) - Jrad_FLD(i,b)) / (Jrad_FLD(i+1,b) + Jrad_FLD(i, b) + 1e-300) ;
                    double flux_limiter = 0;
                    if (R <= 2)
                        flux_limiter = 2 / (3 + std::sqrt(9 + 10*R*R)) ;
                    else 
                        flux_limiter = 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
                    
                    double D       = surf[i] * flux_limiter * tau_inv;
                    double flux    = - 4. * pi * D * (Jrad_FLD(i+1,b) - Jrad_FLD(i,b));

                    outfileDiagnostic<<'\t'<<flux;
                }
                outfileDiagnostic<<endl;
            }
            
            
            
        }else cout << "Unable to open file" << filenameDiagnostic << endl; 
        outfileDiagnostic.close();
    
    
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
        
        
        for(int s=0; s<num_species; s++) 
            monitor<<species[s].prim[5].speed<<'\t';
             //col 4:      total internal energy densities of all species added up,  unit = erg/cm^3
        
        double etot = 0;
        for(int b=0; b<num_bands_out; b++) 
            etot += Jrad_FLD(5,b)*4*pi/c_light; 
        monitor<<etot<<'\t';    //col 3  total radiative energy densities of all bands added up,  unit = erg/cm^3
        
        etot=0;
        for(int s=0; s<num_species; s++) 
            etot += species[s].prim[5].density*species[s].cv*species[s].prim[5].temperature;
        monitor<<etot<<'\t';     //col 4:      total internal energy densities of all species added up,  unit = erg/cm^3
        
        etot = 0;
        for(int b=0; b<num_bands_out; b++) 
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
        for(int b=0; b<num_bands_out; b++) 
            monitor<<Jrad_FLD(5,b)*4*pi/c_light<<'\t'; 
        
        
        //Convergence measures for radiative quantities
        for(int b=0; b<num_bands_out; b++)  {
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

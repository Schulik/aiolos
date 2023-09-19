
/**
 * io.cpp
 * 
 * This file contains input/output routines.
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"


/**
 * Reads the species data such as name, mass, adiabatic index, charge, initial density in ghost cell, opacity files etc.
 * 
 * Does a number of convoluted operations to allow for different file formats. Contains also the code to interpolate opacities on-the-fly. 
 * This happens assuming that the *.opa, or *.op2 or *.op3 files have an arbitrary wavelength spacing, and opacity points given at any wavelength are exact. 
 * The mean opacities are then computed as integrals over the bands chosen by the user. Note that those means are not Planck-means, but purely data-driven averages with the 
 * weighting function per wavelength being 1. If correct Planck, Solar and Rosseland means need to be used, then those have to be already contained in the opa files, or tables need to be specified.
 * 
 * @param[in] filename Species file to read from
 * @param[in] species_index This species number
 */
int c_Species::read_species_data(string filename, int species_index) {
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //
    // Read and interpret *.spc first
    //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    
    //
    // If the opacity model is 'P', physical, then we read in opacity data from files
    //
    opacity_avg_solar            = Eigen::VectorXd::Zero(num_bands_in); //num_cells * num_bands
    opacity_avg_planck           = Eigen::VectorXd::Zero(num_bands_out); //num_cells * num_bands
    opacity_avg_rosseland        = Eigen::VectorXd::Zero(num_bands_out); //num_cells * num_bands
    
    ifstream file(filename);
    string line;
    if(!file) {
        cout<<"Couldnt open species file "<<filename<<"!"<<endl;
        
    }
    
    int found = 0;
    double temp_static_charge = 0.;
    this->num_opacity_datas = -1;
    
    if(debug > 0) cout<<"          In read species Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        std::vector<string> stringlist = stringsplit(line," ");

        if(stringlist[0].find("@") != string::npos) {
            
            if(std::stoi(stringlist[1]) == species_index) {
                
                found = 1;
                
                this->speciesname             = stringlist[2];
                this->mass_amu         = std::stod(stringlist[3]);
                this->degrees_of_freedom = std::stod(stringlist[4]);
                temp_static_charge       = std::stod(stringlist[5]);
                this->initial_fraction = std::stod(stringlist[6]);
                this->density_excess   = std::stod(stringlist[7]);
                this->is_dust_like     = std::stod(stringlist[8]);
                
                this->opacity_data_string = stringlist[9];
                
                if(stringlist.size() == 11) //Assume that additional to an opacity table in opacity_data_string a solar opacity wavelength file has been provided
                    this->opacity_corrk_string = stringlist[10];
                
                this->inv_mass = 1./(mass_amu*amu);
                
                if(debug > 0)
                    cout<<"Found species called "<<speciesname<<" with a mass of "<<mass_amu<<" dof "<<degrees_of_freedom<<" gamma "<<gamma_adiabat<<" and initial_fraction = "<<initial_fraction<<endl;
                
            }

        }
    
    }
    
    if(std::abs(temp_static_charge) < 0.1) {
        this->static_charge = 0;
    }
    else {
        
        this->static_charge = std::floor(temp_static_charge);
    }
    
    if(found == 0) {
            if(debug > 0)
                cout<<"WARNING: Species number "<<species_index<<" not found in parameterfile!"<<endl;
    }
    if(found > 1) {
            cout<<"WARNING: Species number "<<species_index<<" defined more than once in parameterfile!"<<endl;
    }
    
    file.close();
    
    if(debug > 0) cout<<"         Species readin finished."<<endl;
    
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //
    // Read and interpret *.opa, *.op2, *.op3, and *.aiopa files
    //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    
    if(base->opacity_model == 'T' || base->opacity_model == 'K') { 
        //
        // Tabulated data from opacit averages
        //
        
        string opacityinputfile = "inputdata/" + opacity_data_string;
        std::vector<string> stringending = stringsplit(opacityinputfile,".");
        
        if(stringending[1].compare("aiopa") == 0) {
            cout<<"Tabulated opacities chosen, reading file = "<<opacityinputfile<<endl;
        }
        else
            cout<<"Invalid file = "<<opacityinputfile<<" should have format *.aitab ."<<endl;
        
        read_opacity_table(opacityinputfile);
    }
    
    if(base->opacity_model == 'P' || base->opacity_model == 'M' || base->opacity_model == 'C' || base->opacity_model == 'K' ) {
        
        
        if( base->opacity_model == 'K' ) //Run the routine with one argument further on in the species list 
            opacity_data_string = opacity_corrk_string;
        
        cout<<"P or M or C or K opacity chosen & enough data to read in files. Reading file = "<<"inputdata/"<<opacity_data_string<<endl;
        //
        // Start reading opacity data
        //
        
        string opacityinputfile = "inputdata/" + opacity_data_string;
        ifstream file2( opacityinputfile);
        string line2;
        
        int num_readin_columns = 1;
        std::vector<string> stringending = stringsplit(opacityinputfile,".");
        
        if(debug > 0)
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
        
        if(debug > 0)
            cout<<" num_readin_columns = "<<num_readin_columns<<" data_count = "<<data_count<<endl;
        
        Eigen::MatrixXd file_opacity_data = Eigen::MatrixXd::Zero(num_opacity_datas, num_readin_columns + 1);   //Saves the data in format wl,opa,opa,opa as it is in the file
        
        data_count = 0;
        
        while(std::getline( file2, line2 )) {
            std::vector<string> stringlist = stringsplit(line2," ");
            
            file_opacity_data(data_count, 0) = std::stod(stringlist[0]); // Wavelength
            file_opacity_data(data_count, 1) = std::stod(stringlist[1]); // Opacity solar
            if(num_readin_columns >= 2) file_opacity_data(data_count, 2) = std::stod(stringlist[2]); // Opacity planck
            if(num_readin_columns >= 3) file_opacity_data(data_count, 3) = std::stod(stringlist[3]); // Opacity planck
            
            if(debug > 1) {
                cout<< "file_opacity_data = "<<file_opacity_data(data_count, 0)<<" "<<file_opacity_data(data_count, 1);
                if(num_readin_columns >= 2) cout<<" "<<file_opacity_data(data_count, 2);
                if(num_readin_columns >= 3) cout<<" "<<file_opacity_data(data_count, 3);
                cout<<endl;
            }
            data_count++;
        }
        
         if(debug > 0)
            cout<<"Before fileclose"<<endl;
        
        file2.close();
        
        if(debug > 0)
            cout<<"After fileclose"<<endl;
        //
        // Interpolate the possibly non-uniformely binned opacity data on a uniform grid for fast interpolation later
        //
        double minDeltaL = 9999999999.;
        double deltaL;
        // Determine smallest distance
        for(int j = 0; j < num_opacity_datas-1; j++) {
            deltaL    = file_opacity_data(j+1, 0) - file_opacity_data(j, 0); // Wavelength
            minDeltaL = (deltaL < minDeltaL)?deltaL:minDeltaL;
            
            if(debug > 1)
                cout<<" found deltaL = "<<deltaL<<" at j = "<<j<<endl;
        }
        
        
        // Create grid with resolution as smallest distance (artificial high-res grid, that is later used for averaging)
        int num_tmp_lambdas = (file_opacity_data(num_opacity_datas-1,0) - file_opacity_data(0,0))/minDeltaL + 2;    
        
        if(debug>0)
            cout<<" found mindeltaL = "<<minDeltaL<<" num_tmp_lambdas "<<num_tmp_lambdas<<endl;
        //if(debug>1)
        opacity_data = Eigen::MatrixXd::Zero(num_tmp_lambdas, num_readin_columns + 1);   
    
        // Find high-res opacity in read-in opacity 
        for (int col = 1; col < num_readin_columns+1; col++) {
            
            if(debug > 1)
                cout<<" Doing column "<<col<<endl;
            
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
                
                //Boundaries of wl grid
                if(i==0)
                    opacity_data(i,col) = file_opacity_data(0,col);
                else if(i==num_tmp_lambdas-1)
                    opacity_data(i,col) = file_opacity_data(num_opacity_datas-1,col);
                else {
                    
                    
                    //Interpolation on highres grid
                    double m  =  std::log10(file_opacity_data(j+1,col) / file_opacity_data(j,col) ) / std::log10(lmax/lmin);
                    opacity_data(i,col) = pow(10., std::log10(file_opacity_data(j,col)) + m * std::log10(wl/lmin) );
                }
                
            }
            
            //cout<<" ~~~~~~~~~~~~~~~ Read in opacities for col_num = "<<num_readin_columns<<endl;
        }
        
        
        //
        // Compute *the* average opacity per band (planck/rosseland/flux means have to be done on the fly)
        //
        if(debug > 1) {
            cout<<" minDeltaL = "<<minDeltaL<<" num_tmp_lambdas = "<<num_tmp_lambdas<<endl;
            if(debug > 2)
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
                (*opacity_avg)[b] = base->minimum_opacity;
                
                //cout<<" DEBUG band b = "<<b<<" lmin/lmax = "<<lmin<<"/"<<lmax<<endl;
                
                double wl;
                for(int i = 0; i < num_tmp_lambdas; i++) {
                    double wl = opacity_data(i,0);
                    
                    //cout<<" i = "<<i<<" wl = "<<wl<<" opa = "<<opacity_data(i,1)<<endl;
                    
                    if( wl < lmax && wl > lmin) {
                        
                        (*opacity_avg)[b]+= opacity_data(i,col); 
                        wlcount++;
                    }
                    if(i==0) {
                        
                        //TODO: When data boundary and band boundary do not coincide
                    }
                    else if(i==num_tmp_lambdas-2) {
                        //TODO: When data boundary and band boundary do not coincide
                    }
                }
                
                if(wlcount > 1) {
                    (*opacity_avg)[b] /= (double)wlcount;
                }
                else 
                    cout<<" Band "<<b<<" has wlcount==0! lmin/wl/lmax//opacity_data(0,0/1); = "<<lmin<<"/"<<wl<<"/"<<lmax<<"//"<<opacity_data(0,0)<<"/"<<opacity_data(0,1)<<" wlcount = "<<wlcount<<endl;
                
                if(debug > 1)
                    cout<<" DEBUG. Band "<<b<<" has wlcount=="<<wlcount<<". lmin/wl/lmax//opacity_data(0,0/1); = "<<lmin<<"/"<<wl<<"/"<<lmax<<"//"<<opacity_data(0,0)<<"/"<<opacity_data(0,1)<<" opa_found = "<<(*opacity_avg)[b]<<endl;
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
                if((*opacity_avg)[b] < 0.) {
                    
                    cout<<"Band b = "<<b<<"is faulty! opacity = "<<(*opacity_avg)[b]<<" wlcount ="<<wlcount<<endl;
                    cout<<" Band b, lmin "<<lmin<<" opa_data(-1,0) = "<<opacity_data(num_tmp_lambdas,0)<<endl;
                    cout<<" Band b, lmax "<<lmax<<" opa_data(0,0) = "<<opacity_data(0,0)<<endl;
                }
            }
        }
        
        if(debug > 1) {
            cout<<"Pos before additional columns."<<endl;
        }
        
        if(num_readin_columns < 3) 
            for(int b = 0; b < num_bands_out; b++)  {
                
                opacity_avg_rosseland(b) = const_opacity;
                if(debug> 1)
                    cout<<" DEBUG b_out = "<<b<<" avg_rosseland = "<<opacity_avg_rosseland(b)<<" const_opa = "<<const_opacity<<endl;
                
            }
        
        if(num_readin_columns < 2)
            for(int b = 0; b < num_bands_out; b++)  {
                opacity_avg_planck(b) = const_opacity;
                if(debug> 1)
                    cout<<" DEBUG b_out = "<<b<<" avg_planck = "<<opacity_avg_rosseland(b)<<" const_opa = "<<const_opacity<<endl;
                
            }
            
        if(debug > 1) {
            cout<<"Pos after additional columns."<<endl;
        }
            
        for(int b = 0; b < num_bands_in; b++)  if(std::isnan(opacity_avg_solar(b) )) opacity_avg_solar(b) = base->minimum_opacity;
        for(int b = 0; b < num_bands_out; b++) if(std::isnan(opacity_avg_planck(b) )) opacity_avg_planck(b) = base->minimum_opacity;
        for(int b = 0; b < num_bands_out; b++) if(std::isnan(opacity_avg_rosseland(b) )) opacity_avg_rosseland(b) = base->minimum_opacity;
        
        for(int b = 0; b < num_bands_in; b++)  if((opacity_avg_solar(b)<1e-10 )) opacity_avg_solar(b) = base->minimum_opacity;
        for(int b = 0; b < num_bands_out; b++) if((opacity_avg_planck(b)<1e-10 )) opacity_avg_planck(b) = base->minimum_opacity;
        for(int b = 0; b < num_bands_out; b++) if((opacity_avg_rosseland(b)<1e-10 )) opacity_avg_rosseland(b) = base->minimum_opacity;
        
        if(this->this_species_index ==7) {
            if(num_bands_in >= 10) 
                opacity_avg_solar(10) = base->minimum_opacity;
        }
        
        //Done! Now plot debug stuff so that we're sure the opacities do what they should.
        cout<<"Done reading in data for species = "<<species_index<<". Wavelength grids in/out = ";
        for(int b = 0; b <= num_bands_in; b++) cout<<base->l_i_in[b]<<"/";
        cout<<" ||| ";
        for(int b = 0; b <= num_bands_out; b++) cout<<base->l_i_out[b]<<"/";
        cout<<endl<<endl;
        
        cout<<"        avg opacities solar = "<<endl;
        for(int b = 0; b < num_bands_in; b++) cout<<"<"<<base->l_i_in[b]<<" - "<<base->l_i_in[b+1]<<" mum > "<<opacity_avg_solar(b)<<endl; //"/";
        cout<<endl;
        
        cout<<"        avg opacities plnck = ";
        for(int b = 0; b < num_bands_out; b++) cout<<opacity_avg_planck(b)<<"/";
        cout<<endl;
        
        cout<<"        avg opacities rossl = ";
        for(int b = 0; b < num_bands_out; b++) cout<<opacity_avg_rosseland(b)<<"/";
        cout<<endl;
        
    }
     else {
         //No other opacity mode needed for now. Do nothing. Table reading mode has been moved to lines 120ff.
    }
    
    return 0;
}
 
/**
 * Read in pressure and temperature dependent opacity tables from *.aiopa files.
 * 
 * Aiopa files currently have the hardcoded format: Arbitrary number of P-T grids for Irradiation means, One P-T grid for Planck-mean, One P-T grid for outgoing Rosseland mean band value.
 * 
 * @param[in] tablename Filename of the table file.
 */
void c_Species::read_opacity_table(string tablename) {
    
    int num_bands_out = base->num_bands_out;
    int num_bands_in  = base->num_bands_in;
    
    cout<<"Attempting to open "<<tablename<<"!"<<endl;
    
    ifstream file(tablename);
    string line;
    if(!file) {
        cout<<"Couldnt open opacity table "<<tablename<<"!"<<endl;
    }
    //simulation_parameter tmp_parameter = {"NaN",0,0.,0,"NaN"};
    this->num_opacity_datas = -1;
    
    if(debug >= 1) cout<<"          In read Opacities Pos1. Trying to 2*read num_bands_out + num_bands_in ="<<2*num_bands_out+num_bands_in<<" opacity blocks."<<endl;
    if(num_bands_in == num_bands_out)
        cout<<"WARNING: num_bands_in == num_bands_out, this is the default. Check if you really want this. Might break opacity reading."<<endl;
    
    //Get the size of the table first
    //
    std::getline( file, line );
    //cout<<"Got the line!"<<endl;
    std::vector<string> stringlist = stringsplit(line," ");
    
    //cout<<"String is split!"<<stringlist[0]<<endl;
    opa_pgrid_size = std::stod(stringlist[0]);
    opa_tgrid_size = std::stod(stringlist[1]);
    
    //cout<<" Identified grid size ="<<opa_pgrid_size<<" / "<<opa_tgrid_size<<endl;
    
    if(opa_pgrid_size <= 0 || opa_tgrid_size > 9999) {
        cout<<" ERROR with pressure grid size = "<<opa_pgrid_size<<" in read_opacity_table!"<<endl;
    }
    if(opa_tgrid_size <= 0 || opa_tgrid_size > 9999) {
        cout<<" ERROR with temperature grid size = "<<opa_tgrid_size<<" in read_opacity_table!"<<endl;
    }
    
    //Now allocate memory and read data
    //
    opa_pgrid = Eigen::VectorXd::Zero(opa_pgrid_size);
    opa_tgrid = Eigen::VectorXd::Zero(opa_tgrid_size);
    
    opa_pgrid_log = Eigen::VectorXd::Zero(opa_pgrid_size);
    opa_tgrid_log = Eigen::VectorXd::Zero(opa_tgrid_size);
    
    opa_grid_solar     = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_in);
    opa_grid_rosseland = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_out);
    opa_grid_planck    = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_out);
    
    opa_grid_solar_log     = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_in);
    opa_grid_rosseland_log = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_out);
    opa_grid_planck_log    = Eigen::VectorXd::Zero(opa_pgrid_size * opa_tgrid_size * num_bands_out);
    
    //cout<<" After memory allocation. Sizes are "<<opa_pgrid_size * opa_tgrid_size * num_bands_in<<" and "<<opa_pgrid_size * opa_tgrid_size * num_bands_out<<endl;
    
    //One empty line
    std::getline( file, line );
    
    //Read pgrid
    //
    std::getline( file, line );
    stringlist = stringsplit(line," ");
    for(int i=0; i<opa_pgrid_size; i++) {
        opa_pgrid(i) = std::stod(stringlist[i]);
        opa_pgrid_log(i) = std::log10(opa_pgrid(i));
        //cout<<" opa_pgrid = "<<opa_pgrid(i)<<endl;
    }
    
    //Read tgrid
    //
    std::getline( file, line );
    stringlist = stringsplit(line," ");
    for(int i=0; i<opa_tgrid_size; i++) {
        opa_tgrid(i) = std::stod(stringlist[i]);
        opa_tgrid_log(i) = std::log10(opa_tgrid(i));
        //cout<<" opa_tgrid = "<<opa_tgrid(i)<<endl;
    }
    
    // Planck two-temperature mean-opacity grid
    //
    std::getline( file, line );
    for(int bi = 0; bi<1; bi++) { //for(int bi = 0; bi<num_bands_in; bi++) { TODO: REMOVE THIS LATER!
        
        //std::getline( file, line );
        for(int i=0; i<opa_pgrid_size; i++) {
            std::getline( file, line );
            stringlist = stringsplit(line," ");
            
            for(int j=0; j<opa_tgrid_size; j++) {
                //cout<<" opas_ = "<<std::stod(stringlist[j])<<" assigned to "<< j + i * opa_tgrid_size + bi * opa_pgrid_size * opa_tgrid_size<<" list pos "<<j<<endl;
                opa_grid_solar( j + i * opa_tgrid_size + bi * opa_pgrid_size * opa_tgrid_size )     = std::stod(stringlist[j]);
                opa_grid_solar_log( j + i * opa_tgrid_size + bi * opa_pgrid_size * opa_tgrid_size ) = std::log10(std::stod(stringlist[j]));
            }
        }
        
    }
    cout<<" finshed opas "<<endl;
    
    // Planck single-temperature mean-opacity grid
    //
    std::getline( file, line );
    for(int bo = 0; bo<num_bands_out; bo++) {
        
        //std::getline( file, line );
        for(int i=0; i<opa_pgrid_size; i++) {
            std::getline( file, line );
            stringlist = stringsplit(line," ");
            
            for(int j=0; j<opa_tgrid_size; j++) {
                //cout<<" opap_ = "<<std::stod(stringlist[j])<<" assigned to "<< j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size<<" list pos "<<j<<endl;
                opa_grid_planck( j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size )     = std::stod(stringlist[j]);
                opa_grid_planck_log( j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size ) = std::log10(std::stod(stringlist[j]));
            }
        }
    }
    
    cout<<" finshed opap "<<endl;
    // Rosseland mean-opacity grid
    //
    std::getline( file, line );
    for(int bo = 0; bo<num_bands_out; bo++) {
        
        //std::getline( file, line );
        for(int i=0; i<opa_pgrid_size; i++) {
            std::getline( file, line );
            stringlist = stringsplit(line," ");
            
            for(int j=0; j<opa_tgrid_size; j++) {
                //cout<<" opar_ = "<<std::stod(stringlist[j])<<" assigned to "<< j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size<<endl;
                opa_grid_rosseland( j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size )     = std::stod(stringlist[j]);
                opa_grid_rosseland_log( j + i * opa_tgrid_size + bo * opa_pgrid_size * opa_tgrid_size ) = std::log10(std::stod(stringlist[j]));
            }
        }
        
    }
    
    cout<<" Successfully finished reading opacity table. "<<endl;
    //cout<<"Showing read data: pgrid / tgrid = "<<endl<<opa_pgrid<<endl<<"///"<<endl<<opa_tgrid<<endl;
    //cout<<"Showing read data: solar opas = "<<endl<<opa_grid_solar<<endl;
    //cout<<"Showing read data: planck opas = "<<endl<<opa_grid_planck<<endl;
    //cout<<"Showing read data: planck rosseland = "<<endl<<opa_grid_rosseland<<endl;
    
    //char a;
    //cin>>a;
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

/**
 * Read parameter: Looks for the value of one parameter in a file
 * 
 * @param[in] filename The filename that should contain all parameters and their values
 * @param[in] variablename the name of one individual parameter in the file
 * @param[in] debug Display debug information.
 * @return The value of variablename in filename as given by the template type. User is responsible of type safety here.
 */
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


/**
 * Read parameter: Looks for the value of one parameter in a file
 * 
 * @param[in] filename The filename that should contain all parameters and their values
 * @param[in] variablename the name of one individual parameter in the file
 * @param[in] debug Display debug information.
 * @param[in] default_ Default value for variablename, which is returned if variablename is not found in filename.
 * @return The value of variablename in filename as given by the template type. User is responsible of type safety here.
 */
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


/**
 * Write the main output file. Each species generates their own output file of a fixed column number (as opposed to diagnostic files, which change their column numbers depending on simulation setup, i.e. num species, bands in, bands out).
 * 
 * File contains mostly hydrodynamic and thermodynamic information, cfl timestep limitations, gravity, heating and cooling functions.
 * Users can adjust whether they want ghost cell information in their files or not.
 * More details in cheat_sheet.ods in the main aiolos github.
 * 
 * @param[in] timestepnumber Timestamp index on the file (e.g. 0, 1, 2... -1)
 */
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
        
        //double hydrostat2 = 0., hydrostat3 = 0.;
        
        //Print left ghost stuff
        outfile<<base->x_i12[0]<<'\t'<<u[0].u1<<'\t'<<u[0].u2<<'\t'<<u[0].u3<<'\t'<<flux[0].u1<<'\t'<<flux[0].u2<<'\t'<<flux[0].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[0].pres<<'\t'<<u[0].u2/u[0].u1<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<base->phi[0]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
        
        //Print the domain
        base->enclosed_mass_tmp[0] = 0.;
        for(int i = 2; i <= num_cells; i++) {
            base->enclosed_mass_tmp[i] = base->enclosed_mass_tmp[i-1] +  4. * 3.141592 * (pow(base->x_i[i],3.)-pow(base->x_i[i-1],3.) )/3. * u[i].u1;
        }
        
        for(int i=1; i<= num_cells; i++) {
            
            double balance1 = ((flux[i-1].u1 * base->surf[i-1] - flux[i].u1 * base->surf[i]) / base->vol[i] + (source[i].u1 +source_pressure[i].u1));
            double balance2 = ((flux[i-1].u2 * base->surf[i-1] - flux[i].u2 * base->surf[i]) / base->vol[i] + (source[i].u2 +source_pressure[i].u2));
            double balance3 = ((flux[i-1].u3 * base->surf[i-1] - flux[i].u3 * base->surf[i]) / base->vol[i] + (source[i].u3 +source_pressure[i].u3));
            
            //hydrostat = flux[i-1].u2/base->dx[i] ; //hydrostat2 + hydrostat3 ; 
            //hydrostat2 = flux[i].u2/base->dx[i];//pressure[i+1] - pressure[i];
            //hydrostat3 = source[i].u2;//0.5 * (u[i].u1 + u[i+1].u1) * (phi[i+1] - phi[i]);
            
            double Jtot = 0.;
            double Stot = 0.;
            for(int b=0; b<num_bands_in; b++) {
                Stot += base->S_band(i,b);
            }
            for(int b=0; b<num_bands_out; b++) {
                Jtot += base->Jrad_FLD(i,b)/c_light*4.*pi;
            }
            
            //outfile<<base->x_i12[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<balance1<<'\t'<<balance2<<'\t'<<balance3<<'\t'<<prim[i].pres<<'\t'<<u[i].u2/u[i].u1<<'\t'<<prim[i].temperature <<'\t'<<timesteps_cs[i]<<'\t'<<base->cflfactor/timesteps[i]<<'\t'<<prim[i].sound_speed<<'\t'<<timesteps_de[i]<<'\t'<<u_analytic[i]<<'\t'<<base->alphas_sample(i)<<'\t'<<base->phi[i]<<'\t'<<base->enclosed_mass_tmp[i]<<'\t'<<Jtot<<'\t'<<Stot<<'\t'<<base->friction_sample(i)<<endl;
            
            outfile<<base->x_i12[i]<<'\t'<<u[i].u1<<'\t'<<u[i].u2<<'\t'<<u[i].u3<<'\t'<<flux[i].u1<<'\t'<<flux[i].u2<<'\t'<<flux[i].u3<<'\t'<<balance1<<'\t'<<balance2<<'\t'<<balance3<<'\t'<<prim[i].pres<<'\t'<<u[i].u2/u[i].u1<<'\t'<<prim[i].temperature <<'\t'<<timesteps_cs[i]<<'\t'<<base->cflfactor/timesteps[i]<<'\t'<<prim[i].sound_speed<<'\t'<<timesteps_de[i]<<'\t'<<u_analytic[i]<<'\t'<<base->alphas_sample(i)<<'\t'<<phi_s[i]<<'\t'<<base->enclosed_mass_tmp[i]<<'\t'<<-dG(i)+0.*dGdT(i)*prim[i].temperature <<'\t'<<dS(i)<<'\t'<<base->cell_optical_depth(i,0)<<endl;
        } 
        //cout<<"bonus cooling info at output time: "<<-dG(num_cells/2)+dGdT(num_cells/2)*prim[num_cells/2].temperature<<endl;
        //Print right ghost stuff
        outfile<<base->x_i12[num_cells+1]<<'\t'<<u[num_cells+1].u1<<'\t'<<u[num_cells+1].u2<<'\t'<<u[num_cells+1].u3<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<prim[num_cells+1].pres<<'\t'<<u[num_cells+1].u2/u[num_cells+1].u1<<'\t'<<prim[num_cells+1].temperature<<'\t'<<'-'<<'\t'<<base->phi[num_cells+1]<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<'\t'<<'-'<<endl;
   
        cout<<"    Sucessfully written file "<<filename<<" for species = "<<speciesname<<" t = "<<base->globalTime<<" dt = "<<base->dt<<", cfl = "<<base->cflfactor<<" steps = "<<base->steps<<endl;
    }
    else cout << "Unable to open file" << filename << endl; 
    outfile.close();

}

/**
 * Write the diagnostic file. Has a setup-dependent column number  i.e. depends on num species, bands in, bands out, as opposed to the main output files which have a fixed number of columns.
 * 
 * File contains information about all radiative bands, incoming radiation, total heating function, emitted radiation, optical depths per cell, integrated optical depths, opacities, radiative temperatures per band, outgoing fluxes per band, convective fluxes. More details in cheat_sheet.ods in the main aiolos github.
 * 
 * @param[in] outputnumber Timestamp index on the file (e.g. 0, 1, 2... -1)
 */
void c_Sim::print_diagnostic_file(int outputnumber) {
    
    //
    //
    // Print radiation + species diagnostics file. Column numbers in this file will vary depending on setup
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
                    Jtot += Jrad_FLD(i,b)/c_light*4.*pi;  //This actually outputs the energy density E_rad from c E_rad = 4pi J
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
                    outfileDiagnostic<<'\t'<<dS_band(i,b) + dS_band_special(i,b);
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
                        outfileDiagnostic<<'\t'<<species[s].opacity_twotemp(i,b);
                    }
                }
                
                for(int b=0; b<num_bands_out; b++) {
                    
                    for(int s=0; s<num_species; s++) {
                        outfileDiagnostic<<'\t'<<species[s].opacity_planck(i,b);
                    }
                    for(int s=0; s<num_species; s++) {
                        outfileDiagnostic<<'\t'<<species[s].opacity(i,b);
                    }
                }
                
                //Radiative and species temperatures
                outfileDiagnostic<<'\t'<<pow(Jtot/sigma_rad*c_light/4., 0.25);      //col 5: The radiative temperature 
                
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
                    double tau_inv = 1. / (dx * rhokr) ;
                
                    //double tau_inv = 0.5 / (dx * (total_opacity(i,b) + total_opacity(i+1,b))) ;
                    double R       = xi_rad * tau_inv * std::abs(Jrad_FLD(i+1,b) - Jrad_FLD(i,b)) / (Jrad_FLD(i, b) + 1e-300) ;
                    double flux_limiter = 0;
                    if (R <= 2)
                        flux_limiter = 2 / (3 + std::sqrt(9 + 10*R*R)) ;
                    else 
                        flux_limiter = 10 / (10*R + 9 + std::sqrt(81 + 180*R)) ;
                    
                    double D       = surf[i] * flux_limiter * tau_inv;
                    double flux    = - 4. * pi * D * (Jrad_FLD(i+1,b) - Jrad_FLD(i,b));

                    outfileDiagnostic<<'\t'<<flux;
                }
                
                //Convective fluxes
                if(use_convective_fluxes) { 
                    
                    for (int s=0; s < num_species; s++) {
                        outfileDiagnostic<<'\t'<<species[s].lconvect.at(i);
                    }
                        
                        
                }
                
                outfileDiagnostic<<endl;
            }
            
            
            
        }else cout << "Unable to open file" << filenameDiagnostic << endl; 
        outfileDiagnostic.close();
    
    
}

/**
 * Print monitor file. Used to monitor, mass momentum and energy conservation in a cartesian box without any forces acting. Currently disabled output.
 * 
 * @param[in] num_steps Current step number in the main loop.
 */
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

/**
 * Writing the execution directory and the parameter file with a timestap into the execution_log.txt
 * 
 * @param[in] dir directory string
 * @param[in] par parameter string
 */
void c_Sim::write_into_execution_log(string dir, string par, string spcfile) {
    
    string exlog = "execution_log.txt";
    ofstream exlogstream(exlog, std::ios_base::app);
    
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H-%M-%S");
    //oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    auto str = oss.str();
    
    exlogstream<<str<<": "<<dir<<par;
    
    if(spcfile.compare("default.spc")!=0)
        exlogstream<<"   with "<<spcfile;

    if(intent.compare("---")!=0)
       exlogstream<<" intent: "<<intent<<endl;
    else
       exlogstream<<endl;
    
    exlogstream.close();
}

/**
 * Reads string lines from files and interprets them as photoreactions if they start with '$' or as thermoreactions if they start with '%'
 * 
 * @param[in] dir directory string
 * @param[in] filename the file containing the list of photoreactions and thermoreactions. Can be the species file.
 */
void c_Sim::interpret_chem_reaction_list(string dir, string filename) {
    
    ifstream file(filename);
    string line;
    if(!file) {
        cout<<"Couldnt open reaction file "<<filename<<"!"<<endl;
    }
    
    int found = 0;
    double temp_static_charge = 0.;
    double ns = num_species;
    
    if(debug > 0) cout<<"          In read reactions Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        std::vector<string> stringlist = stringsplit(line," ");

        if(stringlist[0].find("$") != string::npos) { //Photoreactions
            
            if(debug>=2)
                cout<<"Starting to interpret a photoreaction with string = "<<line<<endl;
            
            //int starting_band;
            std::vector<int> e_stoch_i; //Reactant/educt species indices
            std::vector<int> p_stoch_i; //Product species indices
            std::vector<double> p_stoch;
            double branching;
            double energy_threshold;
            int band_number;
            
            std::vector<string> stringlist2 = stringsplit(line,"|");
            std::vector<string> reactionstring = stringsplit(stringlist2[0],"->");
            std::vector<string> eductstring    = stringsplit(reactionstring[0]," +");
            std::vector<string> productstring  = stringsplit(reactionstring[1]," +");
            
            eductstring[0].erase(eductstring[0].begin()); //Remove $ from reactant list
            //cout<<" searching in "<<eductstring[0]<<" for "<<stringsplit(eductstring[0]," ")[1]<<endl;
            
            e_stoch_i.push_back(get_species_index( stringsplit(eductstring[0]," ")[1], 1) );
            energy_threshold = std::stod(stringsplit(eductstring[1]," ")[0] ); 
            //cout<<"energy_threshold = "<<energy_threshold<<endl;
            
            for(string elm : productstring) {
                //cout<<" checking product elm = "<<elm<<endl;
                
                std::vector<string> element = stringsplit(elm," ");
                double stoch = std::stod(element[0]);
                int    prod  = get_species_index( element[1], 1);
                
                //cout<<"    adding for stoch/prod "<<stoch<<"/"<<prod<<endl;
                
                p_stoch.push_back(stoch);
                p_stoch_i.push_back(prod);
            }
            
            //for(string elm : stringlist2)
            //    cout<<" "<<elm<<endl;
            //cout<<endl;
            
            branching = std::stod(stringlist2[1]);
            //cout<<" branching = "<<branching<<endl;
            
            if(stringlist2[2].find("A") != string::npos) {
                //cout<<" starting search for band numbers..."<<endl;
                band_number = find_closest_band(energy_threshold);
            }
            else {
                //cout<<" converted bandnumber = "<<std::stod(stringlist2[2])<<endl;
                band_number = std::stod(stringlist2[2]);
            }
            
            //
            //TODO: Check the data we found before passing it on.
            //
            int checkspassed = 1;
            
            for(int elm : e_stoch_i)
                if(elm==-1) {
                    cout<<"REACTANT NOT FOUND!"<<endl;
                    checkspassed = 0;
                }
            for(int elm : p_stoch_i)
                if(elm==-1) {
                    cout<<"REACTANT NOT FOUND!"<<endl;
                    checkspassed = 0;
                }        
            
            /*cout<<" Photo: "<<band_number<<"/"<<energy_threshold<<"/"<<"/"<<branching;
            for(int elm : e_stoch_i)
                cout<<" "<<elm;
            for(int elm : p_stoch_i)
                cout<<" "<<elm;
            cout<<endl; */
            
            if(checkspassed) {
                //cout<<" Passing data to a photoreaction constructor..."<<endl;
                photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, band_number, e_stoch_i, p_stoch_i, {1.}, p_stoch, branching, energy_threshold )); 
            }
                
            else
                cout<<"Erroneous reaction. Reaction ignored. Reaction string ="<<line<<endl;

            //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, num_bands_in-3, {3}, {10,2}, {1.}, {1.,1.}, 1., 24.6 )); 
        }
        
        if(stringlist[0].find("%") != string::npos) { //Thermoreactions
            
            if(debug>2)
                cout<<"Starting to interpret a thermoreaction with string = "<<line<<endl;
            
            std::vector<int> e_stoch_i;
            std::vector<int> p_stoch_i;
            std::vector<double> e_stoch;
            std::vector<double> p_stoch;
            int is_mtype = 0;
            int is_reverse;
            double a; 
            double b; 
            double c;
            
            std::vector<string> stringlist2 = stringsplit(line,"|");
            std::vector<string> reactionstring;
            
            if(stringlist2[0].find("->") != string::npos) {
                is_reverse = 0;
                reactionstring = stringsplit(stringlist2[0],"->");
            } 
            else if(stringlist2[0].find("<-") != string::npos) {
                is_reverse = 1;
                reactionstring = stringsplit(stringlist2[0],"<-");
            }
            
            if(stringlist2[0].find(" M ") != string::npos)
                is_mtype = 1;
            
            std::vector<string> eductstring    = stringsplit(reactionstring[0]," +");
            std::vector<string> productstring  = stringsplit(reactionstring[1]," +");
            
            eductstring[0].erase(eductstring[0].begin()); //Remove % from reactant list
            
            //cout<<" searching in "<<eductstring[0]<<" for "<<stringsplit(eductstring[0]," ")[0]<<endl;
            
            for(string elm : eductstring) {
                //cout<<" checking educt elm = "<<elm<<endl;
                
                if(elm.find(" M ") == string::npos) { //If its not an M, then we add it to the regular list
                    
                    //cout<<"in elm .. "<<endl;
                    std::vector<string> element = stringsplit(elm," ");
                    double stoch = std::stod(element[0]);
                    int    prod  =  get_species_index( element[1], 1);
                    
                    //cout<<"    adding for stoch/prod "<<stoch<<"/"<<prod<<endl;
                    
                    e_stoch.push_back(stoch);
                    e_stoch_i.push_back(prod);
                }
            }
            for(string elm : productstring) {
                //cout<<" checking product elm = "<<elm<<endl;
                
                if(elm.find(" M ") == string::npos) { //If its not an M, then we add it to the regular list
                    
                    //cout<<"in elm .. "<<endl;
                    std::vector<string> element = stringsplit(elm," ");
                    double stoch = std::stod(element[0]);
                    int    prod  =  get_species_index( element[1], 1);
                    
                    //cout<<"    adding stoch/prod "<<stoch<<"/"<<prod<<endl;
                    
                    p_stoch.push_back(stoch);
                    p_stoch_i.push_back(prod);
                }
            }
            
            a = std::stod(stringsplit(stringlist2[1]," ")[0]);
            b = std::stod(stringsplit(stringlist2[1]," ")[1]);
            c = std::stod(stringsplit(stringlist2[1]," ")[2]);
            
            //
            // TODO: More checks
            //
            int checkspassed = 1;
            
            for(int elm : e_stoch_i)
                if(elm==-1) {
                    cout<<"REACTANT NOT FOUND!"<<endl;
                    checkspassed = 0;
                }
            for(int elm : p_stoch_i)
                if(elm==-1) {
                    cout<<"REACTANT NOT FOUND!"<<endl;
                    checkspassed = 0;
                }   
            
            /*
            cout<<" Thermo: "<<is_reverse<<"/"<<is_mtype<<"/"<<"/"<<a<<"/"<<b<<"/"<<c;
            for(int elm : e_stoch_i)
                cout<<" "<<elm;
            for(int elm : p_stoch_i)
                cout<<" "<<elm;
            cout<<endl;*/
            
            if(checkspassed) {
                //cout<<" Passing data to a thermoreaction constructor..."<<endl;
                reactions.push_back(c_reaction(is_reverse, is_mtype, ns, e_stoch_i, p_stoch_i, e_stoch, p_stoch, a, b, c ));
            }
            else
                cout<<"Erroneous reaction. Reaction ignored. Reaction string = "<<line<<endl;
        }
    
    }
    
    file.close();
    
    
}

/**
 * Reads string lines from files and looks for mentions of resonant pairs. Successful finds are put into the resonant_pair_matrix as constant multiplier
 * 
 * @param[in] dir directory string
 * @param[in] filename the file containing the list of photoreactions and thermoreactions. Can be the species file.
 */
void c_Sim::find_resonant_pairs(string dir, string filename) {

    ifstream file(filename);
    string line;
    if(!file) {
        cout<<"Couldnt open reaction file "<<filename<<"!"<<endl;
    }
    
    int found = 0;
    double temp_static_charge = 0.;
    double ns = num_species;
    
    if(debug > 0) cout<<"          In read reactions Pos1"<<endl;
    
    while(std::getline( file, line )) {
    
        std::vector<string> stringlist = stringsplit(line," ");

        if(stringlist[0].find("&") != string::npos) {
            
            int    a  =  get_species_index( stringlist[1], 0);
            int    b  =  get_species_index( stringlist[2], 0);
            
            if(a != -1 && b != -1) {
                cout<<"Found resonant pair with species "<<a<<" and "<<b<<endl;
                resonant_pair_matrix(a,b) = 10.;
            }
                
        }
        
    }
    
    file.close();
}

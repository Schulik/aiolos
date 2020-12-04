///////////////////////////////////////////////////////////
//
//
//  helpers.cpp
//
// This file contains routines doing small jobs.
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
    
    //Set CFLfactor to safe value once finding the radiative equilibrium is over
    if(globalTime > CFL_break_time)
        cflfactor = 10.;
    
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
std::vector<double> np_ones(int size) {

    return std::vector<double>(size, 1.0) ;
}
//
// Return a 1-D array of a constant value, identical to the numpy function
//
std::vector<double> np_somevalue(int size, double set_value) {

     return std::vector<double>(size, set_value) ;

}

double delta_ij(int i, int j) {
    if(i==j)
        return 1.;
    else 
        return 0.;
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

    
//
//
// Compute Planck integral in a quick way
//
//

double compute_planck_function_integral(double lmin, double lmax, double temperature) {
    
    int num_steps=3;
    double dloggrid = pow(lmax/lmin, 1./((double)num_steps));
    double l1;
    double l2;
    double l_avg;
    double l_avginv;
    double expfactor;
    double lam_db = h_planck*c_light/(kb*temperature)/angstroem;
    double prefactor;
    double tempresult = 0;
    
    l1        = lmin ;//* pow(dloggrid,(double)i);
    l2        = lmax ;//l1 * dloggrid;
    l_avg     = 0.5*(l1+l2);
    l_avginv  = 1./l_avg;
    expfactor = h_planck*c_light/(l_avg*angstroem*kb*temperature);
            
    return sigma_rad2*(l2-l1)/(std::exp(expfactor)-1.)*l_avginv*l_avginv*l_avginv*l_avginv*l_avginv;
}


double compute_planck_function_integral2(double lmin, double lmax, double temperature) {
    
    int num_steps=10;
    double dloggrid = pow(lmax/lmin, 1./((double)num_steps));
    double l1;
    double l2;
    double l_avg;
    double l_avginv;
    double expfactor;
    double lam_db = h_planck*c_light/(kb*temperature)/angstroem;
    double prefactor;
    double tempresult = 0;
    //cout<<"    In compute_planck lmin/lmax/dloggrid = "<<lmin<<"/"<<lmax<<"/"<<dloggrid<<" ";
    
    if(lmin > 3.*lam_db)
        return 2.*c_light*kb*temperature * (1./lmin/lmin/lmin-1./lmax/lmax/lmax)/3.;
    else {
                
        for(int i=0; i<num_steps; i++) {
            //cout<<"---";
            l1        = lmin * pow(dloggrid,(double)i);
            l2        = l1 * dloggrid;
            l_avg     = 0.5*(l1+l2);
            l_avginv  = 1./l_avg;
            expfactor = h_planck*c_light/(l_avg*angstroem*kb*temperature);
            //prefactor = 2*h_planck*c_light*c_light/pow(l_avg,5.)/pow(angstroem,4.);
            
            tempresult += sigma_rad2*(l2-l1)/(std::exp(expfactor)-1.)*l_avginv*l_avginv*l_avginv*l_avginv*l_avginv;
            //tempresult += sigma_rad2*(l2-l1)/pow(l_avg,5.)/(std::exp(expfactor)-1.);
            //tempresult += (l1-l2)/(std::exp(-expfactor)-1.);
            //cout<<" "<<l1<<"/"<<l2<<" "<<prefactor*(l2-l1)/(std::exp(-expfactor)-1.);
        }
        //cout<<endl;
        
        return tempresult;
    }
    
}

double c_Sim::compute_planck_function_integral3(double lmin, double lmax, double temperature) {
    
    double power_min;
    double power_max;
    double lT_min = lmin * temperature;
    double lT_max = lmax * temperature;
    double m;
    int imin;
    int imax;
    
    //
    // Lower power
    //
    if(lT_min < planck_matrix(0,0)) {
        m    = planck_matrix(0,1) / planck_matrix(0,0);
        
        power_min = planck_matrix(0,1) + m * lT_min;
    }
    else {
        imin = std::log(lT_min/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imin+1,1) - planck_matrix(imin,1)) / (planck_matrix(imin+1,0)-planck_matrix(imin,0));
        
        power_min = planck_matrix(imin,1) + m * (lT_min - planck_matrix(imin,0));
    }
    
    //
    // Upper power
    //
    if(lmax * temperature > planck_matrix(num_plancks-1,0)) {
        power_max = 1.;
    }
    else {
        
        imax = std::log(lT_max/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imax+1,1) - planck_matrix(imax,1)) / (planck_matrix(imax+1,0)-planck_matrix(imax,0));
        
        power_max = planck_matrix(imax,1) + m * (lT_max - planck_matrix(imax,0));
    }
    
    if(debug > 1)
        cout<<endl<<"Integral3, lmin/lmax/t = "<<lmin<<"/"<<lmax<<"/"<<temperature<<" imin/imax = "<<imin<<"/"<<imax<<" P(imin)/P(imax) = "<<power_min<<"/"<<power_max<<endl; 
    
    return power_max - power_min;
    
}

   


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Compute analytic solution to the wind problem
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void c_Species::compute_analytic_solution() {
    
    //cout<<"In compute_analytic_solution."<<endl;
    
    for(int i=1;i<=num_cells;i++) {
        
        if(prim[num_cells].sound_speed < 0)
            cout<<"Negative sound speed in compute analytic!"<<endl;
        
        bondi_radius  = base->planet_mass/(2.*prim[num_cells].sound_speed*prim[num_cells].sound_speed);
        double rrc    = base->x_i12[i]/bondi_radius;
        double D      = pow(rrc,-4.) * std::exp(4.*(1.-1./rrc)-1. );
        
        if(base->x_i12[i] < bondi_radius) {
            u_analytic[i] = 0.;//prim[num_cells].sound_speed * std::sqrt( - gsl_sf_lambert_W0(-D) ); 
        }
        else {
            
            u_analytic[i] = 0.;//prim[num_cells].sound_speed * std::sqrt( - gsl_sf_lambert_Wm1(-D) ); 
        }
    
    
    }
    
}


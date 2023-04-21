/**
 * helpers.cpp
 * 
 * This file contains routines doing small jobs.
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// CFL Timestep
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * Computes the cfl timestep based on sound speed, velocity, and a custom internal energy-change criterion.
 * 
 * @return The stable cfl timestep dt in s
 */
double c_Sim::get_cfl_timestep() {
    
    //
    // Compute heuristic radiative timestep
    //
    double maxde = 0;
    
    for(int s = 0; s < num_species; s++) {
        for(int i=num_cells-1; i>0; i--)  {
            species[s].de_e[i] = std::abs(species[s].primlast[i].internal_energy - species[s].prim[i].internal_energy)/species[s].prim[i].internal_energy;
            
            species[s].timesteps_de[i] = dt / species[s].de_e[i] * energy_epsilon;
            
            maxde = std::max(species[s].de_e[i], maxde) ;
            if(debug >= 1 && globalTime > 1e-1)
                cout<<" steps "<<steps<<" species "<<s<<" i = "<<i<<" de/e = "<<species[s].de_e[i]<<" de/e/cflfactor = "<<species[s].de_e[i]/cflfactor<<endl;
        }
    }
    
    if(debug >= 1 && globalTime > 1e-1) {
        char a;
        cin>>a;
    }
    
    timestep_rad2 = dt / maxde * energy_epsilon;

    //
    // Compute individual max wave crossing timesteps per cell
    //  t = delta x / v = delta x / momentum / density
    //
    double minstep = 0.;
    
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
        cflfactor = 0.9;
    
    //Invert and apply CFL secutiry factor
    cfl_step = cflfactor / minstep;
    
    if(do_hydrodynamics)
        return min(cfl_step, dt*max_timestep_change);
    else {
        double ddt = min(timestep_rad2, dt*max_timestep_change);
        return min(ddt, dt_max);
    }
        
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Helper functions
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


std::vector<AOS> init_AOS(int num) {   
    return std::vector<AOS>(num);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////// Helper functions with numpy functionality
///////////////////////////////////////////////////////////////////////

std::vector<double> np_zeros(int size) { return std::vector<double>(size, 0.0) ;}
std::vector<double> np_ones(int size) { return std::vector<double>(size, 1.0) ;}
std::vector<double> np_somevalue(int size, double set_value) { return std::vector<double>(size, set_value) ; }

std::vector<int> inp_zeros(int size) { return std::vector<int>(size, 0.0) ;}
std::vector<int> inp_ones(int size) { return std::vector<int>(size, 1.0) ;}
std::vector<int> inp_somevalue(int size, int set_value) { return std::vector<int>(size, set_value) ; }

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

/**
 * Various ways to compute the Planck integral normalised to 0 and 1. 
 * 
 * Used in initialising the top-of-atmosphere fluxes in init_and_bounds.cpp, as well as in radiation.cpp for the self-radiation.
 * The exact details of when and why four versions were implemented are lost to time, so better not touch this.
 */
double compute_planck_function_integral(double lmin, double lmax, double temperature) {
    
    //int num_steps=3;
    //double dloggrid = pow(lmax/lmin, 1./((double)num_steps));
    double l1;
    double l2;
    double l_avg;
    double l_avginv;
    double expfactor;
    //double lam_db = h_planck*c_light/(kb*temperature)/angstroem;
    //double prefactor;
    //double tempresult = 0;
    
    l1        = lmin ;//* pow(dloggrid,(double)i);
    l2        = lmax ;//l1 * dloggrid;
    l_avg     = 0.5*(l1+l2);
    l_avginv  = 1./l_avg;
    expfactor = h_planck*c_light/(l_avg*angstroem*kb*temperature);
            
    return sigma_rad2*(l2-l1)/(std::exp(expfactor)-1.)*l_avginv*l_avginv*l_avginv*l_avginv*l_avginv;
}


double compute_planck_function_integral2(double lmin, double lmax, double temperature) {
    
    int num_steps=10000;
    double dloggrid = pow(lmax/lmin, 1./((double)num_steps));
    double l1;
    double l2;
    double l_avg;
    double l_avginv;
    double expfactor;
    //double lam_db = h_planck*c_light/(kb*temperature)/angstroem;
    //double prefactor;
    double tempresult = 0;
    //cout<<"    In compute_planck lmin/lmax/dloggrid = "<<lmin<<"/"<<lmax<<"/"<<dloggrid<<" ";
            
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

double c_Sim::compute_planck_function_integral3(double lmin, double lmax, double temperature) {
    
    double power_min;
    double power_max;
    double lT_min;
    double lT_max;

    if (num_bands_out == 1)
        return 1 ;
    
    if(temperature < 2.71) {
        lT_min = lmin * 2.71;
        lT_max = lmax * 2.71;
        
    } else {
        lT_min = lmin * temperature;
        lT_max = lmax * temperature;
    }
    
    double m;
    int imin = 0;
    int imax = num_plancks;
    
    //if();
    //if(debug > 1)
    //int temp_imin = std::log(lT_min/planck_matrix(0,0)) / std::log(lT_spacing);
    //if( steps >= 3927 && temp_imin < 0 )
    //    cout<<"Planck Integral3, lmin/lmax/t = "<<lmin<<"/"<<lmax<<"/"<<temperature<<" lT_min / P00 = "<<lT_min<<" / "<<planck_matrix(0,0)<<" imin = "<<temp_imin<<endl;
    
    //
    // Lower power
    //
    if(lT_min < planck_matrix(0,0)) {
        m    = planck_matrix(0,1) / planck_matrix(0,0);
        
        power_min = planck_matrix(0,1) + m * lT_min;
        
        if(lT_max < planck_matrix(0,0)) //Do this only in the lowermost band
            return 1.;
    }
    else {
        imin = std::log(lT_min/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imin+1,1) - planck_matrix(imin,1)) / (planck_matrix(imin+1,0)-planck_matrix(imin,0));
        
        power_min = planck_matrix(imin,1) + m * (lT_min - planck_matrix(imin,0));
    }
    
    if(debug > 1)
        cout<<" imin/imax = "<<imin;
    
    //
    // Upper power
    //
    if(lmax * temperature > planck_matrix(num_plancks-1,0)) {
        power_max = 1.;
        
        if(lT_min > planck_matrix(num_plancks-1,0)) //Do this only in the uppermost band
            return 1;
    }
    else {
        
        imax = std::log(lT_max/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imax+1,1) - planck_matrix(imax,1)) / (planck_matrix(imax+1,0)-planck_matrix(imax,0));
        
        power_max = planck_matrix(imax,1) + m * (lT_max - planck_matrix(imax,0));
    }
    
    if(debug > 1)
        cout<<" / "<<" P(imin="<<imin<<")/P(imax="<<imax<<") = "<<power_min<<"/"<<power_max<<" = "<<power_max-power_min<<endl; 
    
    return power_max - power_min;
    
}

double c_Sim::compute_planck_function_integral4(double lmin, double lmax, double temperature) {
    
    double power_min;
    double power_max;
    double lT_min;
    double lT_max;

    if (num_bands_in == 1)
        return 1 ;
    
    if(temperature < 2.71) {
        lT_min = lmin * 2.71;
        lT_max = lmax * 2.71;
        
    } else {
        lT_min = lmin * temperature;
        lT_max = lmax * temperature;
    }
    
    double m;
    int imin = 0;
    int imax = num_plancks;
    
    //if();
    //if(debug > 1)
    //int temp_imin = std::log(lT_min/planck_matrix(0,0)) / std::log(lT_spacing);
    //if( steps >= 3927 && temp_imin < 0 )
    //    cout<<"Planck Integral3, lmin/lmax/t = "<<lmin<<"/"<<lmax<<"/"<<temperature<<" lT_min / P00 = "<<lT_min<<" / "<<planck_matrix(0,0)<<" imin = "<<temp_imin<<endl;
    
    //
    // Lower power
    //
    if(lT_min < planck_matrix(0,0)) {
        m    = planck_matrix(0,1) / planck_matrix(0,0);
        
        power_min = planck_matrix(0,1) + m * lT_min;
        
        if(lT_max < planck_matrix(0,0)) //Do this only in the lowermost band
            return 1.;
    }
    else {
        imin = std::log(lT_min/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imin+1,1) - planck_matrix(imin,1)) / (planck_matrix(imin+1,0)-planck_matrix(imin,0));
        
        power_min = planck_matrix(imin,1) + m * (lT_min - planck_matrix(imin,0));
    }
    
    if(debug > 1)
        cout<<" imin/imax = "<<imin;
    
    //
    // Upper power
    //
    if(lmax * temperature > planck_matrix(num_plancks-1,0)) {
        power_max = 1.;
        
        if(lT_min > planck_matrix(num_plancks-1,0)) //Do this only in the uppermost band
            return 1;
    }
    else {
        
        imax = std::log(lT_max/planck_matrix(0,0)) / std::log(lT_spacing);
        m    = (planck_matrix(imax+1,1) - planck_matrix(imax,1)) / (planck_matrix(imax+1,0)-planck_matrix(imax,0));
        
        power_max = planck_matrix(imax,1) + m * (lT_max - planck_matrix(imax,0));
    }
    
    if(debug > 1)
        cout<<" / "<<" P(imin="<<imin<<")/P(imax="<<imax<<") = "<<power_min<<"/"<<power_max<<" = "<<power_max-power_min<<endl; 
    
    return power_max - power_min;
    
}

/**
 * Look for a species name in the list of species and return its index.
 * 
 * @param[in] name species name string, as read in from the *.spc file
 * @return Integer number between 0 and s-1
 */
int c_Sim::get_species_index(const string name) {
    
    for(int s = 0; s<num_species; s++) {
        //cout<<" checking speciesname ="<<species[s].speciesname<<" while looking for "<<name;
        //cout<<" resulting in "<<species[s].speciesname.compare(name)<<endl;
        if(species[s].speciesname.compare(name)==0)
            return s;
    }
    
    return -1;
}


/**
 * Compute analytic solution to the wind problem
 * Update: We do not use this function anymore as too many users have trouble getting the gsl lambert_W function. 
 *         If direct in-code comparison is required, comment gsl_lambertW back in and include the gsl library.
 * 
 */
void c_Species::compute_analytic_solution() {
    
    /*
    for(int i=1;i<=num_cells;i++) {
        
        if(prim[num_cells].sound_speed < 0)
            cout<<"Negative sound speed in compute analytic!"<<endl;
        
        bondi_radius  = base->planet_mass/(2.*prim[num_cells].sound_speed*prim[num_cells].sound_speed);
        double rrc    = base->x_i12[i]/bondi_radius;
        double D      = pow(rrc,-4.) * std::exp(4.*(1.-1./rrc)-1. );
        
        if(base->x_i12[i] < bondi_radius) {
            u_analytic[i] = prim[num_cells].sound_speed * std::sqrt( - gsl_sf_lambert_W0(-D) ); 
        }
        else {
            
            u_analytic[i] = prim[num_cells].sound_speed * std::sqrt( - gsl_sf_lambert_Wm1(-D) ); 
        }
    
    
    }*/
    
}

void c_Species::init_analytic_wind_solution() {
    
    /*
    double sonic_radius = base->init_sonic_radius;
    
    double ufinal = 0;
    
    for(int i=1;i<=num_cells;i++) {
        double sound_speed  = prim[i].sound_speed; // std::sqrt(base->planet_mass/(2. * sonic_radius) );
        double rrc    = base->x_i12[i]/sonic_radius;
        double D      = pow(rrc,-4.) * std::exp(4.*(1.-1./rrc)-1. );
        
        if(base->x_i12[i] < sonic_radius) {
            ufinal = sound_speed * std::sqrt( - gsl_sf_lambert_W0(-D) ); 
        }
        else {
            
            ufinal = sound_speed * std::sqrt( - gsl_sf_lambert_Wm1(-D) ); 
        }
        
        prim[i].speed = ufinal;
    }
 
    eos->compute_conserved(&(prim[0]), &(u[0]), num_cells);
    */
}

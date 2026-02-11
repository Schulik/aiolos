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
    
    int most_restrictive_cell = 0;
    //
    // Compute heuristic radiative timestep
    //
    double maxde = 0;
    double minstep = 0.;
    double diffstep = 1e99;

    int cnstr_spc = -1;
    int cnstr_cell= -1;
    double max_temper = 0;
    double max_mach = 0;
    
    for(int s = 0; s < num_species; s++) {

        for(int i=num_cells-1; i>0; i--)  {
            species[s].de_e[i] = std::abs(species[s].primlast[i].internal_energy - species[s].prim[i].internal_energy)/species[s].prim[i].internal_energy;
            
            species[s].timesteps_de[i] = 1e-50;
            if(species[s].de_e[i] > 0)
                species[s].timesteps_de[i] = dt / species[s].de_e[i];
            
            maxde = std::max(species[s].de_e[i], maxde) ;
        }
    }
    
    timestep_rad2 = dt / maxde * energy_epsilon;

    //
    // Compute individual max wave crossing timesteps per cell
    //  t = delta x / v = delta x / momentum / density
    //

    max_snd_crs_time=0;
    for(int s=0; s < num_species; s++) {
        
        int start = 1;
        if(s == -100) {
           start = ignore_electron_cfl_cell;
           start = std::min(ignore_electron_cfl_cell, num_cells);
        }

        species[s].snd_crs_time = 0;
        for(int i=start; i<=num_cells; i++) {
            
            //Computing the inverse timesteps first
            species[s].timesteps[i]    = std::abs(species[s].prim[i].speed / dx[i]); 
            species[s].timesteps_cs[i] = species[s].prim[i].sound_speed / dx[i];
            
            if(s== e_idx) {
                double f    = 1.;
                double flim = mix_p3;
                if(solver == HydroSolver::mix) {
                    f = get_electron_fraction(i);
                    f = 1.-1./std::exp( f*f/flim/flim );
                    f = std::max(f,1e-10);
                }
                species[s].finalstep[i] = std::sqrt(species[s].timesteps[i]*species[s].timesteps[i] + f * species[s].timesteps_cs[i]*species[s].timesteps_cs[i]);

                if(solver == HydroSolver::implicitelectrons) {
                    species[s].finalstep[i] /= cflfactor_electron;
                    //species[s].finalstep[i] = 1e-20 ; //Naively this should be the right approach, but there are numerical imbalances
                }

                species[s].finalstep[num_cells] = 1e-20;
            }
            else
                species[s].finalstep[i]    = std::sqrt(species[s].timesteps[i]*species[s].timesteps[i] + species[s].timesteps_cs[i]*species[s].timesteps_cs[i] ) ;
            
            species[s].snd_crs_time += 2.* dx[i] / species[s].prim[i].sound_speed ;
            
            if(species[s].finalstep[i] > minstep) {
                minstep = species[s].finalstep[i] ;
                cnstr_spc = s;
                cnstr_cell= i;

                max_temper = species[s].prim[i].temperature;
                max_mach   = std::abs(species[s].prim[i].speed/species[s].prim[i].sound_speed);
            }
            //minstep = std::max(minstep, species[s].finalstep[i]) ;
        }
        max_snd_crs_time = std::max(max_snd_crs_time, species[s].snd_crs_time) ;
    }
    
    //Set CFLfactor to safe value once finding the radiative equilibrium is over
    if(globalTime > CFL_break_time)
        cflfactor = 0.9;
    
    //Invert and apply CFL secutiry factor
    cfl_step = cflfactor / minstep;
    
    if(diffusivity_style > -1) {
        for(int s=0; s < num_species; s++) {
                for(int i=1; i<=num_cells; i++) {
                    //Get diffusive timestep
                    //    cout<<" s / cfl_step / diff_step = "<<s<<" / "<<cfl_step<<" / "<<diffstep<<endl;
                    diffstep = std::min(diffstep, std::min(0.4*species[s].diffusive_timestep(i), 1e99 ));
                }
        }
    }

      
    double final_dt = 0;
    if(do_hydrodynamics) {
        final_dt = min( std::min(cfl_step, diffstep), dt*max_timestep_change);

        if(steps%480==0) {
        cout<<"       max limiting cell: "<<cnstr_cell<<" s= "<<species[cnstr_spc].speciesname<< " => dt= "<<cfl_step<<" dt_diff "<<diffstep<<" dt_energy "<<timestep_rad2<<"  "<<" total dt "<< final_dt<<" max_T = "<<max_temper<<" max_mach "<<max_mach<<" steps= "<<steps<<" t= "<<globalTime<<endl;
        }
        return final_dt;

    } else {
        double ddt = min(std::min(timestep_rad2, diffstep), dt*max_timestep_change);
        final_dt = min(ddt, dt_max);

        if(steps%480==0) {
        cout<<"       max limiting cell: "<<cnstr_cell<<" s= "<<species[cnstr_spc].speciesname<<" dt_diff "<<diffstep<<" dt_energy "<<timestep_rad2<<"  "<<" total dt "<< final_dt<<" max_T = "<<max_temper<<" max_mach "<<max_mach<<" steps= "<<steps<<" t= "<<globalTime<<endl;
        }
    }
    return final_dt;
}


/**
 * Computes the cfl timestep based on sound speed, velocity, and a custom internal energy-change criterion.
 * New part in the new routine: reshuffle the nonsensical double-inversion of timesteps, so that it becomes easier to understand and update
 * with the left and right timesteps.
 * Also make sure the CFL safetyfactor works, because in the old one it mostly didn't.
 * 
 * @return The stable cfl timestep dt in s
 */
double c_Sim::get_cfl_timestep2() {
    
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
    double s_l =0;
    double s_r =0;
    
    max_snd_crs_time=0;
    for(int s=0; s < num_species; s++) {
        
        species[s].snd_crs_time = 0;
        for(int i=1; i<=num_cells; i++) {
            
            //Computing the inverse timesteps first
            //species[s].timesteps[i]    = std::abs(species[s].prim[i].speed / dx[i]); 
            //species[s].timesteps_cs[i] = species[s].prim[i].sound_speed / dx[i];
            s_l = std::sqrt(species[s].prim_l[i].speed*species[s].prim_l[i].speed + species[s].prim_l[i].sound_speed*species[s].prim_l[i].sound_speed);
            s_r = std::sqrt(species[s].prim_r[i].speed*species[s].prim_r[i].speed + species[s].prim_r[i].sound_speed*species[s].prim_r[i].sound_speed);
            //species[s].finalstep[i]    = species[s].timesteps[i] + species[s].timesteps_cs[i] ;
            species[s].finalstep[i]    = std::max(s_l,s_r)/dx[i] ;
            minstep = std::max(minstep, species[s].finalstep[i]) ;
            
            if(i>1 && i<num_cells-1) {
                species[s].finalstep[i]    = std::max(s_l,s_r)/dx[i-1] ; //Take neighbouring cells into account
                minstep = std::max(minstep, species[s].finalstep[i]) ;
                
                species[s].finalstep[i]    = std::max(s_l,s_r)/dx[i+1] ;
                minstep = std::max(minstep, species[s].finalstep[i]) ; 
            }
            
            
            species[s].snd_crs_time += 2.* dx[i] / species[s].prim[i].sound_speed ;
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

/**
 * Loops through the waveband limits and finds the closest low-energy-band limit to energy_threshold
 * 
 * @param[in] energy_threshold the photon energy
 * @return band number 
 */
int c_Sim::find_closest_band(double energy_threshold) {
    
    for(int b=0; b<num_bands_in; b++) {
       double photon_energy = 1.24/( l_i_in[b + 1] ); 
       
       if(photon_energy < energy_threshold)
           return b - 1;
    }
    
    return num_bands_in - 1;    
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
int c_Sim::get_species_index(const string name, const int verbose=0) {
    
    std::vector<string> stringlist = stringsplit(name," ");

    for(auto ss: stringlist) {
        for(int s = 0; s<num_species; s++) {

            //Old debugging block, keep in case things break
            /*for(auto ss: stringlist) {
                cout<<" strnglist element "<<ss<<endl;
            }

            for(int s = 0; s<num_species; s++) {
            cout<<" species index = "<<species[s].this_species_index<<endl;
            cout<<" species charge = "<<species[s].static_charge<<endl;
            cout<<" species fraction = "<<species[s].initial_fraction<<endl;
            cout<<" speciesmass = "<<species[s].mass_amu<<endl;
            //cout<<" speciesname = "<<species[s].speciesname<<endl;
            }*/

            //cout<<" checking speciesname["<<s<<"] = "<<endl;
            //cout<<species[s].speciesname<<endl;
            //cout<<species[s].speciesname<<" while looking for "<<stringlist[i]<<" name "<<name<<endl;
            //cout<<" resulting in "<<species[s].speciesname.compare(name)<<endl;
            //cout<<" resulting in "<<species[s].speciesname.compare(stringlist[i])<<endl;

            if(species[s].speciesname.compare(ss)==0) {
                //if(verbose==1)
                //    cout<<" Found species index for "<<ss<<" = "<<species[s].speciesname<<endl;
                return s;
            }             
        }
        
    }
    
    cout<<" Couldn't find species index for searchlist = "<<name<<endl;
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

double c_Sim::get_electron_fraction(int j) {
	double f =0;
	double tot_press = 0;
	double e_press=0;

	for(int s=0; s<num_species; s++) {
            tot_press += species[s].prim[j].pres;
        }
        if (e_idx > -1)
		e_press = species[e_idx].prim[j].pres;

	return e_press/tot_press;
}

std::string cnstWidth( int value, int width )
 {
      std::ostringstream results;
      results.fill( ' ' );
     results.setf( std::ios_base::internal, std::ios_base::adjustfield );
     results << std::setw( value < 0 ? width + 1 : width ) << value;
      return results.str();
 }
 
 void c_Sim::empty_reaction_table() {
     
     for(int j=0; j<num_cells+2; j++) {
         for(int r=0; r<num_reactions; r++){
             reaction_rate_table(j,r) = 0.;
         }
    }
}

void c_Sim::find_reactionrates_relating_to_species(int cell, string target_speciesname) {
    
    int found = 0;
    int speciesindex = 0;
    
    for(c_reaction& reaction : reactions) {
        int sw = 0;
        //for(int s=0; s<num_species; s++) {
            //double dndt_old;

        for(int& ei : reaction.educts) {
                if( species[ei].speciesname == target_speciesname) {
                    sw = 1;
                    found++;
                    speciesindex = ei;
                }
        }
        for(int& pi : reaction.products) {
                if( species[pi].speciesname == target_speciesname) {
                    sw = 1;
                    found++;
                    speciesindex = pi;
                }
        }
        if(found==1) {
            cout<<" n = "<<species[speciesindex].prim[cell].number_density<<" "<<species[speciesindex].speciesname;
            found++;
        }
        if(sw==1) {
            cout<<" "<<reaction.reaction_number<<" "<<reaction.dndt_old<<" ";
        }
            //reaction.set_reac_number(cnt);
            //reaction.set_base_pointer(this);
            //cnt++;
            //num_reactions++;
    //}
    //for(int rr=0; rr<num_reactions; rr++) {
    }
    if(found > 0)
        //cout<<endl;
        cout<<" stp "<<steps<<endl;
}



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
// Precondition matrix and rhs rows for num_species^2 matrices, e.g. those used in chem and heat exchange routines
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
Vector_t c_Sim::return_preconditioned_LU_solution(const Matrix_t &matrix, const Vector_t &rhs, const Vector_t &orig_vector, Eigen::PartialPivLU<Matrix_t>& LUobject, int cell ) {

    Vector_t result = Vector_t(num_species);
    result.setZero();
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Preconditioning block
    ////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    
    Matrix_t rescl_matrix = matrix;

    Vector_t rescl_x   = Vector_t(num_species);
    Vector_t rescl_r   = Vector_t(num_species);
    Vector_t rescl_rhs = rhs;
    Matrix_t diag_x;
    Matrix_t diag_r;
    Matrix_t Aix ;

    int use_preconditioning = 1;
    if(use_preconditioning) {
        
        for(int s=0; s<num_species; s++) {
            rescl_x(s) = 1/orig_vector(s);
        }

        diag_x = rescl_x.asDiagonal();
        Aix    = rescl_matrix * diag_x;
        
        for(int si=0; si<num_species; si++) {
            double rowmax = 1;
            for(int sj=0; sj<num_species; sj++) {
                rowmax = std::max(rowmax, 1./(std::fabs(Aix(si,sj))+std::fabs(rhs(si)))   );
            }
            rescl_r(si) = rowmax;
        }
        diag_r = rescl_r.asDiagonal();

        rescl_matrix      = diag_r * Aix;
        rescl_rhs       = diag_r * rhs;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    // End Preconditioning
    ////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    //Solve and remap
    LUobject.compute(rescl_matrix) ;
    result.noalias() = LUobject.solve(rescl_rhs);
    result = diag_x * result;


    if( (steps==459) && (cell==10)){
        //cout<<"in precondition matrix steps "<<steps<<endl<<matrix<<endl<<rhs<<endl<<result<<endl;
    }

    return result;

}






//
// This routine constructs a T_mean after the hydro step and is therefore using hydro variables.
// It will be ultimately used to correct negative temperatures which can occur from the computation of E-Ekin at large mach number
//
void c_Sim::update_T_mean(int j, int flag) {

        double avgT_nom   = 0;
        double avgT_denom = 0;
        int sw = 0;

        for(int si=0; si<num_species; si++) {
            double tt = species[si].u[j].u3 - 0.5 * species[si].u[j].u2*species[si].u[j].u2/species[si].u[j].u1;
                   tt /= (species[si].u[j].u1*species[si].cv);
            //cout<<" in update T_mean "<<si<<" "<<tt<<endl;

            if( (tt>0) && (!std::isnan(tt))  ) { //Ignore broken contributions
                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
                avgT_denom += species[si].u[j].u1 * species[si].cv;
            }
            if( tt<0  ) { 
                sw =1;
            }

        }
        T_mean[j] = avgT_nom/avgT_denom;
    
        if(std::isnan(T_mean[j])) {
            cout<<" T_mean is NaN in j= "<<j<<" steps "<<steps<<" flag "<<flag<<" ";
            for(int si=0; si<num_species; si++) {
                cout<<species[si].prim[j].temperature;
                if(std::isnan(species[si].prim[j].temperature))
                    cout<<"("<<species[si].speciesname<<")";
                cout<<" ";
 /*                AOS      tmp  = species[si].u[j];
                AOS_prim tmpp = species[si].prim[j];
                double tt = tmp.u3 - 0.5 * tmp.u2*tmp.u2/tmp.u1;
                tt /= (tmp.u1*species[si].cv);
                cout<<" in update T_mean "<<si<<" "<<tt<<" u = "<<tmp.u1<<" "<<tmp.u2<<" "<<tmp.u3<<" prim = "<<tmpp.internal_energy<<" "<<tmpp.temperature<<" "<<tmpp.pres<<" "<<tmpp.sound_speed<<" "<<tmpp.speed<<endl;
                 */
            }
            cout<<endl;

         /*    char a;
            cin>>a; */
        }
        if(sw==1) {
            cout<<" T_mean contains negatives! j= "<<j<<" steps "<<steps<<" flag "<<flag<<" ";
            for(int si=0; si<num_species; si++) {
                cout<<species[si].prim[j].temperature;
                if(species[si].prim[j].temperature<0)
                    cout<<"("<<species[si].speciesname<<")";
                cout<<" ";
            }
            cout<<endl;


/*             char a;
            cin>>a; */
        }
}

//
// Similar to update_T_mean but no debug functionality
//
double c_Sim::return_T_mean(int j) {

    double avgT_nom   = 0;
        double avgT_denom = 0;

        for(int si=0; si<num_species; si++) {
            double tt = species[si].prim[j].temperature;

            if( (tt>0) && (!std::isnan(tt))  ) { //Ignore broken contributions
                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
                avgT_denom += species[si].u[j].u1 * species[si].cv;
            }

        }
        return avgT_nom/avgT_denom;
}


//
// Similar to update_T_mean but no debug functionality
//
double c_Sim::return_T_mean(int j, Vector_t passed_temps) {

    double avgT_nom   = 0;
    double avgT_denom = 0;
    int cnt=0;

        for(int si=0; si<num_species; si++) {
            double tt = passed_temps(si);

            if( (tt>0) && (!std::isnan(tt))  ) { //Ignore broken contributions
                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
                avgT_denom += species[si].u[j].u1 * species[si].cv;
                cnt++;
            }
        }
        if(cnt==0)
            return -1;
        return avgT_nom/avgT_denom;
}





//
// Compute total internal energy in cell j
//
double c_Sim::return_e_total(int j) {

        double avgT_nom   = 0;

        for(int si=0; si<num_species; si++) {
            double tt = species[si].prim[j].temperature;
            if( (tt>0) && (!std::isnan(tt))  ) { //Ignore broken contributions
                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
            }
        }
        return avgT_nom;
}

//
// Compute total internal energy in cell j
//
double c_Sim::return_e_total(int j, Vector_t passed_temps) {

        double avgT_nom   = 0;

        for(int si=0; si<num_species; si++) {
            double tt = passed_temps(si);
            if( (tt>0) && (!std::isnan(tt))  ) { //Ignore broken contributions
                avgT_nom   += species[si].u[j].u1 * species[si].cv * tt;
            }
        }
        return avgT_nom;
}


//
// Same function, but from a different age and forgotten
//
double c_Sim::return_total_e(int j) {
    double tmp = 0;

    for(int si=0; si<num_species; si++) {
        //tmp += species[si].u[j].u1 * species[si].prim[j].internal_energy;
        tmp += species[si].u[j].u1 * species[si].cv * species[si].prim[j].temperature;
    }

    return tmp;
}



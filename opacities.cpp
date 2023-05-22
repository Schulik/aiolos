///////////////////////////////////////////////////////////
//
//
//  opacities.cpp
//
// This file contains routines computing the attenuation of solar radiation, the resulting heating function S,
// and routines for initializing and computing non-constant opacities.
//
//
//
///////////////////////////////////////////////////////////

#define EIGEN_RUNTIME_NO_MALLOC

#include <cassert>
#include "aiolos.h"

#define c_max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define c_min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


/**
 * Linear interpolation
 */
double lint(double xa, int N, double* X, double* RI) {
    int i2=0, i1=0, i;
    if(xa < X[0])
        xa = X[0];
    else if (xa > X[N])
        xa = X[N];

    // Find position between the two 
    for(i = 0; i < N; i++) {
        if(xa >= X[i]) {
            i1 = i;
            i2 = i+1;
        }
    }
    //Left and right boundaries are known, interpolate linearly
    return (RI[i2] - RI[i1]) / (X[i2] - X[i1]) * (xa - X[i1]) + RI[i1];
}

/**
 * Loglinear interpolation
 */
double logint(double xa, int N, double* X, double* RI) {
    int i2=0, i1=0, i;
    if(xa < X[0])
        xa = X[0];
    else if (xa > X[N])
        xa = X[N];

    // Find position between the two 
    for(i = 0; i < N; i++) {
        if(xa >= X[i]) {
            i1 = i;
            i2 = i+1;
        }
    }
    //Left and right boundaries are known, interpolate logly
    return std::pow(10, (std::log10(RI[i2]) - std::log10(RI[i1])) / (std::log10(X[i2]) - std::log10(X[i1])) * (std::log10(xa) - std::log10(X[i1])) + std::log10(RI[i1]) );
}

/**
 * Wrapper for individual species opacity updates
 */
void c_Sim::update_opacities() {
 
    if(debug > 1)
        cout<<" Starting update_opacities()..."<<endl;
    //
    // Compute kappa based on a species individual density, temperature, irradiation temperature.
    // Initial values are valid no matter whether we use the photochemistry module or not.
    //
    for(int s=0; s<num_species; s++)
        species[s].update_opacities(); //now to be found in opacities.cpp
}

/**
 * The central source for all opacity information, for all in/out and high-energy bands
 * Values for opacity_avg_solar, opacity_avg_planck, opacity_avg_rosseland are read from file and initialized in io.cpp, c_Species::read_species_data(..    )
 * Values for optical depths are computed in update_dS(), after photochemistry, as this potentially changes densities
 * 
 * Opacity models currently implemented, encoded via their corresponding chars:
 * U: User defined opacity. Go into user_opacity() and specify your custom algorithm.
 * C: Constant opacities. Modifiers for solar, planck and rosseland opas exist individually.
 * P: 'Physical': Constant with simple pressure-broadening parameterization powerlaw above 0.1 bars
 * F: Freedman model. Only the fit for Rosseland opacities is currently included.
 * M: Malygin2014 (Planck Rosseland Gas opa)/Semenov2003 (Planck Rosseland Dust opa) combined opacities. Solar opas taken from *opa files.
 * D: Dust only model, given a certain dust size.
 * T: Tabulated opacities. Takes p-T dependent solar, planck and rosseland data from *aiopa files.
 * K: 'Kombined' opacities. Planck and Rosseland is p-T dependent from *aiopa data and solar are p-T-constant, but spectrally resolved from *opa files.
 */
void c_Species::update_opacities() {

    if (base->opacity_model == 'U') {
        // User-defined opacities
        user_opacity() ;
    } 
    else if(base->opacity_model == 'P') {
        
        if(is_dust_like == 0) {
            
            //is_gas_like: Either constant or read-in opacities, with additional pressure broadening as suggested by Robinson & Catling 2013 Nat.Geo
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b) = base->const_opacity_solar_factor * opacity_avg_solar(b);// * (1. + pressure_broadening_factor * pow(prim[j].pres/1e5, pressure_broadening_exponent)); 
                    //cout<<base->l_i_in[b]<<"   "<<opacity_avg_solar(b)<<endl;
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)         = base->const_opacity_rosseland_factor * opacity_avg_rosseland(b);// * (1. + pressure_broadening_factor * pow(prim[j].pres/1e5, pressure_broadening_exponent)); 
                    opacity_planck(j,b)  = base->const_opacity_planck_factor * opacity_avg_planck(b);// * (1. + pressure_broadening_factor * pow(prim[j].pres/1e5, pressure_broadening_exponent)); 
                }
            }
        } else { //is_dust_like
            
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b)= opacity_avg_solar(b); 
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)        = opacity_avg_rosseland(b);  
                    opacity_planck(j,b) = opacity_avg_planck(b); 
                }
            }
        }
        
        if(debug > 3) {
            cout<<" Physical opacities used, and species["<<this_species_index<<"] is_dust_like ="<<is_dust_like<<endl;
        }
                
    }
    else if (base->opacity_model == 'F') {
        
        if(is_dust_like == 0) {
            
            //is_gas_like: Either constant or read-in opacities, with additional pressure broadening as suggested by Robinson & Catling 2013 Nat.Geo
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b)= base->const_opacity_solar_factor * base->freedman_opacity(prim[j].pres/1e6, prim[j].temperature, 0.);
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)        = base->freedman_opacity(prim[j].pres/1e6, prim[j].temperature, 0.); 
                    opacity_planck(j,b) = base->freedman_opacity(prim[j].pres/1e6, prim[j].temperature, 0.);
                }
                
            }
        }
        else { //is_dust_like
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b)= opacity_avg_solar(b); 
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)        = opacity_avg_rosseland(b);  
                    opacity_planck(j,b) = opacity_avg_planck(b);
                }
                
            }
        }
        
        if(debug > 3) {
            cout<<" Freedman opacities used, and species["<<this_species_index<<"] is_dust_like ="<<is_dust_like<<endl;
        }
        
    }
    else if (base->opacity_model == 'M') {
        
        for(int j=0; j< num_cells+2; j++) {
            
            for(int b=0; b<num_bands_out; b++) {
                double oldp;
                if(base->steps<2)
                    oldp = base->const_opacity_planck_factor * base->opacity_semenov_malygin(0, prim[j].temperature,    prim[j].density, prim[j].pres, this->is_dust_like);
                else
                    oldp = opacity_planck(j,b);
                double newp = base->const_opacity_planck_factor * base->opacity_semenov_malygin(0, prim[j].temperature,    prim[j].density, prim[j].pres, this->is_dust_like);
                opacity_planck(j,b) = std::sqrt(oldp*newp);
                
                double oldopa;
                if(base->steps<2)
                    oldopa = base->const_opacity_rosseland_factor * base->opacity_semenov_malygin(1, prim[j].temperature, prim[j].density, prim[j].pres, this->is_dust_like); 
                else
                    oldopa = opacity(j,b);
                double newopa = base->const_opacity_rosseland_factor * base->opacity_semenov_malygin(1, prim[j].temperature, prim[j].density, prim[j].pres, this->is_dust_like); 
                //if(j>3)
                    //newopa = base->const_opacity_rosseland_factor * base->opacity_semenov_malygin(1, prim[j-1].temperature, prim[j].density, prim[j].pres, this->is_dust_like); 
                opacity(j,b)  = std::sqrt(oldopa*newopa);    
            }
            
            for(int b=0; b<num_bands_in; b++) {
                if(b == num_bands_in-1)
                    opacity_twotemp(j,b) = base->const_opacity_solar_factor * opacity_avg_solar(b); 
                else
                    opacity_twotemp(j,b) = base->const_opacity_solar_factor * opacity_avg_solar(b); 
            }
        }
    }
    else if  (base->opacity_model == 'D') {

        // Toy model due to single grain size.
        if (is_dust_like) {
            const double RHO = 1, BETA=1;
            double size = std::pow(3/(4*M_PI)*mass_amu*amu/RHO,1/3.) ;
            
            double k0 = 3/(4*RHO*size) ;
            double T0 = 1/(2*M_PI * 3.475 * size) ;
            
            double ks ;
            if (BETA == 1) 
                ks = k0 * base->T_star / (T0 + base->T_star) ;
            else
                ks = k0 / (1 + std::pow(base->T_star/T0, -BETA)) ;

            for(int j=0; j< num_cells+2; j++) {
                double k ;
                double T = prim[j].temperature ;

                if (BETA == 1) 
                    k = k0 * T / (T0 + T) ;
                else
                    k = k0 / (1 + std::pow(T/T0, -BETA)) ;

                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b) = ks ;
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)         = k ;
                    opacity_planck(j,b)  = k ;
                }
            }
            
        } else {
            // Choose some default low value for a gas species
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b) = 1e-4 ;
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity(j,b)         = 1e-4 ;
                    opacity_planck(j,b)  = 1e-4 ;
                }
                
            }   
        }
        
        if(debug > 3) {
            cout<<" Simple dust model opacities used, and species["<<this_species_index<<"] is_dust_like ="<<is_dust_like<<endl;
        }
    }
    else if (base->opacity_model == 'T'){ //Tabulated opacities, assuming data comes in cm^2/particle, hence divide by particle mass
        
            for(int j=0; j< num_cells+2; j++) {
                
                for(int b=0; b<num_bands_in; b++) {
                    opacity_twotemp(j,b) = base->const_opacity_solar_factor * interpol_tabulated_opacity( opa_grid_solar_log , b, prim[j].temperature, prim[j].pres) * inv_mass;
                    
                }
                for(int b=0; b<num_bands_out; b++) {
                    opacity_planck(j,b)  = base->const_opacity_planck_factor * interpol_tabulated_opacity( opa_grid_planck_log , b, prim[j].temperature, prim[j].pres) * inv_mass;
                    opacity(j,b)         = base->const_opacity_rosseland_factor * interpol_tabulated_opacity( opa_grid_rosseland_log, b, prim[j].temperature, prim[j].pres) * inv_mass;
                    
                }
                
            }
    } 
    else if (base->opacity_model == 'K'){ //Mixture model: Solar opacities are wavelength-dependent but based on a const-P-T-spectrum whereas the planck and rosseland opas are P-T dependent
    
        for(int j=0; j< num_cells+2; j++) {
            
            for(int b=0; b<num_bands_in; b++) {
                        opacity_twotemp(j,b) = base->const_opacity_solar_factor * opacity_avg_solar(b); // * (1. + pressure_broadening_factor * pow(prim[j].pres/1e5, pressure_broadening_exponent)); 
                        
                        if(base->steps >901e99) {
                                cout<<" b_in = "<<b<<" kS = "<<opacity_twotemp(j,b)<<endl;
                                //cout<<base->l_i_in[b]<<"   "<<opacity_avg_solar(b)<<endl;
                        }
                        
                    }
            
            for(int b=0; b<num_bands_out; b++) {
                        opacity_planck(j,b)  = base->const_opacity_planck_factor * interpol_tabulated_opacity( opa_grid_planck_log , b, prim[j].temperature, prim[j].pres) * inv_mass;
                        opacity(j,b)         = base->const_opacity_rosseland_factor * interpol_tabulated_opacity( opa_grid_rosseland_log, b, prim[j].temperature, prim[j].pres) * inv_mass;
                        
                        if(base->steps > 901e99) {
                            cout<<" b_in = "<<b<<" kP = "<<opacity_planck(j,b)<<endl;
                            cout<<" b_in = "<<b<<" kR = "<<opacity(j,b)<<endl;    
                        }
                        
                    }
                    
                    if(base->steps>901e99) {
                        char a;
                        cin>>a;
                    }
                
        }
    
    } else { //opacity_model == 'C', the default
        
        for(int j=0; j< num_cells+2; j++) {
            
            for(int b=0; b<num_bands_in; b++) {

                opacity_twotemp(j,b)= base->const_opacity_solar_factor * opacity_avg_solar(b); //* (1. + pressure_broadening_factor * pow(prim[j].pres/1e5, pressure_broadening_exponent)); 
                
            }
            for(int b=0; b<num_bands_out; b++) {
                opacity(j,b)        = base->const_opacity_rosseland_factor * const_opacity;
                
                if(prim[j].temperature < 5.e2) {
                    opacity_planck(j,b) = base->const_opacity_planck_factor * const_opacity;
                } else {
                    opacity_planck(j,b) = base->const_opacity_planck_factor * const_opacity;// * pow((prim[j].temperature/5.e2),4.);
                }
                
            }
        }
        
        if(debug > 3) {
            cout<<" Constant opacities used, and species["<<this_species_index<<"] is_dust_like ="<<is_dust_like<<endl;
        }
    }
    
    if(debug >= 1 && base->steps <= 2) {
        
        for(int j=0; j< num_cells+2; j++) {
            
            if(j==num_cells+1) {
                
                for(int b=0; b<num_bands_in; b++) {
                    
                    cout<<"steps "<<base->steps<<" j/b/s = "<<j<<"/"<<b<<"/"<<this_species_index<<" opas = "<<opacity_twotemp(j,b)<<" solar_fac = "<<base->const_opacity_solar_factor<<" temper = "<< prim[j].temperature<<" p ="<<prim[j].pres<<endl;
                    
                }
                
                for(int b=0; b<num_bands_out; b++) {
                    
                    cout<<"steps "<<base->steps<<" j/b/s = "<<j<<"/"<<b<<"/"<<this_species_index<<" opap = "<<opacity_planck(j,b)<<" opar = "<<opacity(j,b)<<endl;
                    
                }
                
                cout<<"Stopping to inspect last opas in j =="<<j<<endl;
                //char a;
                //cin>>a;
            }
            
            
        }
    }

}


/**
 * P-T-z dependent opacity fit from Freedman2014.
 * 
 * @param[in] P pressure in dyne/cm^2 i.e. 1e6 dyne/cm^2 = 1 bar
 * @param[in] T temperature in K
 * @param[in] Z metallicity in log(Fe/Fe_solar), i.e. Z=0 is solar
 * @return Rosseland gas opacity (i.e. dust-free) in cm^2/g
 */
double c_Sim::freedman_opacity(double P, double T, double _Z) 
{
        /* Setup the coefficients */
        const double c[8]  = 
            {0, 10.602, 2.882, 6.09e-15, 2.954, -2.526, 0.843, -5.490} ;
        /* Extra digits have been added to ensure continuity */
        double d[6] ;
        if (T < 800) {
            d[0] = -14.051 ;
            d[1] =  3.055  ;
            d[2] =  0.024  ;
            d[3] =  1.877  ;
            d[4] = -0.4451766673621382 ;
            d[5] =  0.8321 ;
        }
        else {
            d[0] =  82.241 ;
            d[1] = -55.456 ;
            d[2] =  8.753389704734504 ;
            d[3] =  0.7048 ;
            d[4] = -0.0414 ;  
            d[5] =  0.8321 ;
        }
        double lgP = std::log10(P) ;
        if(lgP < -2.5) //Safety cap to avoid pole at low pressures
            lgP = -2.5;
        double lgT = std::log10(T) ;
        double lK0 = c[1] * std::atan(lgT - c[2])
            - (c[3] / (lgP + c[4])) * std::exp((lgT - c[5])*(lgT - c[5]))
            + c[6] * _Z + c[7] ;
        double lK1 = d[0] + d[1]*lgT + d[2]*lgT*lgT
            + (d[3] + d[4]*lgT)*lgP
            + d[5] * _Z * (0.5 + std::atan(5*lgT - 12.5)/pi) ;
        
        return std::pow(10, lK0) + std::pow(10, lK1);
}


/**
 * Initializes the combined Semenov03/Malygin14 opacities from my phd
 * 
 * allocates the small double 5x6 arrays for Semenov Dust
 * the double 126x100 for the Malygin Gas opacities
 * And sets the pointers to the allocated memory on the local RAM
 */
void c_Sim::init_malygin_opacities()
{
    int c, err=0;
    char dataname[256];
    FILE *datafile;
    
    sprintf(dataname, "%s%s", "inputdata/", "data_gasopacity_malygin.txt");
    fflush (stdout);
    fprintf(stdout,"The chosen datafile is %s. Also we're in the new function now.\n",dataname);
    fflush (stdout);
    
    datafile = fopen (dataname, "r");
    if (datafile == NULL) {
        cout<<"OPACITYLAW == 2 file data_gasopacity_malygin.txt not found please provide the good path for the file."<<endl;
    }
    fprintf(stdout,"After open.");
    fflush(stdout);
  
    if (datafile == NULL) {
        fprintf (stderr, "WARNING ! Opacity file not found. Stopping. Searched for %s\n",dataname); 
        exit (1);
    }
  
  	err += fscanf(datafile, "%i", &opacity_gas_rows);
	err += fscanf(datafile, "%i", &opacity_gas_cols);
	
	fprintf(stdout,"Initializing opacitydata now with %i rows and %i cols.\n", opacity_gas_rows, opacity_gas_cols);
  	fflush(stdout);

	//opacity_rows & cols: global variables. TO DO if necessary: read their length in dynamically
	//allocate memory for *opa_tscale, *opa_pscale, *opa_ross, *opa_planck;
	opa_gas_pscale = (double *)malloc(sizeof(double) * opacity_gas_rows); //Pressure grid
	opa_gas_tscale = (double *)malloc(sizeof(double) * opacity_gas_cols); //Temperature grid
	
	opa_gas_ross   = (double *)malloc(sizeof(double) * opacity_gas_rows * opacity_gas_cols); //Rosseland opacities
	opa_gas_planck = (double *)malloc(sizeof(double) * opacity_gas_rows * opacity_gas_cols); //Planck opacities
	
	if( (opa_gas_pscale == NULL) || (opa_gas_tscale == NULL) || (opa_gas_ross == NULL) || (opa_gas_planck == NULL)) {
	 	fprintf (stderr, "WARNING ! Opacityfiles failed to initialize. \n"); 
	}
	
	for(c = 0; c < opacity_gas_rows; c++) {
	  	err += fscanf(datafile, "%lf", &opa_gas_pscale[c]);
	}
	for(c = 0; c < opacity_gas_cols; c++) {
	  	err += fscanf(datafile, "%lf", &opa_gas_tscale[c]);
	}
  	for (c = 0; c < opacity_gas_rows * opacity_gas_cols; c++) {
	  err += fscanf(datafile, "%lf", &opa_gas_ross[c]);
  	}
  	for (c = 0; c < opacity_gas_rows * opacity_gas_cols; c++) {
	  err += fscanf(datafile, "%lf", &opa_gas_planck[c]);
  	}
  
  	fclose (datafile);
  
    if(err > 0)
        cout<<"Had "<<err<<" reading errors in init_malygin_opacities. Check what is wrong."<<endl;
}


/**
 * Gas opacities from Malygin+2014 tables
 * This function uses the data for *opa_gas_tscale, *opa_gas_pscale, *opa_gas_ross, *opa_gas_planck and interpolates them bilinearly on the given grid
 * 
 * @param[in] rosseland 0 or 1, if 0 returns planck mean instead of rosseland mean
 * @param[in] T_gas gas temperature in K
 * @param[in] pressure gas pressure in dyne/cm^2
 * @return Gas rosseland or planck mean in cm^2/g
 */
double c_Sim::get_gas_malygin(int rosseland, double T_gas, double pressure) {
	
  	double tP, tT, denom, tempKappa;	//tempPressure and denominator
  	static double mul1, mul2, mul3, mul4; //temp doubles for the interpolation. static for faster speed
  	static int Plow=0, Phigh=0;
	int Tlow=0, Thigh=0;
	int i;
	double *array;
	
	if(rosseland) 
		  array = &opa_gas_ross[0];
	else
		  array = &opa_gas_planck[0];
    
	tP = pressure;
	tT = T_gas;
	/*
	fprintf(stdout,"Input Dens = %e  temperture = %e    pressure = %e\n",rho, tT, tP);
	fflush(stdout);
	*/
	//Search for value in T grid
	if(tT < opa_gas_tscale[0])
// 	  		tT = opa_gas_tscale[0];
 	  		return 1e-20;				//T too low for gas to play a role
	else if (tT > opa_gas_tscale[opacity_gas_cols-1])
		  tT = opa_gas_tscale[opacity_gas_cols-1];

	if(tP < opa_gas_pscale[0])
	  		tP = opa_gas_pscale[0];
	else if (tP > opa_gas_pscale[opacity_gas_rows-2])
		  tP = opa_gas_pscale[opacity_gas_rows-2];

	for(i = 0; i < opacity_gas_cols; i++) {
			if(tT >= opa_gas_tscale[i]) {
			 	Tlow  = i;
			  	Thigh = i+1;
			}
		}
		
    for(int j = 0; j < opacity_gas_rows-1; j++) {
			if(tP >= opa_gas_pscale[j]) {
			 	Plow  = j;
			  	Phigh = j+1;
			}
		}
	
// 	Alternative grid search in P:
	//fac   = ((double) opacity_gas_rows) / log10( opa_gas_pscale[opacity_gas_rows-1]/opa_gas_pscale[0] );
	//Plow  = (int)(fac*log10(tP/opa_gas_pscale[0] ));
	//Phigh = Plow + 1;
// 	This method is much faster, but should we use a different P-T grid the 7.0 would be wrong
	
// 	fprintf(stdout,"%e %i %i %e %e\n",fac, Plow, Phigh, tP, rho);
// 	fflush(stdout);
	
	//fprintf(stdout,"P = %e found between %e and %e for %i & %i.\n",tP, opa_gas_pscale[Plow], opa_gas_pscale[Phigh],Plow, Phigh);
	if(std::isnan(opa_gas_tscale[Thigh]) )
        fprintf(stdout,"T = %e found as NaN between %e and %e for %i & %i.\n",T_gas, opa_gas_tscale[Tlow], opa_gas_tscale[Thigh],Tlow, Thigh);
	//getchar();
	
	//Interpolate on pressure-temperature grid
	
	double pmax = std::log(opa_gas_pscale[Phigh]);
	double pmin = std::log(opa_gas_pscale[Plow]);
	double pval = std::log(tP);
	
	denom = 1.0/((opa_gas_tscale[Thigh] - opa_gas_tscale[Tlow]) * (pmax - pmin));
	
    mul1 = (opa_gas_tscale[Thigh] - tT) * (pmax - pval);
	mul2 = (tT - opa_gas_tscale[Tlow])  * (pmax - pval);
	mul3 = (opa_gas_tscale[Thigh] - tT) * (pval - pmin);
	mul4 = (tT - opa_gas_tscale[Tlow])  * (pval - pmin);
	tempKappa = denom * ( mul1 * std::log10(array[Tlow + Plow * opacity_gas_cols]) + 
	                      mul2 * std::log10(array[Thigh + Plow * opacity_gas_cols]) 
						+ mul3 * std::log10(array[Tlow + Phigh * opacity_gas_cols]) + 
                          mul4 * std::log10(array[Thigh + Phigh * opacity_gas_cols]) );
	
    if(T_gas > 1.0e50 ) {
        fprintf(stdout,"T = %e, rounded %e, p = %e found between %e and %e for Tlow= %i Tlow= %i, pmin,pval,pmax = %e %e %e kappa = %e.\n",T_gas, tT, pressure , opa_gas_tscale[Tlow], opa_gas_tscale[Thigh],Tlow, Thigh,opa_gas_pscale[Plow],tP,opa_gas_pscale[Phigh],std::pow(10.,tempKappa));
        fprintf(stdout,"mul1 = %e, mul2= %e, mul3 = %e, mul4= %e, denom = %e \n",mul1, mul2, mul3, mul4, denom );
        fprintf(stdout,"Interpolation kappas: %e %e %e %e \n",array[Tlow + Plow * opacity_gas_cols],array[Thigh + Plow * opacity_gas_cols],array[Tlow + Phigh * opacity_gas_cols],array[Thigh + Phigh * opacity_gas_cols]);
    }

     
    /*
// 	return wanted value
	
   fprintf(stdout,"Found kappa: %e\n",tempKappa);
  */
  	return std::pow(10.,tempKappa);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~Dust Opacities after Semenov2003 and Gas after Malygin 2014~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


/**
 * Dust opacities from the Semenov+2003 mix. Calls routine to add the gas opacites from Malygin+2014 tables.
 * This function uses the data for *opa_gas_tscale, *opa_gas_pscale, *opa_gas_ross, *opa_gas_planck and interpolates them bilinearly on the given grid
 * 
 * @param[in] rosseland 0 or 1, if 0 returns planck mean instead of rosseland mean
 * @param[in] temperature gas temperature in K
 * @param[in] rho gas density (assuming it dominates the gas/dust mixture) in g/cm^3
 * @param[in] pressure gas pressure in dyne/cm^2
 * @param[in] caller_is_dust if 1, gas opacity is ignored
 * @return Combined rosseland or planck mean in cm^2/g
 */
double c_Sim::opacity_semenov_malygin(int rosseland, double temperature, double rho, double pressure, int caller_is_dust) 
{
	//The Malygin gas-opacities live as data on a grid, thus we create a grid
	//The Semenov dust-opacities live as coefficients in a analytic function, thus we need only those
	//We want this data to be initialized and a pointer to it already lying around, cuz we'll call the opacity function often

	//Rosseland mean opacities // const static 
    double opa_ross[6][6] = { 
{0.244153190819473e-03,-0.191588103585833e-03,0.229132294521843e-03,0.454088516995579e-06,-0.389280273285906e-08,0.319006401785413e-11},
{0.292517586013911e+01,-0.801305197566973e-01,0.923001367150898e-03,-0.349287851393541e-05,0.587151269508960e-08,-0.371333580736933e-11},
{-0.125893536782728e+02,0.160197849829283e+00,-0.582601311634824e-03,0.110792133181219e-05,-0.105506373600323e-08,0.404997080931974e-12},
{-0.192550086994197e+01,0.354301116460647e-01,-0.127355043519226e-03,0.233773506650911e-06,-0.207587683535083e-09,0.728005983960845e-13},
{0.165244978116957e+01,-0.260479348963077e-02,0.590717655258634e-05,-0.300202906476749e-08,0.803836688553263e-12,-0.916709776405312e-16},
{0.00,0.00,0.00,0.00,0.00,0.00}};

	//Planck mean opacities
      double opa_planck[6][6] = {
{0.423062838413742e-03,-0.191556936842122e-02,0.130732588474620e-02,-0.108371821805323e-04,0.427373694560880e-07,-0.664969066656727e-10},
{-0.173587657890234e+01,0.302186734201724e-01,0.311604311624702e-03,-0.210704968520099e-05,0.472321713029246e-08,-0.371134628305626e-11},
{-0.638046050114383e+01,0.120954274022502e+00,-0.436299138970822e-03,0.784339065020565e-06,-0.681138400996940e-09,0.233228177925061e-12},
{-0.345042341508906e+01,0.664915248499724e-01,-0.240971082612604e-03,0.417313950448598e-06,-0.349467552126090e-09,0.115933380913977e-12},
{0.667823366719512e+01,-0.166789974601131e-01,0.238329845360675e-04,-0.140583908595144e-07,0.417753047583439e-11,-0.503398974713655e-15}, 
{0.00,0.00,0.00,0.00,0.00,0.00}};

    //Water ice evaporation temperature, then metallic iron, orthopyroxene and olivine; // depending on density array 'ro'
    double tt[8][4] = {
		   {1.090E+02,8.350E+02,9.020E+02,9.290E+02},
			{1.180E+02,9.080E+02,9.800E+02,9.970E+02},
			{1.290E+02,9.940E+02,1.049E+03,1.076E+03},
			{1.430E+02,1.100E+03,1.129E+03,1.168E+03},
			{1.590E+02,1.230E+03,1.222E+03,1.277E+03},
			{1.800E+02,1.395E+03,1.331E+03,1.408E+03},
			{2.070E+02,1.612E+03,1.462E+03,1.570E+03},
			{2.440E+02,1.908E+03,1.621E+03,1.774E+03}};

    //The density values for which the evaporation temperatures are given
    double ro[8] = {1.0E-18, 1.0E-16, 1.0E-14, 1.0E-12, 1.0E-10, 1.0E-08, 1.0E-06, 1.0E-04 }; 

    // Global variable(s):
    double *eD;
    double aKext;
    // Local variable(s):

    double T_ev[4], temp[8], dT[5], tmax1, tmax2, tmin, T1, T2, TD;
    double T[5], aKrL, aKrR, aKg_ext = 1e-10, AA, BB, FF, kappa_gas, kappa_dust;
    int KK;
    int i, j, iter;
    int smooth;
    
    
    //Cap rho due to the table limit and adjust pressure accordingly
    if(rho > 1e-4) {
        pressure = pressure * 1e-4/rho;
        rho      = 1e-4;
    }
    
//          printf("Initial conditions in Semenov %e %e \n", temperature, rho);
		
      // Initialization of the output parameters:
      aKext = 0.0;
      if (temperature <= 0.9) 
			return 0.9; //T must be more that few [K]
    //.............................................................................
    // Constant(s):
    //.............................................................................
      //PI = 4.0*atan(1.0);
      
      //Choose what we want
    if(rosseland) 
        eD = opa_ross[0];
    else
        eD = opa_planck[0];
      
    //.............................................................................
    // Interpolation of a matrix of evaporation temperatures for a given density
    //.............................................................................
      for(i = 0; i < 4 ; i++) {
         for(j = 0; j < 8 ; j++) {
            temp[j] = tt[j][i];
         }
         T_ev[i] = lint(rho, 7, &ro[0], &temp[0]); //Binterpolation
      }
      
    //.............................................................................
    //  Set up the evaporation temperature array 'T(1:5)':
    //.............................................................................    
    T[0] = T_ev[0];  //evaporation temperature of water ice,
    T[1] = 275.0;  //evaporation temperature of volatile organics,
    T[2] = 425.0;  //evaporation temperature of refractory organics,
    T[3] = 680.0;  //evaporation temperature of troilite,
    tmax1 = c_max(T_ev[1], T_ev[2]);
    tmax2 = c_max(T_ev[2], T_ev[3]);
    tmin  = c_min(tmax1, tmax2);
    T[4] = tmin;     //average evaporation temperatures of iron,
    //                        !olivine, and orthopyroxene
    //.............................................................................  
    // Determination of temperature regimes where smoothing is necessary 
    //............................................................................
      dT[0] = 5.0e0;  //an interval of T where ice is evaporating,
      dT[1] = 5.0e0;  //an interval of T where vol. organics is evaporating,
      dT[2] = 15.0e0;  //an interval of T where ref. organics is evaporating,
      dT[3] = 5.0e0;  //an interval of T where troilite is evaporating,
      dT[4] = 100.0e0; //an wide interval of T where iron, pyroxe, and
    //                            !olivine are evaporating,
    //.............................................................................
    // Determination of a temperature regime for a given temperature:
    //.............................................................................
      
    KK = 5;  //default value -> gas-dominated opacity,
    if ( temperature <= (T[0]+dT[0]) )  
        KK = 0;    
    for(iter = 1; iter < 5; iter++)
        if ((temperature > ( T[iter-1]+dT[iter-1] )) && (temperature <= ( T[iter]+dT[iter] ))) 
            KK = iter;    
       
    //.............................................................................
    // The dust-dominated opacity:
    //.............................................................................
    // Switch to smoothing of the opacity if a temperature is near an 
    // evaporation temperature of a dust material:
    //.............................................................................
    smooth = 0;
    for(i = 0; i < 5; i++) 
        if (abs(temperature - T[i]) <= dT[i]) 
            smooth = 1; 
    //-----------------------------------------------------------------------------
    // If 'T_in' inside of (Tev-dT, Tev+dT) then start smoothing:
    //-----------------------------------------------------------------------------
    if(smooth == 1) {
			T1 = T[KK] - dT[KK];
        	T2 = T[KK] + dT[KK];
        	TD = temperature - T[KK];
    //-----------------------------------------------------------------------------
    // Calculation of a mean opacity (extinction):
    //-----------------------------------------------------------------------------
      	aKrL = ((((eD[(KK+1)*6 - 1]*T1+eD[(KK+1)*6 - 2])*T1+eD[(KK+1)*6 - 3])*T1+eD[(KK+1)*6 - 4])*T1+eD[(KK+1)*6 - 5])*T1+eD[(KK+1)*6 - 6];

       	if(KK == 4) {
      		aKrR = aKg_ext; 
        } else 
      		aKrR = ((((eD[(KK+2)*6 - 1]*T2+eD[(KK+2)*6 - 2])*T2+eD[(KK+2)*6 - 3])*T2+eD[(KK+2)*6 - 4])*T2+eD[(KK+2)*6 - 5])*T2+eD[(KK+2)*6 - 6];
      
        AA = 0.5e0*(aKrL-aKrR);
        BB = 0.5e0*(aKrL+aKrR);
        FF = pi/2.0/dT[KK];
       	aKext = BB-AA*sin(FF*TD);

    } else  {
    //.............................................................................
    //  Smoothing is not necessary, direct calculation by a fit polynom of
    //  fifth degree: y=a*x^5+b*x^4+c*x^3+d*x^2+e*x+f
    //.............................................................................
         aKext = ((((eD[(KK+1)*6 - 1]*temperature+eD[(KK+1)*6 - 2])*temperature+eD[(KK+1)*6 - 3])*temperature+eD[(KK+1)*6 - 4])*temperature+eD[(KK+1)*6 - 5])*temperature+eD[(KK+1)*6 - 6];	
    }
	
	kappa_dust = aKext;
    
    if(rosseland == 0) {
        if(temperature > 700.)
            kappa_gas  = get_gas_malygin(rosseland, temperature, pressure);
        else
            kappa_gas  = get_gas_malygin(rosseland, 700., pressure) * pow(temperature/700., 5.);
    } else {
        if(temperature > 700.)
            kappa_gas  = get_gas_malygin(rosseland, temperature, pressure);
        else {
            kappa_gas  = get_gas_malygin(rosseland, 700., pressure) * pow(temperature/700., 5.);
    
     // double F = (log T - log T_L) / (log T_U - log T_L);       //Volume mixing rule
     // double S = (1-cos(F*pi))/2.;
     // double kappa = S log kappa_U(R,T) - (1-S) log kappa_L(R,T)
            
            //
            // For now i have commented the freedman opacities out as solution for low temperatures, because of the discontinuity at the intersection temperature
            // in k_R. This seems to trigger oscillations in E_rad in very optically thin regions, which destabilize and lead to E_rad < 0 at some point
            //
//             double k1 = freedman_opacity(pressure, temperature, 0.); 
//             double k2 = get_gas_malygin(rosseland, rho, temperature, pressure);
//             kappa_gas = c_max(k1, k2); 
        }
            
    }
    
    //
    // I have switched off the dust opacities as their instantaneous evaporation leads to large oscillations in the opacities.
    // Those didn't influence the temperatures very much, but its ugly and might lead to numerical trouble in the outflow
    //
	if(rosseland == 0) { //Planck opacities are additive
        
        if(caller_is_dust)
            return kappa_dust;
        else
            return kappa_gas + dust_to_gas_ratio * kappa_dust; 
        
    }
	else{ //Rosseland opacities are not additive
        if(caller_is_dust)
            return kappa_dust;
        else
            //return c_max(kappa_gas, dust_to_gas_ratio * kappa_dust); 
            //return kappa_gas + dust_to_gas_ratio * kappa_dust;
            return kappa_gas + 1.*kappa_dust;
        //return kappa_gas + dust_to_gas_ratio * kappa_dust;
        
    }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//This function produces a kappa-landscape
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/**
 * This function produces kappa(T) at constant density or pressure, for any opacity function desired.
 * Just needs to be called once, e.g. at the beginning of c_Sim.execute() in advection.cpp
 */
void c_Sim::kappa_landscape()
{
  	int i,j;
	const int arraysize = 500;
  	double recentpress;
  	//double opa_freedman[arraysize];
  	double opa_ross_semenov[arraysize];
    double opa_planck_semenov[arraysize];
    double opa_planck_twotemp[arraysize];
  	double temperature[arraysize];
	//double density[21]={1e-18,1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,3e-12,6e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-1,1e-0};
    double pressure[12]={1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e18};
	
	printf("In kappa landscape.\n");
	fflush (stdout);
	
	for(j = 0; j < 12; j++) {
	   recentpress = pressure[j];
	   printf("Doing pressure %e \n",pressure[j]);
		fflush (stdout);
  
  		for(i = 0; i < arraysize; i++) {
		  	temperature[i] = 1e+1 + pow( 10.0 , 6.0 * ((double)i/(double)arraysize));
			
			
			fflush (stdout);
            
            //double inv_mass = species[0].inv_mass;
            //opa_freedman[i] = freedman_opacity(pressure , temperature[i], 0.);
//             opa_planck_twotemp[i] = species[0].interpol_tabulated_opacity( species[0].opa_grid_solar , 0, temperature[i], pressure[j]) * inv_mass;
//             opa_planck_semenov[i] = species[0].interpol_tabulated_opacity( species[0].opa_grid_planck , 0, temperature[i], pressure[j]) * inv_mass;
//             opa_ross_semenov[i]   = species[0].interpol_tabulated_opacity( species[0].opa_grid_rosseland, 0, temperature[i], pressure[j]) * inv_mass;
// 			
            //double p_actual = temperature[i] * density[j] * kb * inv_mass;
            //double d_actual = pressure[j] / (kb * temperature[i] * inv_mass);
            //printf("Doing temperature %e and density %e resulting in pressure/bars = %e \n",temperature[i],d_actual, pressure[j]/1e6);
            opa_planck_twotemp[i] = 1.;
            opa_planck_semenov[i] = freedman_opacity(pressure[j] , temperature[i], 0.); //opacity_semenov_malygin(0, temperature[i], d_actual, pressure[j], 0);
            opa_ross_semenov[i]   = freedman_opacity(pressure[j] , temperature[i], 0.); //opacity_semenov_malygin(1, temperature[i], d_actual, pressure[j], 0);
            //printf("Semenov planck done \n");
			//fflush (stdout);
  		}
  		
  		printf("Done press %e \n",pressure[j]);
		fflush (stdout);
  
  		FILE *output;
  		char name[256];
  
  		printf ("Writing opacity landscape for press %e...\n", recentpress);
  		fflush (stdout);
  		sprintf (name, "kappa_freedman_landscape2_%e.dat", recentpress);
  		//output = fopenp (name, "a");
        output = fopen (name, "a");
 		for(i = 0; i < arraysize; i++) {
  				//fprintf (output, "%e\t%e\t%e\t%e\t%e\n",
            
  				fprintf (output, "%e\t%e\t%e\t%e\n",\
	   			temperature[i], opa_ross_semenov[i], opa_planck_semenov[i], opa_planck_twotemp[i]);
 		}
  		fclose (output);
  		fflush (stdout);
	}
}


/**
 * This function uses the data for given on any 1-D array (pseudo 2-D) similar to c_Sim::get_gas_malygin() and interpolates them log-bilinearly on the given grid.
 * Uses 1% higher and lower pressure values and averages them into the final product, to get more stability in regions of wildly oscillating opacities.
 */
double c_Species::interpol_tabulated_opacity(const Eigen::VectorXd& array, int band, double T_gas, double pressure) {
	
  	double tT, denom;	//tempPressure and denominator
  	float tP;
  	static double mul1, mul2, mul3, mul4; //temp doubles for the interpolation. static for faster speed
  	static int Plow=0, Phigh=0;
	int Tlow=0, Thigh=0;
	int i;
    
	tP = pressure;
	tT = T_gas;
    
	//Search for value in T grid
	if(tT < opa_tgrid[0])
 	  		tT = opa_tgrid[0];
 	  		//return 1e-20;				//T too low for gas to play a role
	else if (tT > opa_tgrid[opa_tgrid_size-1])
		  tT = opa_tgrid[opa_tgrid_size-1];

	if(tP < opa_pgrid[0])
	  		tP = opa_pgrid[0];
	else if (tP > opa_pgrid[opa_pgrid_size-1])
		  tP = opa_pgrid[opa_pgrid_size-1];

	//for(i = 0; i < opa_tgrid_size; i++) {
    for(i = opa_tgrid_size-2; i >= 0; i--) {
			if(tT >= opa_tgrid[i]) {
			 	Tlow  = i;
			  	Thigh = i+1;
                break;
			}
        }
        
    if(Thigh == opa_tgrid_size) {
        Thigh = opa_tgrid_size-1;
        Tlow  = opa_tgrid_size-2;
        tT    = opa_tgrid(Thigh);
    }
    
    //for(i = 0; i < opa_pgrid_size; i++) {
    for(i = opa_pgrid_size-2; i >= 0; i--) {
        if(tP >= opa_pgrid[i]) {
			 	Plow  = i;
			  	Phigh = i+1;
                break;
			}
        }
    //
    // Detection of grid boundaries
    //
    if(Phigh == opa_pgrid_size) {
        Phigh = opa_pgrid_size-1;
        Plow  = opa_pgrid_size-2;
        tP    = opa_pgrid(Phigh);
    }
    
    //cout<<" in interpol, Plow/Phigh/Tlow/Thigh = "<<opa_pgrid[Plow]<<"/"<<opa_pgrid[Phigh]<<"/"<<opa_tgrid[Tlow]<<"/"<<opa_tgrid[Thigh]<<"  tT,tP = "<<tT<<"/"<<tP<<endl;
	
	//Get corner values of opacity for interpolation
    //Floats for speed!
 
    float k1 = array(Tlow + Plow * opa_tgrid_size + band * opa_tgrid_size * opa_pgrid_size);
    float k2 = array(Thigh + Plow * opa_tgrid_size + band * opa_tgrid_size * opa_pgrid_size);
    float k3 = array(Tlow + Phigh * opa_tgrid_size + band * opa_tgrid_size * opa_pgrid_size);
    float k4 = array(Thigh + Phigh * opa_tgrid_size + band * opa_tgrid_size * opa_pgrid_size);
    
    float pmax = opa_pgrid_log(Phigh);
	float pmin = opa_pgrid_log(Plow);
    
    denom = ((opa_tgrid(Thigh) - opa_tgrid(Tlow)) * (pmax - pmin));
    denom = 1.0/denom;
    
    float ktotal = 0;
    float muls[3]  = {0.99, 1.00, 1.01};
    
    for(int i=0; i<3; i++) {
        float tpp = tP * muls[i];
        float pval = njuffas_log10f(tpp);
        
        mul1 = (opa_tgrid(Thigh) - tT) * (pmax - pval);
        mul2 = (tT - opa_tgrid(Tlow))  * (pmax - pval);
        mul3 = (opa_tgrid(Thigh) - tT) * (pval - pmin);
        mul4 = (tT - opa_tgrid(Tlow))  * (pval - pmin);
        
        ktotal += denom *     (mul1 * k1
                            + mul2 * k2
                            + mul3 * k3 
                            + mul4 * k4 );
    }
    
  	return std::pow(10., ktotal*0.33333333333333333);
}

#include "main.h"



    
        ///////////////////////////////////////////////////////////
        //
        // Chapter on gravitation
        //
        //
        ///////////////////////////////////////////////////////////



    //
    // Initialize the gravity array with a non-self gravitating solution
    //
    // takes: r as array
    //
    void hydro_run::init_grav_pot() {
        
        for(int i = 1; i <= num_cells; i++) {
            enclosed_mass[i] = planet_mass;      //No self-gravity
            phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            
        }
            

    }

    //
    // Compute the gravity array as a self-gravitating solution
    //
    // takes: r as array
    //    
    void hydro_run::update_mass_and_pot() {
        
        if(debug >= 3)
            cout<<"In update mass. "<<endl;
        
        enclosed_mass[0] = planet_mass;
        
        for(int i = 1; i <= num_cells+1; i++) {
            
            if(debug >= 4) {
                if(i==0)
                    cout<<"In update mass, 0th element."<<endl;
                if(i==num_cells-1)
                    cout<<"In update mass, last element."<<endl;
            }
            
            
            enclosed_mass[i] = enclosed_mass[i-1] +  4. * 3.141592 * (pow(x_i[i],3.)-pow(x_i[i-1],3.) )/3. * u[i].u1; //Straightfoward integration of the poisson eqn
            
            
            if (use_self_gravity == 1) 
                phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            else
                phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[0]);
            
        }
        
        if(debug >= 3)
            cout<<"Done in update mass. "<<endl;
    }
    
    //
    // Initialize the gravity array with a non-self gravitating solution
    //
    // takes: r as single value
    //
    // Options: linear gravity and 1/r gravity with gravitational smoothing from Klahr&Kley 2006
    //
    double hydro_run::get_phi_grav(double &r, double &mass) {
        
        //
        // Deepening the potential if needed
        //
        if(globalTime < rs_time) {
                
            rs_at_moment = 0.5 - (0.5 + rs) * (globalTime/rs_time) ;
        }
        else
            rs_at_moment = rs;
        
        //
        // Linear or 1/r gravity
        //
        if(use_linear_gravity)
            return mass*abs(r);
        else
        {
            if(abs(r-planet_position) > rs_at_moment )
                return -mass/abs(r-planet_position);
            else
                return -mass * (pow(r-planet_position,3.)/pow(rs_at_moment,4.)  - 2.* pow(r-planet_position,2.)/pow(rs_at_moment,3.) + 2./rs_at_moment );
            
        }
            //return -mass/( sqrt( pow(r-planet_position,2.) + rs_at_moment*rs_at_moment) );
            
    }

    
    //
    // p_hydrostatic - hydrostatic reconstruction, computes pressure difference between current pressure and hydrostatic state, after Eq. 16 in KM2016
    //
    // takes:
    // returns:
    //
    double hydro_run::get_p_hydrostatic(AOS &u, double &phi_l, double &phi_r,const int &i) {
        
        double ptemp    = (gamma_adiabat-1.)*(u.u3  - 0.5* u.u2*u.u2 / u.u1 );
        double pfinal   = pressure[i] - 0.5 * u.u1 * (phi_l - phi_r);
        
        /*
        if(switchi > 0) {
            pfinal = pressure[i] - dx[i] * omegaplus[i] * u.u1 * (phi[i-1] - phi[i]) / (dx[i-1] + dx[i]);
        }
        else {
            pfinal = pressure[i] - dx[i] * omegaminus[i] * u.u1 * (phi[i+1] - phi[i]) / (dx[i+1] + dx[i]);
        }*/
        
        
        //if self.debug > 2:
         //   print(" In P_HYDROSTATIC : phi_l-phi_r = " + repr(phi_l) + "-" + repr(phi_r) + " p,pfinal = " + repr([p,pfinal]))

        if (pfinal < 0.) {
            if(suppress_warnings == 0) {
                char a;
                cout<<"Warning: Negative pressure computed and replaced in cell "<<i<<" at time= "<<globalTime<<" iter = "<<steps<<endl;
                cout<<" cell "<<i;
                cout<<" u =  "<<u.u1<<"/"<<u.u2<<"/"<<u.u3;
                cout<<" phil-phir "<<(phi_l-phi_r)<<" p= "<<pressure[i]<< " ptemp= "<<ptemp<<" pfinal= "<<pfinal<<endl;
                cin>>a;
            }
                
            return pressure[i];
        }
            
        return pfinal;
    }
    
    
    
    //
    // p_hydrostatic_nonuniform - hydrostatic reconstruction, after Eq. A.4 in KM2016
    //
    // takes:
    // returns:
    //
    double hydro_run::get_p_hydrostatic_nonuniform(const int &i, const int &plusminus) {
        
        //double ptemp    = (gamma_adiabat-1.)*(u.u3  - 0.5* u.u2*u.u2 / u.u1 );
        double pfinal;   //= pressure[i] - 0.5 * u.u1 * (phi_l - phi_r);
        double phi_l = 0;
        double phi_r = 0;
        
        /*if(i==0)
            pfinal = pressure[0] - 0.5 * u[0].u1 * (phi[1] - phi[0]);
        
        else if(i==num_cells+1)
            pfinal = pressure[num_cells+1] - 0.5 * u[num_cells+1].u1 * (phi[num_cells] - phi[num_cells+1]);
        
        else {
            */
            if(plusminus == -1) {
                pfinal = pressure[i] - dx[i] * omegaplus[i] * u[i].u1 * (phi[i-1] - phi[i]) / (dx[i-1] + dx[i]);
            }
            else {
                pfinal = pressure[i] - dx[i] * omegaminus[i] * u[i].u1 * (phi[i+1] - phi[i]) / (dx[i+1] + dx[i]);
            }
            
        //}
        
        //if self.debug > 2:
         //   print(" In P_HYDROSTATIC : phi_l-phi_r = " + repr(phi_l) + "-" + repr(phi_r) + " p,pfinal = " + repr([p,pfinal]))

        if (pfinal < 0.) {
            if(suppress_warnings == 0) {
                cout<<"Warning: Negative pressure computed and replaced in nonuniform cell "<<i<<" at time= "<<globalTime<<" iter = "<<steps<<endl;
                cout<<" cell "<<i;
                cout<<" u =  "<<u[i].u1<<"/"<<u[i].u2<<"/"<<u[i].u3;
                cout<<" phil-phir "<<(phi_l-phi_r)<<" p= "<<pressure[i]<<" pfinal= "<<pfinal<<endl;
                char a ;
                cin>>a;
            }
                
            return pressure[i];
        }
            
        return pfinal;
    }

    //
    // source_grav - numerical source function, after Eq. 11 in KM2016
    //                        
    // takes:   state vector of rho,rhou,E, (hence R^vars) position (\in R^1)  and makes gravity, coriolis force or whatever from it 
    // returns: vector \in R^vars
    //
    AOS hydro_run::source_grav(AOS &u, int &j) {
        
        //if self.debug > 2:
        //    print(" In SOURCE_GRAV: phileft, phiright = " + repr([phileft, phiright]))
    
        //
        // The last regular cells get the uniform grid treatment
        // because omegas cannot be defined in those cells, without defining spacing for the ghost cells, which I dont want to do
        //
        if(j==1 || j==num_cells ) 
            return AOS(0, u.u1, u.u2) * (-1.) * (phi[j+1] - phi[j-1]) / (x_i12[j+1] - x_i12[j-1]); 
        else {
            
            return AOS(0, u.u1, u.u2) * (-1.) * ( omegaplus[j] * (phi[j] - phi[j-1]) / (dx[j] + dx[j-1]) + omegaminus[j] * (phi[j+1] - phi[j]) / (dx[j+1] + dx[j]) );
            
        }
        
        //return  AOS(0, u.u1, u.u2) * (-1.) * (phi[j+1] - phi[j-1]) / (x_i12[j+1] - x_i12[j-1]);
    }

    
        
        
        
    
        ///////////////////////////////////////////////////////////
        //
        // Chapter on radiative fluxes
        //
        //
        ///////////////////////////////////////////////////////////
        
    void hydro_run::update_radiation() {
        
        //
        // Update opacities in case thats needed.
        //
        //for(int i = 0; i < num_cells; i++) {
        //   opacity = some_function(density, temperature);
        //}
        
        opticaldepth[num_cells-1] = 1e-6; // Some start value for the top of the atmosphere. TODO: Improve this into some physically motivated estiamte.
        
        for(int i = num_cells-2; i>=0; i--)  {
            opticaldepth[i] = opticaldepth[i+1] + opacity[i] * u[i].u1 * (x_i12[i+1] - x_i12[i]);
            
            radiative_flux[i] = pow(temperature[i], 4.) / ( 3./4. * ( 2./3. + opticaldepth[i]) ); //Hubeny 1990 relation with zero scattering
            
        }
        
    }

    
    
    
    /*
     * THIS FUNCTION IS NOW IN COVID19 QUARANTINE BECAUSE ITS UTTER GARBAGE AND NOT USEFUEL
     * IT IS BAD AND WHOEVER WROTE IT SHOULD FEEL BAD
     * 
     * 
     * 
    //
    // This function calculates from the radiative fluxes the new internal energy, so that the total energy doesn't explode.
    // Or at least that's the intetion'
    //
        
    AOS hydro_run::source_radflux(int i) {
        
        double temp_eint = u[i].u3 - radiative_flux[i]; //Or something better.
        
        return AOS(0.,0.,temp_eint);
    }

    */

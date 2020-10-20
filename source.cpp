#include "aiolos.h"



    
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
    void c_Sim::init_grav_pot() {
        
        for(int i = 1; i <= num_cells; i++) {
            enclosed_mass[i] = planet_mass;      //No self-gravity
            phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            
        }
        phi[0]           = get_phi_grav(x_i12[1],         planet_mass);
            

    }

    //
    // Compute the gravity array as a self-gravitating solution
    //
    // takes: r as array
    //    
    void c_Sim::update_mass_and_pot() {
        
        double species_dens_sum;
        
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
            
            //TODO: PUTIN loop over all species densities
            species_dens_sum = 0.;
            for(int s = 0; s < num_species; s++)
                species_dens_sum += species[s].u[i].u1;
        
            
            
            enclosed_mass[i] = enclosed_mass[i-1] +  4. * 3.141592 * (pow(x_i[i],3.)-pow(x_i[i-1],3.) )/3. * species_dens_sum; //Straightfoward integration of the poisson eqn
            
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
    double c_Sim::get_phi_grav(double &r, double &mass) {
        
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
            return -mass*abs(domain_max - r);
        else
        {
            if(abs(r-planet_position) > rs_at_moment )
                return -mass/abs(r-planet_position);
            else
                return -mass * (pow(r-planet_position,3.)/pow(rs_at_moment,4.)  - 2.* pow(r-planet_position,2.)/pow(rs_at_moment,3.) + 2./rs_at_moment );
            
        }
            
    }    
    
    //
    // p_hydrostatic_nonuniform - hydrostatic reconstruction, after Eq. A.4 in KM2016
    //
    // takes:
    // returns:
    //
    double c_Species::get_p_hydrostatic_nonuniform(const int &i, AOS& u, const int &plusminus) {
        
        double pfinal;
    
        if(plusminus == -1) {
                pfinal = pressure[i] - base->dx[i] * base->omegaplus[i] * u.u1 * (base->phi[i-1] - base->phi[i]) / (base->dx[i-1] + base->dx[i]);
        } else {
                pfinal = pressure[i] - base->dx[i] * base->omegaminus[i] * u.u1 * (base->phi[i+1] - base->phi[i]) / (base->dx[i+1] + base->dx[i]);
        }
            
        //if self.debug > 2:
         //   print(" In P_HYDROSTATIC : phi_l-phi_r = " + repr(phi_l) + "-" + repr(phi_r) + " p,pfinal = " + repr([p,pfinal]))

        if (pfinal < 0.) {
            if(suppress_warnings == 0) {
                cout<<"Warning: Negative pressure computed and replaced in nonuniform cell "<<i<<" at time= "<<base->globalTime<<" iter = "<<base->steps<<endl;
                cout<<" cell "<<i<<" case /pm = "<<plusminus<<endl;
                cout<<" u =  "<<u.u1<<"/"<<u.u2<<"/"<<u.u3<<endl;
                cout<<" phil-phir "<<(base->phi[i+plusminus] - base->phi[i])<<" p= "<<pressure[i]<<" pfinal= "<<pfinal<<endl;
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
    AOS c_Species::source_grav(AOS &u, int &j) {
        
        //if self.debug > 2:
        //    print(" In SOURCE_GRAV: phileft, phiright = " + repr([phileft, phiright]))
    
        //
        // The last regular cells get the uniform grid treatment
        // because omegas cannot be defined in those cells, without defining spacing for the ghost cells, which I dont want to do
        //
        if(j==0 || j==num_cells+1) 
            return AOS(0, u.u1, u.u2) * (-1.) * (base->phi[j+1] - base->phi[j-1]) / (base->x_i12[j+1] - base->x_i12[j-1]); 
        else {
            
            return AOS(0, u.u1, u.u2) * (-1.) * ( base->omegaplus[j] * (base->phi[j] - base->phi[j-1]) / (base->dx[j] + base->dx[j-1]) + base->omegaminus[j] * (base->phi[j+1] - base->phi[j]) / (base->dx[j+1] + base->dx[j]) );
            
        }
        
    }

    
        ///////////////////////////////////////////////////////////
        //
        // Chapter on radiative fluxes
        //
        //
        ///////////////////////////////////////////////////////////
        
    void c_Species::update_radiation(std::vector<AOS>& u) {
        
        //
        // Update opacities in case thats needed.
        //
        //for(int i = 0; i < num_cells; i++) {
        //   opacity = some_function(density, temperature);
        //}
        
        opticaldepth[num_cells-1] = 1e-6; // Some start value for the top of the atmosphere. TODO: Improve this into some physically motivated estiamte.
        
        for(int i = num_cells-2; i>=0; i--)  {
            opticaldepth[i] = opticaldepth[i+1] + opacity[i] * u[i].u1 * (base->x_i12[i+1] - base->x_i12[i]);
            
            radiative_flux[i] = pow(temperature[i], 4.) / ( 3./4. * ( 2./3. + opticaldepth[i]) ); //Hubeny 1990 relation with zero scattering
            
        }
        
    }

    
    
        ///////////////////////////////////////////////////////////
        //
        // Chapter on frictional momentum exchange
        //
        //
        ///////////////////////////////////////////////////////////
        
    void c_Sim::compute_friction_step() {

        // Just do the basic update here
        for(int s=0; s < num_species; s++) {
            for(int j=0; j < num_cells+1; j++)
                species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
        }
        
        //Dummy function so far
        /*
        
        if(num_species == 2) {

            //Apply analytic solutions ...
        
        }
        else {
            
            //Build friction matrix
            
            //Iterate
            
            for(int s = 0; s < num_species; s++)
                //Add results to momenta to obtain new momenta
            
        }
        */
        
    }

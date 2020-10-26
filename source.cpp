#include <cassert>

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
    double c_Species::get_dp_hydrostatic(const int &i, const int &plusminus) {
        
        double dp_final;
    
        if(plusminus == -1) {
                dp_final = - base->dx[i] * base->omegaplus[i] * (base->phi[i-1] - base->phi[i]) / (base->dx[i-1] + base->dx[i]);
        } else {
                dp_final = - base->dx[i] * base->omegaminus[i] * (base->phi[i+1] - base->phi[i]) / (base->dx[i+1] + base->dx[i]);
        }
            
        return dp_final;
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
        assert(j > 0 && j <= num_cells) ;

        double dphidr_p = (base->phi[j] - base->phi[j-1]) / (base->dx[j] + base->dx[j-1]) ;
        double dphidr_m = (base->phi[j+1] - base->phi[j]) / (base->dx[j+1] + base->dx[j]) ;

        return AOS(0, u.u1, u.u2) * (-1.) * ( base->omegaplus[j] * dphidr_p  + base->omegaminus[j] * dphidr_m);
        
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
            
            radiative_flux[i] = pow(prim[i].temperature, 4.) / ( 3./4. * ( 2./3. + opticaldepth[i]) ); //Hubeny 1990 relation with zero scattering
            
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
        //for(int s=0; s < num_species; s++) {
        //    for(int j=0; j < num_cells+1; j++)
        //        species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
        //}
        
        
        
        if(num_species == 2 && friction_solver == 0) {
            
            if(debug > 0) cout<<"in friction, pos1"<<endl;

            //Apply analytic solutions ...
            double alpha=0, eps=0, f1=0, f2=0;
            double v1b=0, v2b=0, v1a=0, v2a=0, dv1=0, dv2=0;
            
            for(int j=0; j <= num_cells+1; j++){
                
                v1b = species[0].prim[j].speed;
                v2b = species[1].prim[j].speed;
                
                alpha = alpha_collision; //*(1e-2/x_i12[j]); // (f[0]+f[1])/(mu0*f[0]+mu1*f[1]) * k_b T/(m_i * b_i)
                eps   = species[0].u[j].u1 / species[1].u[j].u1;
                f1    = (1. + dt*alpha)/(1. + dt*alpha*(1.+eps)) ;
                f2    = dt*alpha*eps / (1. + dt * alpha);
                //f2    = 1. / (1. + 1./ (dt * alpha));
                
                v2a    = (v2b + v1b * f2 ) * f1;
                v1a    = v1b - (v2a - v2b) / eps;
                
                species[0].prim[j].speed = v1a;
                species[1].prim[j].speed = v2a;
            }
            
            if(debug > 0) cout<<"in friction, pos2"<<endl;
            
            for(int s=0; s<num_species; s++)
                species[s].eos->compute_conserved(&(species[s].prim[0]), &(species[s].u[0]), num_cells+2);        
            
            if(debug > 0) {
                char a;
                double dekin1 = 0.5*species[0].u[num_cells+1].u1 * (v1a - v1b) * (v1a - v1b) ; 
                double dekin2 = 0.5*species[1].u[num_cells+1].u1 * (v2a - v2b) * (v2a - v2b) ;
                double dvrel1 = (v1a - v1b)/v1b;
                double dvrel2 = (v2a - v2b)/v2b;
                cout<<"Relative differences in velocities 1/2 = "<<dvrel1<<" / "<<dvrel2<<" differences in Ekin = "<<dekin1<<" / "<<dekin2<<endl;
                cin>>a;
            }
        
        }
        else if(num_species > 1 && friction_solver == 1) {
            
            //friction_matrix_T   = Eigen::ArrayXXf::Zero(num_species, num_species);
            //friction_matrix_M   = Eigen::ArrayXXf::Zero(num_species, num_species);
            //friction_vec_input  = Eigen::ArrayXXf::Zero(num_species, 1);
            //friction_vec_output = Eigen::ArrayXXf::Zero(num_species, 1);
            
            for(int j=0; j <= num_cells+1; j++){
                
                if(debug > 0) cout<<"In friction, just starting cell ="<<j<<endl;
                
                // Construct input vector consisting of old velocities
                for(int si=0; si<num_species; si++)
                    friction_vec_input(si) = species[si].prim[j].speed;
                
                if(debug > 0) cout<<"    Before friciton coeffs."<<endl;
                
                // Compute or set friction coefficients
//                 friction_coefficients = Eigen::MatrixXd::Constant(num_species, num_species, alpha_collision); 
                
                for(int si=0; si<num_species; si++) 
                    for(int sj=0; sj<num_species; sj++)
                    {
                       if(si==sj)
                           friction_coefficients(si,sj) = 0.;
                       else if(si > sj)
                           friction_coefficients(si,sj) = alpha_collision;
                       else
                           friction_coefficients(si,sj) = alpha_collision  * species[sj].u[j].u1 / species[si].u[j].u1 ;
                        
                    }
                
                if(debug >= 0 && j==0 && steps == 0) cout<<"    Computed coefficient matrix ="<<endl<<friction_coefficients<<endl;
                if(debug >= 0 && j==0 && steps == 0) cout<<"    velocities ="<<endl<<friction_vec_input<<endl;
                if(debug >= 0 && j==0 && steps == 0) cout<<"    velocities[0] = "<<friction_vec_input(0)<<endl;
                if(debug >= 0 && j==0 && steps == 0) cout<<"    rho[0] = "<<species[0].u[j].u1<<endl;
                
                // Build friction matrices
                for(int si=0; si<num_species; si++) 
                    for(int sj=0; sj<num_species; sj++)
                    {
                        friction_matrix_M(si,sj) = - friction_coefficients(si,sj)*(1.-delta_ij(si,sj));
                        
                        for(int k=0; k<num_species; k++)
                            if(k != si)
                                friction_matrix_M(si,sj) += friction_coefficients(si,k) * delta_ij(si,sj);
                    }
                    
                if(debug >= 0 && j==0 && steps == 1) cout<<"    Final M ="<<endl<<friction_matrix_M<<endl;
                if(debug >= 0 && j==0 && steps == 1) cout<<"    ALPHA12 ="<<friction_coefficients(0,1)<<endl;
                if(debug >= 0 && j==0 && steps == 1) cout<<"    alpha21 ="<<friction_coefficients(1,0)<<endl;
                
                if(debug >= 0 && j==0 && steps == 1) cout<<"    M*v ="<<endl<<friction_matrix_M*friction_vec_input<<endl;
                
                // Build total matrix
                friction_matrix_T   = identity_matrix + friction_matrix_M * dt;
                
                if(debug >= 0 && j==0 && steps == 1) cout<<"    Final T ="<<endl<<friction_matrix_T<<endl;
                
                // Iterate
                friction_vec_output = friction_matrix_T.inverse() * friction_vec_input;
                
                if(debug >= 0 && j==0 && steps == 1) cout<<"    Final T inverse ="<<endl<<friction_matrix_T.inverse()<<endl;
                
                if(debug >= 0 && j==0 && steps == 1) cout<<"    v output ="<<endl<<friction_vec_output<<endl;
                
                if(debug > 0) {
                    cout<<"    Before final asignment."<<endl;
                    char a;
                    cin>>a;
                }
                
                
                for(int si=0; si<num_species; si++) {
                    species[si].prim[j].speed = friction_vec_output(si);
                }
                    
            }
            
            for(int si=0; si<num_species; si++) 
                species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
            
        }
        
        
    }

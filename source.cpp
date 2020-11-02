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
    
    
    //
    // Friction solver 0
    //
    
    void c_Sim::compute_friction_analytical() {

        if(debug > 0) cout<<"in analytic friction, num_species = "<<num_species<<endl;
        
        // Just do the basic update here
        //for(int s=0; s < num_species; s++) {
        //    for(int j=0; j < num_cells+1; j++)
        //        species[s].u[j] = species[s].u[j] + species[s].dudt[0][j]*dt ;
        //}
        
        if(num_species == 2) {
            
                

                //Apply analytic solutions ...
                double alpha=0, eps=0, f1=0, f2=0;
                double v1b=0, v2b=0, v1a=0, v2a=0;
                
                for(int j=0; j <= num_cells+1; j++){
                    
                    v1b = species[0].prim[j].speed;
                    v2b = species[1].prim[j].speed;
                    
                    alpha = alpha_collision * (1e-1/x_i12[j]); // (f[0]+f[1])/(mu0*f[0]+mu1*f[1]) * k_b T/(m_i * b_i)
                    eps   = species[0].u[j].u1 / species[1].u[j].u1;
                    f1    = (1. + dt*alpha)/(1. + dt*alpha*(1.+eps)) ;
                    f2    = dt*alpha*eps / (1. + dt * alpha);
                    
                    v2a    = (v2b + v1b * f2 ) * f1;
                    v1a    = v1b - (v2a - v2b) / eps;
                    
                    species[0].prim[j].speed = v1a;
                    species[1].prim[j].speed = v2a;
                    
                    species[0].prim[j].internal_energy += dt * alpha * (species[1].mass_amu/(species[0].mass_amu+species[1].mass_amu))  * pow( v1a - v2a, 2.);
                    species[1].prim[j].internal_energy += dt * alpha * eps * (species[0].mass_amu/(species[0].mass_amu+species[1].mass_amu))  * pow( v1a - v2a, 2.);
                    
                    
                }
                
                
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
        if(num_species == 3) {
                
                double det;
                double v1a=0, v2a=0, v3a=0, v1b=0, v2b=0, v3b=0;
                Eigen::Matrix3d a = Eigen::Matrix3d::Zero();
                Eigen::Vector3d dens_vector(0,0,0);
                
                if(debug > 0) cout<<"    Before radial loop. a ="<<a<<endl;
                
                for(int j=0; j <= num_cells+1; j++){
                
                    if(debug > 1) cout<<"    Before species loop."<<endl;
                    
                    for(int si=0; si<num_species; si++) {
                        dens_vector(si)        = species[si].u[j].u1; 
                    }
                    
                    v1b = species[0].prim[j].speed;
                    v2b = species[1].prim[j].speed;
                    v3b = species[2].prim[j].speed;
                    
                    if(debug > 1) cout<<"    Before coeff loop."<<endl;
                    
                    for(int si=0; si<num_species; si++) 
                        for(int sj=0; sj<num_species; sj++)
                        {
                        if(si==sj)
                            a(si,sj) = 0.;
                        else if(si > sj)
                            a(si,sj) = friction_coeff_mask(si,sj) * alpha_collision;
                        else
                            a(si,sj) = friction_coeff_mask(si,sj) * alpha_collision  * dens_vector(sj) / dens_vector(si);
                        }
                    
                    if(debug > 1) cout<<"    Before radial loop."<<endl;
                    /*
                    det = - ( a(1, 3) +dt* a(1, 3) *a(2, 1) + a(1, 2)* a(2, 3) +dt* a(1, 3)* a(2, 3)) * a(3, 1); 
                    det += (-a(1, 3) * a(2, 1) - a(2, 3) -dt* a(1, 2)* a(2, 3) - dt* a(1, 3)* a(2, 3)) * a(3, 2); 
                    det += (-a(1, 2) * a(2, 1) + (1 + dt* (a(1, 2) + a(1, 3)))* (1 +dt* (a(2, 1) + a(2, 3)))) *(1 +dt *(a(3, 1) + a(3, 2)));
                    
                    v1a =  v1b*(-a(2, 3) * a(3, 2) + (1 + dt* (a(2, 1) + a(2, 3))) * (1 +  dt * (a(3, 1) + a(3, 2))));
                    v1a += v2b*( a(1, 2) + dt * a(1, 2) * a(3, 1) + dt * a(1, 2) * a(3, 2) + a(1, 3)* a(3, 2));
                    v1a += v3b*(a(1, 3) + dt * a(1, 3)* a(2, 1) + a(1, 2) * a(2, 3) + dt * a(1, 3)* a(2, 3));
                    
                    v2a  = v1b*(a(2, 1) + dt * a(2, 1) * a(3, 1) + a(2, 3)* a(3, 1) + dt * a(2, 1) * a(3, 2));
                    v2a += v2b*( -a(1, 3) * a(3, 1) + (1 + dt * (a(1, 2)* + a(1, 3))) * (1 + dt * (a(3, 1) + a(3, 2) ) ))  ;
                    v2a += v3b*( a(1, 3) * a(2, 1) + a(2, 3) + dt * a(1, 2) * a(2, 3) + dt * a(1, 3) * a(2, 3));
                    
                    v3a  = v1b*(a(3, 1) + dt * a(2, 1) * a(3, 1) + dt * a(2, 3) * a(3, 1) + a(2, 1) * a(3, 2) );
                    v3a += v2b*(a(1, 2)* a(3, 1) + a(3, 2) + dt * a(1, 2) * a(3, 2) + dt * a(1, 3) * a(3, 2) );
                    v3a += v3b*(-a(1, 2)* a(2, 1) + (1 + dt* (a(1, 2) + a(1, 3)) ) * (1 + dt * (a(2, 1) + a(2, 3)) ) );
                    */
                    
                    det = - ( a(0, 2) + dt* a(0, 2) * a(1, 0) + a(0, 1) * a(1, 2) +dt* a(0, 2)* a(1, 2)) * a(2, 0); 
                    det += (-a(0, 2) * a(1, 0) - a(1, 2) -dt* a(0, 1)* a(1, 2) - dt* a(0, 2)* a(1, 2)) * a(2, 1); 
                    det += (-a(0, 1) * a(1, 0) + (1. + dt* (a(0, 1) + a(0, 2) ))* (1. +dt* (a(1, 0) + a(1, 2)))) *(1. +dt *(a(2, 0) + a(2, 1)));
 /*Det
    ([a, {1, 3}] + t [a, {1, 3}] [a, {2, 1}] + [a, {1, 2}] [a, {2, 3}] + t [a, {1, 3}] [a, {2, 3}]) * [a, {3, 1}] + 
   (-[a, {1, 3}] [a, {2, 1}] - [a, {2, 3}] - t [a, {1, 2}] [a, {2, 3}] - t [a, {1, 3}] [a, {2, 3}])     * [a, {3, 2}] + 
   (-[a, {1, 2}] [a, {2, 1}] + (1 + t ([a, {1, 2}] + [a, {1, 3}])) (1 + t ([a, {2, 1}] + [a, {2, 3}]))) * (1 + t ([a, {3, 1}] + [a, {3, 2}]))
                    */
                    
                    v1a =  v1b*(-a(1, 2) * a(2, 1) + (1. + dt* (a(1, 0) + a(1, 2))) * (1. +  dt * (a(2, 0) + a(2, 1))));
                    v1a += v2b*( a(0, 1) + dt * a(0, 1) * a(2, 0) + dt * a(0, 1) * a(2, 1) + a(0, 2)* a(2, 1));
                    v1a += v3b*(a(0, 2) + dt * a(0, 2)* a(1, 0) + a(0, 1) * a(1, 2) + dt * a(0, 2)* a(1, 2));
                    
                    v2a  = v1b*(a(1, 0) + dt * a(1, 0) * a(2, 0) + a(1, 2)* a(2, 0) + dt * a(1, 0) * a(2, 1));
                    v2a += v2b*( -a(0, 2) * a(2, 0) + (1. + dt * (a(0, 1)* + a(0, 2))) * (1. + dt * (a(2, 0) + a(2, 1) ) ))  ;
                    v2a += v3b*( a(0, 2) * a(1, 0) + a(1, 2) + dt * a(0, 1) * a(1, 2) + dt * a(0, 2) * a(1, 2));
                    
                    v3a  = v1b*(a(2, 0) + dt * a(1, 0) * a(2, 0) + dt * a(1, 2) * a(2, 0) + a(1, 0) * a(2, 1) );
                    v3a += v2b*(a(0, 1)* a(2, 0) + a(2, 1) + dt * a(0, 1) * a(2, 1) + dt * a(0, 2) * a(2, 1) );
                    v3a += v3b*(-a(0, 1)* a(1, 0) + (1. + dt* (a(0, 1) + a(0, 2)) ) * (1. + dt * (a(1, 0) + a(1, 2)) ) );
                    
                    v1a /= det;
                    v2a /= det;
                    v3a /= det;
                    
                    species[0].prim[j].speed = v1a;
                    species[1].prim[j].speed = v2a;
                    species[2].prim[j].speed = v3a;
                    
                    species[0].prim[j].internal_energy += dt * a(0,1) * (species[1].mass_amu/(species[0].mass_amu+species[1].mass_amu))  
                                                                     * pow( v1a - v2a, 2.);
                    species[0].prim[j].internal_energy += dt * a(0,2) * (species[2].mass_amu/(species[0].mass_amu+species[2].mass_amu))  
                                                                     * pow( v1a - v3a, 2.);
                                                                     
                    species[1].prim[j].internal_energy += dt * a(1,0) * (species[0].mass_amu/(species[1].mass_amu+species[0].mass_amu))  
                                                                     * pow( v1a - v2a, 2.);
                    species[1].prim[j].internal_energy += dt * a(1,2) * (species[2].mass_amu/(species[1].mass_amu+species[2].mass_amu))  
                                                                     * pow( v2a - v3a, 2.);
                    
                    species[2].prim[j].internal_energy += dt * a(2,0) * (species[0].mass_amu/(species[2].mass_amu+species[0].mass_amu))  
                                                                     * pow( v1a - v3a, 2.);
                    species[2].prim[j].internal_energy += dt * a(2,1) * (species[1].mass_amu/(species[2].mass_amu+species[1].mass_amu))  
                                                                     * pow( v2a - v3a, 2.);
                    ///TODO: Internal energy update
                    
            }
        
            if(debug > 0) {
                    char b;
                    double dekin1 = 0.5*species[0].u[num_cells+1].u1 * (v1a - v1b) * (v1a - v1b) ; 
                    double dekin2 = 0.5*species[1].u[num_cells+1].u1 * (v2a - v2b) * (v2a - v2b) ;
                    double dvrel1 = (v1a - v1b)/v1b;
                    double dvrel2 = (v2a - v2b)/v2b;
                    cout<<"alpha_collision ="<<alpha_collision<<" det = "<<det<<" a(0,0)="<<a(0,0)<<" matrix a = "<<a<<endl;
                    cout<<"Relative differences in velocities 1/2 = "<<dvrel1<<" / "<<dvrel2<<" differences in Ekin = "<<dekin1<<" / "<<dekin2<<endl;
                    cin>>b;
            }
            
    }
    
    if(debug > 0) cout<<"in friction, pos2"<<endl;
            
        for(int s=0; s<num_species; s++){
                species[s].eos->update_p_from_eint(&(species[s].prim[0]), num_cells+2);
                species[s].eos->compute_conserved(&(species[s].prim[0]), &(species[s].u[0]), num_cells+2);        
        }
        
}

//
// Friction solver 2
//
void c_Sim::compute_friction_numerical_dynamic() {
        
    if(debug > 0) cout<<"in numerical, dense friction, num_species = "<<num_species<<endl;
            
    double alpha_local;
    double coll_b;
    Eigen::VectorXd dens_vector(num_species);
    Eigen::VectorXd numdens_vector(num_species);
    Eigen::VectorXd mass_vector(num_species);
    
    for(int j=0; j <= num_cells+1; j++){
                
        for(int si=0; si<num_species; si++) {
                    friction_vec_input(si) = species[si].prim[j].speed;
                    dens_vector(si)        =  species[si].u[j].u1; 
                    numdens_vector(si)     =  species[si].prim[j].number_density; 
                    mass_vector(si)        =  species[si].mass_amu;
            }
        
        
                
        for(int si=0; si<num_species; si++) {
            for(int sj=0; sj<num_species; sj++) {
                    
                    if(collision_model == 'C') {
                        alpha_local = alpha_collision;
                    }
                    else {
                        coll_b      = 1./(1.4142*pow(numdens_vector(sj),0.3333333333333));     // Update this, if charged species collide
                        alpha_local = (numdens_vector(si) + numdens_vector(sj))/(numdens_vector(si)*mass_vector(sj) + numdens_vector(si)*mass_vector(sj) );
                        //cout<<"   in coll, alpha_local_prefactor = "<<alpha_local<<" ";
                        alpha_local *= kb * species[si].prim[j].temperature /(mass_vector(si) * coll_b);
                        //cout<<" alpha_local = "<<alpha_local<<" T = "<<species[si].prim[j].temperature<< " n="<<numdens_vector(si)<<" mass="<<mass_vector(si)<<" b = "<<coll_b<<" n^1/3 = "<<pow(numdens_vector3(sj),0.3333333333333)<<" n^-1/3"<<(1./pow(numdens_vector3(sj),0.3333333333333))<<endl;
                        
                    }
                    
                    if(si==0 && sj == 1)
                        alphas_sample(j) = alpha_local;
                            
                
                    if(si==sj)
                           friction_coefficients(si,sj) = 0.;
                    else if(si > sj)
                           friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local;
                    else
                           friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local  * dens_vector(sj) / dens_vector(si) ;
                    
                    if(si==0 && sj==3 && j==7 && steps == 1)
                        cout<<"    si/sj=0/3, alpha_local = "<<alpha_local<<" friction_coeff = "<<friction_coefficients(si,sj)<<" dens ratio "<<(dens_vector(sj) / dens_vector(si))<<endl;
                    //if(friction_coefficients(si,sj) < 1e-20)
                    //    friction_coefficients(si,sj) = 0.;
                    }
                    
                }
                
                if(debug >= 0 && j==7 && steps == 1) cout<<"    Computed coefficient matrix ="<<endl<<friction_coefficients<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    velocities ="<<endl<<friction_vec_input<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    velocities[0] = "<<friction_vec_input(0)<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    rho[0] = "<<species[0].u[j].u1<<endl;
                
                friction_matrix_M   = (friction_coefficients * unity_vector).asDiagonal();
                friction_matrix_M  -= friction_coefficients;
                friction_matrix_T   = identity_matrix + friction_matrix_M * dt;
                if(debug >= 0 && j==7 && steps == 1) cout<<"    T = "<<friction_matrix_T<<endl;
                
                //friction_vec_output = friction_matrix_T.partialPivLu().solve(friction_vec_input) ;
                LU.compute(friction_matrix_T) ;
                friction_vec_output = LU.solve(friction_vec_input);
                //friction_matrix_T   = friction_matrix_T.inverse();
                //friction_vec_output = friction_matrix_T * friction_vec_input;
                
                //friction_vec_output = friction_matrix_T.partialPivLu().solve(friction_vec_input) ;
                
                if(debug >= 0 && j==7 && steps == 1) cout<<"    T.inv = "<<friction_matrix_T<<endl;
                
                
                if(debug >= 0 && j==7 && steps == 1) cout<<"    v_inp = "<<friction_vec_input<<endl;
                if(debug >= 0 && j==7 && steps == 1) cout<<"    v_out = "<<friction_vec_output<<endl;
                
                //*/
                // Update new speed and internal energy
                //
                
                for(int si=0; si<num_species; si++)
                    species[si].prim[j].speed = friction_vec_output(si);
                
                for(int si=0; si<num_species; si++) {
                    double temp = 0;
                    
                    for(int sj=0; sj<num_species; sj++) {
                        
                        double tempsj = dt * friction_coefficients(si,sj) * (species[sj].mass_amu/(species[sj].mass_amu+species[si].mass_amu))  * pow( friction_vec_output(si) - friction_vec_output(sj), 2.);
                        
                        if(si != sj)
                            temp += tempsj;
                        
                        if(debug >= 0 && j==6 && steps == 1) cout<<"    e+tmpsj = "<<tempsj<<" at j = "<<j<<" for species si = "<<si<<" sj = "<<sj<<" dv = "<<friction_vec_output(si) - friction_vec_output(sj)<<" coeff = "<<friction_coefficients(si,sj)<<endl;
                        
                        if(debug >= 0 && j==7 && steps == 1) cout<<"    e+tmpsj = "<<tempsj<<" at j = "<<j<<" for species si = "<<si<<" sj = "<<sj<<" dv = "<<friction_vec_output(si) - friction_vec_output(sj)<<" coeff = "<<friction_coefficients(si,sj)<<endl;
                    }
                    species[si].prim[j].internal_energy += temp;
                    
                    
                    //
                    // Update scheme 2
                    //
                    
                    //update_cons_prim_after_friction(&species[si].u[j], &species[si].prim[j], friction_dEkin(si), friction_vec_output(si), species[si].mass_amu, species[si].gamma_adiabat, species[si].cv);
                    
                }
                
            }
            
            //
            // Update scheme 1
            //
            for(int si=0; si<num_species; si++) {
                species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
                species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
            }
    
}
    
//
// Friction solver 1
//

//EIGEN_DONT_INLINE 
 void c_Sim::compute_friction_numerical_static() {
        
    if(debug > 0) cout<<"in numerical, dense friction8, num_species = "<<num_species<<endl;
    
    //if(num_species==12) {
        double alpha_local;
        double coll_b;
        
        for(int j=1; j <= num_cells; j++){
                    
                for(int si=0; si<num_species; si++) {
                        friction_vec_input3(si) = species[si].prim[j].speed;
                        dens_vector3(si)        =  species[si].u[j].u1; 
                        numdens_vector3(si)     =  species[si].prim[j].number_density; 
                        mass_vector3(si)        =  species[si].mass_amu;
                }
                
                    
                for(int si=0; si<num_species; si++) 
                    for(int sj=0; sj<num_species; sj++)
                    {
                        
                        if(collision_model == 'C') {
                            alpha_local = alpha_collision;
                        }
                        else {
                            coll_b      = 1./(1.4142*numdens_vector3(sj));     // Update this, if charged species collide
                            alpha_local = (numdens_vector3(si) + numdens_vector3(sj))/(numdens_vector3(si)*mass_vector3(sj) + numdens_vector3(si)*mass_vector3(sj) );
                            alpha_local *= kb * species[si].prim[j].temperature /(mass_vector3(si) * coll_b);
                        }
                        
                        
                        if(si==0 && sj == 1)
                            alphas_sample(j) = alpha_local;
                            
                        if(si==sj)
                            friction_coefficients3(si,sj) = 0.;
                        else if(si > sj)
                            friction_coefficients3(si,sj) = friction_coeff_mask(si,sj) * alpha_local;
                        else
                            friction_coefficients3(si,sj) = friction_coeff_mask(si,sj) * alpha_local  * dens_vector3(sj) / dens_vector3(si);
                    }
                    
                if(debug > 0) cout<<"   before matrix operations. "<<endl;
                    
                //friction_matrix_M = friction_coefficients.rowwise().sum().asDiagonal();
                friction_matrix_M3   = (friction_coefficients3 * unity_vector3).asDiagonal();
                friction_matrix_M3  -= friction_coefficients3;
                    
                friction_matrix_T3   = identity_matrix3;
                friction_matrix_T3  += dt * friction_matrix_M3 ;
                    
                if(debug > 1) cout<<"    in friction 8, j = "<<j<<" before inverse, T = "<<endl<<friction_matrix_T3<<endl<<" M = "<<endl<<friction_matrix_M3<<endl<<" A = "<<endl<<friction_coefficients3<<endl<<" v input = "<<friction_vec_input3<<endl;
                    
                //friction_vec_output3.noalias() = friction_matrix_T3.inverse().eval() * friction_vec_input3;
                //friction_matrix_T3   = friction_matrix_T3.inverse().eval();
                //friction_vec_output3.noalias() = friction_matrix_T3 * friction_vec_input3;
                //friction_vec_output3 = friction_matrix_T3 * friction_vec_input3;
                //friction_vec_output3           = friction_vec_input3;
                
                LU3.compute(friction_matrix_T3) ;
                friction_vec_output3 = LU3.solve(friction_vec_input3);
                
                for(int si=0; si<num_species; si++) {
                    double temp = 0;
                    
                    for(int sj=0; sj<num_species; sj++) {
                        temp += dt * friction_coefficients3(si,sj) * (species[sj].mass_amu/(species[sj].mass_amu+species[si].mass_amu))  * pow( friction_vec_output3(si) - friction_vec_output3(sj), 2.);
                    }
                    species[si].prim[j].internal_energy += temp;
                }
            
            
        } //End radial loop
            //
            // Update scheme 1
            //
            for(int si=0; si<num_species; si++) {
                species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
                species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
            }
        
    //}
    
}
    
    
        ///////////////////////////////////////////////////////////
        //
        // Chapter on frictional momentum exchange
        //
        //
        ///////////////////////////////////////////////////////////
        
    void c_Sim::compute_friction_numerical_sparse() {

        if(debug > 0) cout<<"in numerical, sparse friction, num_species = "<<num_species<<endl;
            
            double alpha_local;
            Eigen::VectorXd dens_vector(num_species);
            for(int j=0; j <= num_cells+1; j++){
                
                for(int si=0; si<num_species; si++) {
                    friction_vec_input(si) = species[si].prim[j].speed;
                    dens_vector(si)        =  species[si].u[j].u1; 
                }
                    
                
                //if(debug > 0) cout<<"    Before friciton coeffs."<<endl;
                
                // Compute or set friction coefficients
//                 friction_coefficients = Eigen::MatrixXd::Constant(num_species, num_species, alpha_collision); 
                alpha_local = alpha_collision * (1e-1/x_i12[j]);
                
                for(int si=0; si<num_species; si++) 
                    for(int sj=0; sj<num_species; sj++)
                    {
                       if(si==sj)
                           friction_coefficients(si,sj) = 0.;
                       else if(si > sj)
                           friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local;
                       else
                           //friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local  * dens_vector2(sj,j) / dens_vector2(si,j) ;
                           friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local  * dens_vector(sj) / dens_vector(si) ;
                       //friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local  * species[sj].u[j].u1 / species[si].u[j].u1 ;
                       
                        
                    }
                
                if(debug >= 1 && j==0 && steps == 0) cout<<"    Computed coefficient matrix ="<<endl<<friction_coefficients<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    velocities ="<<endl<<friction_vec_input<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    velocities[0] = "<<friction_vec_input(0)<<endl;
                if(debug >= 1 && j==0 && steps == 0) cout<<"    rho[0] = "<<species[0].u[j].u1<<endl;
                
                //Eigen::MatrixXd friction_matrix_M2(num_species,num_species);
                //friction_matrix_M = friction_coefficients.rowwise().sum().asDiagonal();
                friction_matrix_M = (friction_coefficients * unity_vector).asDiagonal();
                friction_matrix_M.noalias()          -= friction_coefficients;
                
                if(debug >= 1 && j==0 && steps == 1) {
                    cout<<"    Final M ="<<endl<<friction_matrix_M<<endl;
                    //cout<<"    Final M2 ="<<endl<<friction_matrix_M2<<endl;
                }
                    
                if(debug >= 1 && j==0 && steps == 1) cout<<"    ALPHA12 ="<<friction_coefficients(0,1)<<endl;
                if(debug >= 1 && j==0 && steps == 1) cout<<"    alpha21 ="<<friction_coefficients(1,0)<<endl;
                if(debug >= 1 && j==0 && steps == 1) cout<<"    M*v ="<<endl<<friction_matrix_M*friction_vec_input<<endl;
                
                // Build total matrix
                friction_matrix_T.noalias()           = identity_matrix + friction_matrix_M * dt;
                
                if(debug >= 1 && j==0 && steps == 1) cout<<"    Final T ="<<endl<<friction_matrix_T<<endl;
                
                // Iterate
                friction_vec_output.noalias()         = friction_matrix_T.inverse() * friction_vec_input;
                
                if(debug >= 1 && j==0 && steps == 1) cout<<"    Final T inverse ="<<endl<<friction_matrix_T.inverse()<<endl;
            
                if(debug >= 1 && j==0 && steps == 1) cout<<"    v output ="<<endl<<friction_vec_output<<endl;
                
                if(debug > 0) {
                    cout<<"    Before final asignment."<<endl;
                    char a;
                    cin>>a;
                }
                
                //
                // Update new speed and internal energy
                //
                
                for(int si=0; si<num_species; si++)
                    species[si].prim[j].speed = friction_vec_output(si);
                
                for(int si=0; si<num_species; si++) {
                    friction_dEkin(si) = 0.;
                    
                    for(int sj=0; sj<num_species; sj++) {
                            double temp = 0;
                            
                            //if(si != sj){
                                 temp += dt * friction_coefficients(si,sj) * (species[sj].mass_amu/(species[sj].mass_amu+species[si].mass_amu))  * pow( friction_vec_output(si) - friction_vec_output(sj), 2.);
                                
                                //friction_dEkin(si) += dt * species[si].u[j].u1 * friction_coefficients(si,sj)  * ( friction_vec_output(si) - friction_vec_output(sj)) * (species[si].mass_amu*friction_vec_output(si) +species[sj].mass_amu*friction_vec_output(sj) ) /(species[sj].mass_amu+species[si].mass_amu);
                            //}
                            
                            species[si].prim[j].internal_energy += temp;
                    }
                    
                    
                    //
                    // Update scheme 2
                    //
                    
                    //update_cons_prim_after_friction(&species[si].u[j], &species[si].prim[j], friction_dEkin(si), friction_vec_output(si), species[si].mass_amu, species[si].gamma_adiabat, species[si].cv);
                    
                }
                
            }
            
            //
            // Update scheme 1
            //
            for(int si=0; si<num_species; si++) {
                species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
                species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
            }
                
            
        
        
        
}

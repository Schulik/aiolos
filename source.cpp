///////////////////////////////////////////////////////////
//
//
//  source.cpp
//
// This file contains routines pertaining to the gravitational and friction sources.
//
//
//
///////////////////////////////////////////////////////////
#define EIGEN_RUNTIME_NO_MALLOC


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
        
        for(int i = 1; i <= num_cells+1; i++) {
            enclosed_mass[i] = planet_mass;      //No self-gravity
            phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            
            if(use_tides == 1)
                phi[i] -=0.5* 3.*G*star_mass*x_i12[i]*x_i12[i]/pow(planet_semimajor*au,3.);
            
        }
        phi[0]           = get_phi_grav(x_i12[1],         planet_mass);
        //phi[1]           = get_phi_grav(x_i12[2],         planet_mass);
        
        rhill = planet_semimajor*au * std::pow(planet_mass / 3. / star_mass, 0.333333333333333333333);
    }

    //
    // Compute the gravity array as a self-gravitating solution
    //
    // takes: r as array
    //    
    void c_Sim::update_mass_and_pot() {
        
        double species_dens_sum;
        
        if(debug >= 1)
            cout<<"In update mass. "<<endl;
        
        enclosed_mass[0] = planet_mass;
        enclosed_mass[1] = planet_mass;
        
        for(int i = 2; i <= num_cells; i++) {
            
            if(debug >= 1) {
                if(i>=0)
                    cout<<"In update mass in "<<i<<"th element."<<endl;
                if(i==num_cells-1)
                    cout<<"In update mass, last element."<<endl;
            }
            
            //Loop over all species densities
            species_dens_sum = 0.;
            for(int s = 0; s < num_species; s++)
                species_dens_sum += species[s].u[i].u1;
            
            double xn3 = x_i[i];
            xn3 = xn3*xn3*xn3;
            double xl3 = x_i[i-1];
            xl3 = xl3*xl3*xl3;
            
            enclosed_mass[i] = enclosed_mass[i-1] +  4. * 3.141592 * (xn3 - xl3)/3. * species_dens_sum; //Straightfoward integration of the poisson eqn
            
            if (use_self_gravity == 1) 
                phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            else
                phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[0]);
            
            //if use tidal gravity
            if(use_tides == 1) {
                double d3 = planet_semimajor*au;
                d3 = d3*d3*d3;
                phi[i] -= 0.5* 3.*G*star_mass*x_i12[i]*x_i12[i]/d3;
            }
                
            //            -G*mass/abs(r-planet_position);
        }
        
        if(debug >= 1)
            cout<<"Done in update mass. "<<endl;
    }
    
    //
    // Initialize the gravity array with a non-self gravitating solution
    //
    // takes: r as single value
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
            return -G*mass*abs(domain_max - r);
        else
        {
            if(abs(r-planet_position) > rs_at_moment )
                return -G*mass/abs(r-planet_position);
            else
                return -G*mass * (pow(r-planet_position,3.)/pow(rs_at_moment,4.)  - 2.* pow(r-planet_position,2.)/pow(rs_at_moment,3.) + 2./rs_at_moment );
            
        }
            
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

        double dphidr_p = (phi_s[j] - phi_s[j-1]) / (base->dx[j] + base->dx[j-1]) ;
        //double dphidr_p = (base->phi[j]*K_zzf[j] - base->phi[j-1]*K_zzf[j-1]) / (base->dx[j] + base->dx[j-1]) ;
        double dphidr_m = (phi_s[j+1] - phi_s[j]) / (base->dx[j+1] + base->dx[j]) ;
        //double dphidr_m = (base->phi[j+1]*K_zzf[j+1] - base->phi[j]*K_zzf[j]) / (base->dx[j+1] + base->dx[j]) ;

        return AOS(0, u.u1, u.u2) * (-1.) * ( base->omegaplus[j] * dphidr_p  + base->omegaminus[j] * dphidr_m);
        
    }

    

    ///////////////////////////////////////////////////////////
    //
    // Chapter on frictional momentum exchange
    //
    //
    ///////////////////////////////////////////////////////////
    
void c_Sim::compute_drag_update() {
        if(alpha_collision > 0 && num_species > 1) {
            if(friction_solver == 0)
                compute_friction_analytical();
            else
                compute_friction_numerical();
        }
    }
    
    //
    // Friction solver 0
    //
    
void c_Sim::compute_friction_analytical() {

        if(debug > 0) cout<<"in analytic friction, num_species = "<<num_species<<endl;
        
        if(num_species == 2) {
        
            //Apply analytic solutions ...
            double alpha=0, eps=0, f1=0, f2=0;
            double v1b=0, v2b=0, v1a=0, v2a=0;
            double coll_b;
            double mumass, meanT, m0, m1;
            
            for(int j=0; j <= num_cells+1; j++){
                
                v1b = species[0].prim[j].speed;
                v2b = species[1].prim[j].speed;
                
                if(collision_model == 'C') {
                    alpha = alpha_collision;
                } else {
                    // Physical model       
                    m0  = species[0].mass_amu*amu;
                    m1  = species[1].mass_amu*amu;
                    mumass = m0 * m1 / (m0 + m1);
                    meanT  = (m1 * species[0].prim[j].temperature + m0 * species[1].prim[j].temperature) / (m0 + m1);             
                    if (species[0].is_dust_like || species[1].is_dust_like) {
                        // Use Epstein drag law (Hard sphere model)
                        double RHO_DUST = 1;
                        double s0=0, s1=0 ;
                        if (species[0].is_dust_like)
                            s0 = std::pow(3*m0/(4*M_PI*RHO_DUST), 1/3.) ;
                        if (species[0].is_dust_like)
                            s1 = std::pow(3*m1/(4*M_PI*RHO_DUST), 1/3.) ;
                        
                        double A = M_PI*(s0+s1)*(s0+s1) ;
                        double v_th = std::sqrt(8*kb*meanT/(M_PI * mumass)) ;

                        alpha = 4/3. * species[1].prim[j].density * v_th * A/(m0+m1) ;
                    } else {
                        // Gas-gas collisions. Used binary diffusion coefficient                        
                        //ntot = species[0].prim[j].number_density + species[1].prim[j].number_density;
                        //fi   = species[0].prim[j].number_density / ntot;
                        //fj   = species[1].prim[j].number_density / ntot;
                        //mtot = species[0].mass_amu*amu + species[1].mass_amu*amu;
                        //mui  = species[0].mass_amu*amu / mtot;
                        //muj  = species[1].mass_amu*amu / mtot;
                        
                        coll_b      = 5.0e17 * std::pow(meanT, 0.75) ;     // from Zahnle & Kasting 1986 Tab. 1

                        alpha = kb * meanT * species[1].prim[j].number_density / (m0 * coll_b) ; // From Burgers book, or Schunk & Nagy.
                        //alpha_local = kb * meanT /(mumass * coll_b); 
                            
                        //alpha *= (fi + fj)/(fi * muj + fj * mui);
                        alpha *= alpha_collision;
                    }
                }
                
                alphas_sample(j) = alpha;
                
                friction_sample(j) = alphas_sample(j) * (v2b - v1b);
                
                //alpha = alpha_collision * (1e-1/x_i12[j]); // (f[0]+f[1])/(mu0*f[0]+mu1*f[1]) * k_b T/(m_i * b_i)
                eps   = species[0].u[j].u1 / species[1].u[j].u1;
                f1    = (1. + dt*alpha)/(1. + dt*alpha*(1.+eps)) ;
                f2    = dt*alpha*eps / (1. + dt * alpha);
                
                v2a    = (v2b + v1b * f2 ) * f1;
                v1a    = v1b - (v2a - v2b) / eps;
                
                species[0].prim[j].speed = v1a;
                species[1].prim[j].speed = v2a;
                
                //if(debug >= 1 && j==1200 && steps == 10) {
                if(debug >= 1) {
                    //cout<<" v1b, v2b = "<<v1b<<" "<<v2b<<endl<<" v1a, v2a = "<<v1a<<" "<<v2a<<endl;
                    //cout<<" dv1, dv2 = "<<v1a-v1b<<" "<<v2a-v2b<<endl;
                    //cout<<" dp1, dp2 = "<<(v1a-v1b)*species[0].u[j].u1 + (v2a-v2b)*species[1].u[j].u1<<endl;
                    cout<<"    dt*alpha = "<<dt*alpha<<" alpha = "<<alpha<<" eps = "<<species[1].u[j].u1/species[0].u[j].u1<<" dp/|p| = "<<((v1a-v1b)*species[0].u[j].u1 + (v2a-v2b)*species[1].u[j].u1)/(std::sqrt(species[0].u[j].u2*species[0].u[j].u2) + std::sqrt(species[1].u[j].u2*species[1].u[j].u2) )<<" meanT = "<<meanT<<" coll_b = "<<coll_b<<" cell = "<<j<<endl;
                    
                    //char a;
                    //cin>>a;
                }
                
                double v_half = 0.5*(v1a + v1b) - 0.5*(v2a + v2b) ;
                double v_end  = v1a - v2a ;

                species[0].prim[j].internal_energy += dt * alpha * (species[1].mass_amu/(species[0].mass_amu+species[1].mass_amu)) * v_half * v_end ;
                species[1].prim[j].internal_energy += dt * alpha * eps * (species[0].mass_amu/(species[0].mass_amu+species[1].mass_amu)) * v_half * v_end ;
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
                
                det = - ( a(0, 2) + dt* a(0, 2) * a(1, 0) + a(0, 1) * a(1, 2) +dt* a(0, 2)* a(1, 2)) * a(2, 0); 
                det += (-a(0, 2) * a(1, 0) - a(1, 2) -dt* a(0, 1)* a(1, 2) - dt* a(0, 2)* a(1, 2)) * a(2, 1); 
                det += (-a(0, 1) * a(1, 0) + (1. + dt* (a(0, 1) + a(0, 2) ))* (1. +dt* (a(1, 0) + a(1, 2)))) *(1. +dt *(a(2, 0) + a(2, 1)));
                
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
    
    if(debug > 0) {
        cout<<"in friction, pos2 num_cells = "<<num_cells<<endl;  
        char a;
        cin>>a;
    } 
            
    for(int s=0; s<num_species; s++){
            species[s].eos->update_p_from_eint(&(species[s].prim[0]), num_cells+2);
            species[s].eos->compute_conserved(&(species[s].prim[0]), &(species[s].u[0]), num_cells+2);        
    }
        
}

//
// Friction solver 2
//
void c_Sim::compute_friction_numerical() {
    
    Eigen::internal::set_is_malloc_allowed(false) ;
         
    if(debug > 0) cout<<"in numerical, dense friction, num_species = "<<num_species<<endl;
    
    for(int j=0; j <= num_cells+1; j++){

        fill_alpha_basis_arrays(j);
        compute_alpha_matrix(j);
                
        if(debug >= 1 && j==1200 && steps == 10)  {
            cout<<"    Computed coefficient matrix ="<<endl<<friction_coefficients<<endl;
            cout<<"    velocities ="<<endl<<friction_vec_input<<endl;
            cout<<"    velocities[0] = "<<friction_vec_input(0)<<endl;
            cout<<"    rho[0] = "<<species[0].u[j].u1<<endl;
        }
            
            
                
        friction_matrix_T = identity_matrix - friction_coefficients * dt;
        friction_matrix_T.diagonal().noalias() += dt * (friction_coefficients * unity_vector);
        
        LU.compute(friction_matrix_T) ;
        friction_vec_output.noalias() = LU.solve(friction_vec_input);
        
        if(debug >= 1 && j==700 && steps == 10) {
            
            cout<<"    T = "<<endl<<friction_matrix_T<<endl;
            cout<<" a10/a01 = "<<friction_matrix_T(0,1)/friction_matrix_T(1,0)<<" dens(0)/dens(1) = "<<dens_vector(1)/dens_vector(0)<<" dt*a*eps = "<<friction_coefficients(0,1) * dt * dens_vector(0)/dens_vector(1)<<endl;
            cout<<"    v_inp = "<<endl<<friction_vec_input<<endl;
            cout<<"    v_out = "<<endl<<friction_vec_output<<endl;
            cout<<"    dv_ = "<<endl<<friction_vec_output-friction_vec_input<<endl;
            cout<<"    dp  = "<<(friction_vec_output(0)-friction_vec_input(0))*dens_vector(0)<<" + "<<(friction_vec_output(1)-friction_vec_input(1))*dens_vector(1)<<" = "<<(friction_vec_output(0)-friction_vec_input(0))*dens_vector(0) + (friction_vec_output(1)-friction_vec_input(1))*dens_vector(1)<<endl;
            cout<<"      friction_coeff ="<<endl<<friction_matrix_T<<endl;
            cout<<"      friction_coeff ="<<endl<<friction_matrix_T(0,1)/dt<<endl;
            cout<<" dt * alpha01 = "<<endl<<friction_matrix_T(0,1)<<endl;
            cout<<" dt * alpha01 * eps = "<<endl<<friction_matrix_T(0,1) * dens_vector(0)/dens_vector(1)<<endl;
            char a;
            cin>>a;
        }
      
        //*/
        // Update new speed and internal energy
        //
        
        for(int si=0; si<num_species; si++)
            species[si].prim[j].speed = friction_vec_output(si);
        
        for(int si=0; si<num_species; si++) {
            double temp = 0;
            
            for(int sj=0; sj<num_species; sj++) {
                double v_end = friction_vec_output(si) - friction_vec_output(sj) ;
                double v_half = 0.5*(v_end + friction_vec_input(si) - friction_vec_input(sj)) ; 
                                               
                temp +=  dt * friction_coefficients(si,sj) * (species[sj].mass_amu/(species[sj].mass_amu+species[si].mass_amu)) * v_half * v_end ;
                
                if(si==0 && sj == 1) {
                    friction_sample(j) = alphas_sample(j) * (friction_vec_input(si) - friction_vec_input(sj));
                }
            }
            species[si].prim[j].internal_energy += temp;
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


void c_Sim::fill_alpha_basis_arrays(int j) { //Called in compute_friction() in source.cpp
    
    for(int si=0; si<num_species; si++) {
        friction_vec_input(si) =  species[si].prim[j].speed;
        dens_vector(si)        =  species[si].u[j].u1; 
        numdens_vector(si)     =  species[si].prim[j].number_density; 
        mass_vector(si)        =  species[si].mass_amu*amu;
        temperature_vector(si) =  species[si].prim[j].temperature;
        temperature_vector_augment(si) = species[si].prim[j].temperature - dt * (species[si].dS(j) + species[si].dG(j)) / species[si].u[j].u1 / species[si].cv;
    }
}
    
void c_Sim::compute_alpha_matrix(int j) { //Called in compute_friction() and compute_radiation() in source.cpp
        
        double alpha_local;
        double coll_b;
        //double ntot;
        double mtot;
        double mumass, mumass_amu, meanT;
        
        for(int si=0; si<num_species; si++) {
            for(int sj=0; sj<num_species; sj++) {
                    
                    if(collision_model == 'C') {
                        alpha_local = alpha_collision;
                        if (si > sj) alpha_local *= dens_vector(sj) / dens_vector(si) ;
                    }
                    else {
                        // Physical drag law
                        mtot = mass_vector(si) + mass_vector(sj);
                        //mui  = mass_vector(si) / mtot;
                        //muj  = mass_vector(sj) / mtot;
                        mumass = mass_vector(si) * mass_vector(sj) / (mass_vector(si) + mass_vector(sj));
                        mumass_amu = mumass/amu;
                        meanT  = (mass_vector(sj)*temperature_vector(si) + mass_vector(si)*temperature_vector(sj)) / (mass_vector(si) + mass_vector(sj)); //Mean collisional mu and T from Schunk 1980
                        
                        if (species[si].is_dust_like || species[sj].is_dust_like) {
                            // Use Epstein drag law (Hard sphere model)
                            double RHO_DUST = 1;
                            double s0=0, s1=0 ;
                            if (species[si].is_dust_like)
                                s0 = std::pow(3*mass_vector(si)/(4*M_PI*RHO_DUST), 1/3.) ;
                            if (species[sj].is_dust_like)
                                s1 = std::pow(3*mass_vector(sj)/(4*M_PI*RHO_DUST), 1/3.) ;
                            
                            double A = M_PI*(s0+s1)*(s0+s1) ;
                            double v_th = std::sqrt(8*kb*meanT/(M_PI * mumass)) ;

                            alpha_local = 4/3. * dens_vector(sj) * v_th * A / mtot ;
                        } else {
                            // Gas-gas collisions. Used binary diffusion coefficient                        

                            //ntot = numdens_vector(si) + numdens_vector(sj);
                            //fi   = numdens_vector(si) / ntot;
                            //fj   = numdens_vector(sj) / ntot;
                            
                            int ci = (int)(species[si].static_charge*1.01);
                            int cj = (int)(species[sj].static_charge*1.01);
                            double qn = (std::fabs(species[si].static_charge) + std::fabs(species[sj].static_charge));
                            string ccase = "";
                            
                            if(std::abs(ci) + std::abs(cj) == 0) { //n-n collision
                                coll_b      = 5.0e17 * std::sqrt(std::sqrt(meanT*meanT*meanT)) ;     //std::pow(meanT, 0.75) // from Zahnle & Kasting 1986 Tab. 1
                                alpha_local = kb * meanT * numdens_vector(sj) / (mass_vector(si) * coll_b) ; // From Burgers book, or Schunk & Nagy.
                                ccase = " n-n ";
                            }
                            else if(std::abs(ci) == 0 || std::abs(cj) == 0) {//i-n collision
                                
                                if(species[si].mass_amu < 1.1 && species[sj].mass_amu < 1.1) {  //Resonant H,H+ collisions
                                    alpha_local = 2.65e-10 * numdens_vector(sj) * std::sqrt(meanT) * std::pow((1. - 0.083 * std::log10(meanT)), 2.);
                                    
                                    ccase = "resonant i-n ";
                                    
                                } else {

                                    alpha_local = 1*2.21 * 3.141592 * numdens_vector(sj) * mass_vector(sj)/(mass_vector(sj)+mass_vector(si));
                                    alpha_local *= std::sqrt(0.66 / mumass ) * 1e-12 * qn * elm_charge; //0.66 is the polarizability of neutral atomic H                                        
                                    
                                    ccase = "nonresonant i-n ";
                                }
                            }
                            else { //i-i collision
                                alpha_local = 1.27 * species[si].static_charge * species[si].static_charge * species[sj].static_charge * species[sj].static_charge * std::sqrt(mumass_amu) * amu/mass_vector(si);
                                alpha_local *= numdens_vector(sj) / std::sqrt(meanT*meanT*meanT);
                                
                                ccase = " i-i ";
                            }
                            
                            //if(si!=sj)
                            //    cout<<"cell "<<j<<" colliding "<<species[si].speciesname<<" q="<<species[si].static_charge<<" with "<<species[sj].speciesname<<" q="<<species[sj].static_charge<<" case = "<<ccase<<" alpha/n = "<<alpha_local/ numdens_vector(sj)<<" n_sj = "<<numdens_vector(sj)<<endl;
                            
                            //if(photochemistry_level == 1 && (si + sj) == 3)
                            alpha_local *= alpha_collision;
                        }
                    }
                    
                    if(si==0 && sj == 1) {
                        alphas_sample(j) = alpha_local;
                        //cout<<"    spec "<<species[si].name<<" j = "<<j<<" alpha_local = "<<alpha_local<<endl;
                        if(debug > 0 && steps == 1 && ((j==2) || (j==num_cells-2) || (j==num_cells/2)) )
                            cout<<"    spec "<<species[si].speciesname<<" j = "<<j<<" alpha_local = "<<alpha_local<<endl;
                    }
                    
                    // Fill alpha_ij and alpha_ji
                    if(si==sj)
                       friction_coefficients(si,sj) = 0.;
                    else
                       friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local;

                   
                   if(debug >= 1 && j==700 && steps == 10) {
                        cout<<" debug alpha matrix(si="<<si<<",sj="<<sj<<") = "<<friction_coefficients(si,sj)<<" mask = "<<friction_coeff_mask(si,sj)<<" alpha_local = "<<alpha_local<<" 1/eps = "<<(dens_vector(si) / dens_vector(sj))<<" eps ="<<(dens_vector(sj) / dens_vector(si))<<endl;
                   }
              
            }
        }
        //char a;
        //cin>>a;
}


void c_Sim::compute_collisional_heat_exchange_matrix(int j) {
    
    // Get the alpha matrix
    compute_alpha_matrix(j) ;

    // Convert to collision matrix
    for(int si=0; si<num_species; si++) {
        double diag_sum = 0 ;
        for(int sj=0; sj<num_species; sj++) {
            friction_coefficients(si, sj) *= 
                3 * kb / (species[si].cv * (mass_vector(si) + mass_vector(sj))) ;

            diag_sum += friction_coefficients(si, sj) ;
        }
        // Set the diagonal to sum of Ti terms. 
        friction_coefficients(si, si) -= diag_sum ;
    }

}

void c_Sim::compute_collisional_heat_exchange() {
    Eigen::internal::set_is_malloc_allowed(false) ;
    
    if (num_species == 1)
        return ;

    if(debug > 0) 
        cout << "in compute_collisional_heat_exchange, num_species = "
             << num_species << endl;
    
    for(int j=0; j <= num_cells+0; j++){

        fill_alpha_basis_arrays(j);
        compute_collisional_heat_exchange_matrix(j);

        // Solve implicit equation for new temperature
        friction_matrix_T = identity_matrix - friction_coefficients * dt;
        
        LU.compute(friction_matrix_T) ;
        friction_vec_input.noalias() = LU.solve(temperature_vector_augment);
        
        // Set friction_vec_output to dT/dt:
        friction_vec_output = 
            friction_coefficients * friction_vec_input ;

        //
        // Update internal energy (dE = Cv dT/dt * dt)
        //
        for(int si=0; si<num_species; si++) {
            species[si].prim[j].internal_energy += species[si].cv*friction_vec_output(si)*dt;
            species[si].prim[j].temperature     += friction_vec_output(si)*dt;
        }
    }

    //
    // Update conserved quantities
    //
    for(int si=0; si<num_species; si++) {
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }

}


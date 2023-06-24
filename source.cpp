/**
 * source.cpp
 * 
 * This file contains routines pertaining to the gravitation and friction sources.
 * 
 * 
 */

#define EIGEN_RUNTIME_NO_MALLOC


#include <cassert>
#include "aiolos.h"

///////////////////////////////////////////////////////////
//
// Gravitation
//
//
///////////////////////////////////////////////////////////

/**
 * Initialize the gravity array with a non-self gravitating solution and the tidal field.
 * Computes rhill which is output in the beginning of the run.
 */
void c_Sim::init_grav_pot() {

    double rh = planet_semimajor * au * pow(planet_mass / (3.* star_mass ),0.333333333333333333);        

    for(int i = 1; i <= num_cells+1; i++) {
        enclosed_mass[i] = planet_mass;      //No self-gravity
        phi[i]           = get_phi_grav(x_i12[i], enclosed_mass[i]);
            
        if(use_tides == 1) {
        	//if(x_i[i]>0.2*rh)
            		phi[i] -=0.5* 3.*G*init_star_mass*x_i12[i]*x_i12[i]/pow(planet_semimajor*au,3.);
	}
            //phi[i] -=0.5* 3.*G*star_mass*x_i12[i]*x_i12[i]/pow(planet_semimajor*au,3.);
            
    }
    phi[0]           = get_phi_grav(x_i12[1],         planet_mass);
        
    rhill = planet_semimajor*au * std::pow(planet_mass / 3. / star_mass, 0.333333333333333333333);
}


/**
* Compute the integrated mass in shells and the gravity array as a self-gravitating solution
*/
void c_Sim::update_mass_and_pot() {
    
    double species_dens_sum;
    
    if(debug >= 1)
        cout<<"In update mass. "<<endl;
    
    enclosed_mass[0] = planet_mass;
    enclosed_mass[1] = planet_mass; //Ghosts
    
    for(int i = 2; i <= num_cells; i++) {
        
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
            double rh = planet_semimajor * au * pow(planet_mass / (3.* star_mass ),0.333333333333333333);
	    double m0 = init_star_mass;
	    double m1 = star_mass;
            double t0 = ramp_star_mass_t0;
	    double t1 = ramp_star_mass_t1;

	    double a = (m0-m1)/(t0-t1);
	    double b = 0.5*(m0+m1 - a *(t1+t0));

	    double current_star_mass = m0;
	    if(globalTime > t0)
	     	current_star_mass =  a * globalTime + b;
	    if(globalTime > t1)
		current_star_mass = star_mass;

            //if(steps==5)
	   //	cout<<" In tidal grav, t0/t1 = "<<t0<<"/"<<t1<<" m0/current_mass ="<<m0/msolar<<"/"<<current_stellar_mass/msolar<<endl;

            d3 = d3*d3*d3;
	    //if(x_i[i]>0.2*rh)
	            phi[i] -= 0.5* 3.*G*star_mass*x_i12[i]*x_i12[i]/d3;
        }
            
    }
    
    if(debug >= 1)
        cout<<"Done in update mass and pot. "<<endl;
}

//
/**
 * Helper function implementing the actual gravity law used (i.e. $\phi \propto -1/r$ or $\phi \propto +z^2$). Works also in cartesian coordinates
 * Options: linear gravity and 1/r gravity with gravitational smoothing from Klahr&Kley 2006
 * 
 * @param[in] r physical, intercell radii in cm
 * @param[in] mass physical mass inside shells of radii r in g
 * @return   grav potential in erg/g
 */
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

/**
 * Compute the well-balanced source function via the lhs and rhs grav potential differences for any cell-centered index j. After Eq. 11 in KäppeliMishra2016
 * This function uses phi_s, the species-dependent gravitational potential. If the turbulent mixing parameter is 0, then phi_s = phi. If it is non-zero then phi_s regulates
 * how different species feel a weaker field in the lower atmopshere, to adjust their scale height to the dominant species, assumed to be species[0].
 * 
 * @param[in] u Conservative data in cell j
 * @param[in] j Interface number
 * @return    The well-balanced potential difference in cell j
 */
AOS c_Species::source_grav(AOS &u, int &j) {
    assert(j > 0 && j <= num_cells) ;

    double dphidr_p = (phi_s[j] - phi_s[j-1]) / (base->dx[j] + base->dx[j-1]) ;
    double dphidr_m = (phi_s[j+1] - phi_s[j]) / (base->dx[j+1] + base->dx[j]) ;

    return AOS(0, u.u1, u.u2) * (-1.) * ( base->omegaplus[j] * dphidr_p  + base->omegaminus[j] * dphidr_m);
}


/**
 * For scenarios using a changed lower atmospheric gravity field to equalize all scale heights for all species.
 * Switched on using 
 * NOTE: Future applications might require modifying the loop to find homopause.
 * 
 * @param[in] argument Species number. In current use this is identical to this->this_species_index, but the option should be open to point to one dominant species.
 */
void c_Species::update_kzz_and_gravpot(int argument) {
    
    
    
    int homopause_boundary_i = 0;
    double mu = base->species[0].mass_amu;
    double mi = this->mass_amu;
    
    if(argument > 0 && mi > 1e-3) { //Do the kzz buildup for all other species except atomic hydrogen (assumed species[0]) and no electrons (mass = 5e-4).
        
        for(int i=0; i<num_cells+2; i++) {
            double n   = base->species[0].prim[i].number_density;
            double kzz = base->K_zz[i]; //This is a meta-parameter for k_zz/b 
            double par = std::pow(n, 1./1.);
            
            K_zzf[i] = (1. + kzz*par*mu/mi) / (1. + kzz*par);
            if(kzz*n < 1.)
                K_zzf[i] = 1.;
            else
                K_zzf[i] = mu/mi;
            
            double one = 0.99999;
            if(K_zzf[i] > one && K_zzf[i-1] < one) //found homopause
                homopause_boundary_i = i;
        }
            
    }
    else {
        for(int i=0; i<num_cells+2; i++) {
            K_zzf[i] = 1.;
        }
    }
    
    for(int i=0; i<num_cells+2; i++) {
                phi_s[i] = base->phi[i] * K_zzf[i];
    }
    double phicorrection = base->phi[homopause_boundary_i]*(1.-mu/mi); //Correct for the jump in Phi at the homopause
    if(homopause_boundary_i == 0)
        phicorrection = 0;
    
    for(int i=0; i<homopause_boundary_i; i++) 
        phi_s[i] += phicorrection;
    
    //cout<<"s = "<<argument<<" phicorr = "<<phicorrection<<" homopause_i = "<<homopause_boundary_i<<endl;
    //cout<<"Finished updating kzz in species "<<speciesname<<endl;
}



///////////////////////////////////////////////////////////
//
// Friction
//
//
///////////////////////////////////////////////////////////

/**
 * Wrapper function calling the appropriate analytical or numerical functions.
 */
void c_Sim::compute_drag_update() {
	//cout<<"In  drag_update"<<endl;
        if(friction_solver >= 0 && num_species > 1) {
            if(friction_solver == 0)
                compute_friction_analytical();
            else{
		//cout<<"calling friction numerical"<<endl;
                 compute_friction_numerical();
		}
        }
    }
    

/**
 *  Friction solver 0. Only works for two species properly. Three species solution was copied&converted from mathematica matrix inversion, but produces nonsense.
 */
void c_Sim::compute_friction_analytical() {

        if(debug > 0) cout<<"in analytic friction, num_species = "<<num_species<<endl;
        
        if(num_species == 2) {
        
            //Apply analytic solutions
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
                    } else { // Gas-gas collisions
                        coll_b      = 5.0e17 * std::pow(meanT, 0.75) ;     // from Zahnle & Kasting 1986 Tab. 1
                        alpha = kb * meanT * species[1].prim[j].number_density / (m0 * coll_b) ; // From Burgers book, or Schunk & Nagy.
                        alpha *= alpha_collision;
                    }
                }
                
                alphas_sample(j) = alpha;
                friction_sample(j) = alphas_sample(j) * (v2b - v1b);
                
                eps   = species[0].u[j].u1 / species[1].u[j].u1;
                f1    = (1. + dt*alpha)/(1. + dt*alpha*(1.+eps)) ;
                f2    = dt*alpha*eps / (1. + dt * alpha);
                
                v2a    = (v2b + v1b * f2 ) * f1;
                v1a    = v1b - (v2a - v2b) / eps;
                
                species[0].prim[j].speed = v1a;
                species[1].prim[j].speed = v2a;
                
                if(debug >= 1) {
                    //cout<<" v1b, v2b = "<<v1b<<" "<<v2b<<endl<<" v1a, v2a = "<<v1a<<" "<<v2a<<endl;
                    //cout<<" dv1, dv2 = "<<v1a-v1b<<" "<<v2a-v2b<<endl;
                    //cout<<" dp1, dp2 = "<<(v1a-v1b)*species[0].u[j].u1 + (v2a-v2b)*species[1].u[j].u1<<endl;
                    cout<<"    dt*alpha = "<<dt*alpha<<" alpha = "<<alpha<<" eps = "<<species[1].u[j].u1/species[0].u[j].u1<<" dmom/|mom| = "<<((v1a-v1b)*species[0].u[j].u1 + (v2a-v2b)*species[1].u[j].u1)/(std::sqrt(species[0].u[j].u2*species[0].u[j].u2) + std::sqrt(species[1].u[j].u2*species[1].u[j].u2) )<<" meanT = "<<meanT<<" coll_b = "<<coll_b<<" cell = "<<j<<endl;
                }
                
                double v_half = 0.5*(v1a + v1b) - 0.5*(v2a + v2b) ;
                double v_end  = v1a - v2a ;

                species[0].prim[j].internal_energy += dt * alpha * (species[1].mass_amu/(species[0].mass_amu+species[1].mass_amu)) * v_half * v_end ;
                species[1].prim[j].internal_energy += dt * alpha * eps * (species[0].mass_amu/(species[0].mass_amu+species[1].mass_amu)) * v_half * v_end ;
            }
            
            if(debug > 0) {
                double dekin1 = 0.5*species[0].u[num_cells+1].u1 * (v1a - v1b) * (v1a - v1b) ; 
                double dekin2 = 0.5*species[1].u[num_cells+1].u1 * (v2a - v2b) * (v2a - v2b) ;
                double dvrel1 = (v1a - v1b)/v1b;
                double dvrel2 = (v2a - v2b)/v2b;
                cout<<"Relative differences in velocities 1/2 = "<<dvrel1<<" / "<<dvrel2<<" differences in Ekin = "<<dekin1<<" / "<<dekin2<<endl;
                char a;
                cin>>a;
            }
        }
        if(num_species == 3) { //Experimental 3-species solution
                
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
                ///TODO: Internal energy update, if 3-species solver ever works.
                
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

/**
 * Friction solver 1: Compute collision coefficients, build collision matrix and solve in a fast way
 * Collision coefficients are hard-coded in compute_alpha_matrix().
 */
void c_Sim::compute_friction_numerical() {
    
    Eigen::internal::set_is_malloc_allowed(false) ;
         
    if(debug > 0) cout<<"in numerical, dense friction, num_species = "<<num_species<<endl;
    
    for(int j=0; j <= num_cells+1; j++){

        fill_alpha_basis_arrays(j);
        compute_alpha_matrix(j);
                
        if(debug >= 2 && j==1200 && steps == 10)  {
            cout<<"    Computed coefficient matrix ="<<endl<<friction_coefficients<<endl;
            cout<<"    velocities ="<<endl<<friction_vec_input<<endl;
            cout<<"    velocities[0] = "<<friction_vec_input(0)<<endl;
            cout<<"    rho[0] = "<<species[0].u[j].u1<<endl;
        }
        
        friction_matrix_T = identity_matrix - friction_coefficients * dt;
        friction_matrix_T.diagonal().noalias() += dt * (friction_coefficients * unity_vector);
        
        LU.compute(friction_matrix_T) ;
        friction_vec_output.noalias() = LU.solve(friction_vec_input);
        
        if(debug >= 3 && j==700 && steps == 10) {
            
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
      
        //
        // Update new speed and internal energy
        //
        
	if(use_avg_velocity) {
		double tot_mom = 0;
		double tot_dens= 0;
		double avg_velocity = 0;

		for(int si=0; si<num_species; si++) {
			tot_mom += friction_vec_output(si) * species[si].prim[j].density;
			tot_dens += species[si].prim[j].density;
		}

		avg_velocity = tot_mom/tot_dens;

		if(globalTime < avg_velocity_t1) {


			if(globalTime < avg_velocity_t0) {
				for(int si=0; si<num_species; si++) {
					species[si].prim[j].speed = avg_velocity;
				}

				//if(steps%10000==0 && j==num_species/2)
                                //        cout<<" In avg velocity, before t0, avg_v= "<<avg_velocity<<endl;
			}
			else {
				alpha_collision = 1.; //Switch collision alphas to their nominal values
				double a = (1+ 1e-10)/(avg_velocity_t1 - avg_velocity_t0 + 1e-10);
				double b = (1+ 1e-10)/(1-avg_velocity_t1/avg_velocity_t0 + 1e-10);
				double fac = a * globalTime + b;

				//if(steps%50000==0)
				//	cout<<" In avg velocity, between t0 and t1, fac= "<<fac<<endl;

				for(int si=0; si<num_species; si++) {
                                        species[si].prim[j].speed = friction_vec_output(si) * fac + avg_velocity * (1. - fac);
                                }
			}
		}
		else {
			for(int si=0; si<num_species; si++)
		                species[si].prim[j].speed = friction_vec_output(si);
		}


	} else {
        	for(int si=0; si<num_species; si++)
            	species[si].prim[j].speed = friction_vec_output(si);
	}
        
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
    // Update conserved
    //
    for(int si=0; si<num_species; si++) {
        species[si].eos->update_p_from_eint(&(species[si].prim[0]), num_cells+2);
        species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        
    }
    
}

/**
 *  Build all primitives that we need multiple times during the construciton of the alpha matrix
 * 
 * @param[in] j Cell number in which to compute the temp data
 */
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

/**
 * Build collision matrix. Here we also distinguish between gas/dust species and between charged/neutral species in terms of collission coefficients.
 * 
 * * @param[in] j Cell number in which to build the collision matrix
 */
void c_Sim::compute_alpha_matrix(int j) { //Called in compute_friction() and compute_radiation() in source.cpp
        
	if(steps< 4 && j ==100)
		cout<<steps<<" In compute_alpha_matrix j==100"<<endl;
        double alpha_local;
        double coll_b;
        //double mtot;
        double mumass, mumass_amu, meanT;

        for(int si=0; si<num_species; si++) {
            for(int sj=0; sj<num_species; sj++) {
                    
                    if(collision_model == 'C') { //Constant drag law
                        alpha_local = alpha_collision;
                        if (si > sj) alpha_local *= dens_vector(sj) / dens_vector(si) ;
                    } else {
                        if(si==sj) {
                        alpha_local = 0.; 
                        } else {                      // Physical drag laws
                            //mtot = mass_vector(si) + mass_vector(sj);
                            mumass = mass_vector(si) * mass_vector(sj) * inv_totmasses(si,sj); /// (mass_vector(si) + mass_vector(sj));
                            mumass_amu = mumass/amu;
                            meanT  = (mass_vector(sj)*temperature_vector(si) + mass_vector(si)*temperature_vector(sj)) * inv_totmasses(si,sj); // / (mass_vector(si) + mass_vector(sj)); //Mean collisional mu and T from Schunk 1980
                            
                            if (species[si].is_dust_like || species[sj].is_dust_like) { //One of the collision partners is dust
                                // Use Epstein drag law (Hard sphere model)
                                double RHO_DUST = 1;
                                double s0=0, s1=0 ;
                                if (species[si].is_dust_like)
                                    s0 = std::pow(3*mass_vector(si)/(4*M_PI*RHO_DUST), 1/3.) ;
                                if (species[sj].is_dust_like)
                                    s1 = std::pow(3*mass_vector(sj)/(4*M_PI*RHO_DUST), 1/3.) ;
                                
                                double A = M_PI*(s0+s1)*(s0+s1) ;
                                double v_th = std::sqrt(8*kb*meanT/(M_PI * mumass)) ;

                                alpha_local = 4/3. * dens_vector(sj) * v_th * A * inv_totmasses(si,sj) ;
                            } else { // Gas-gas collisions. Used binary diffusion coefficient                        

                                int ci = (int)(species[si].static_charge*1.01);
                                int cj = (int)(species[sj].static_charge*1.01);
                                double qn = (std::fabs(species[si].static_charge) + std::fabs(species[sj].static_charge));
                                //string ccase = ""; //Debug string to make sure the cases are picked right
                                
                                if(std::abs(ci) + std::abs(cj) == 0) { //n-n collision
                                    coll_b      = 5.0e17 * std::sqrt(std::sqrt(meanT*meanT*meanT)) ;     //std::pow(meanT, 0.75) // from Zahnle & Kasting 1986 Tab. 1
                                    alpha_local = kb * meanT * numdens_vector(sj) * species[si].inv_mass / coll_b ; // From Burgers book, or Schunk & Nagy.
                                    //ccase = " n-n ";
                                }
                                else if(std::abs(ci) == 0 || std::abs(cj) == 0) {//i-n collision
                                    
                                    if(species[si].mass_amu < 1.1 && species[sj].mass_amu < 1.1) {  //Resonant H,H+ collisions
                                        alpha_local = 2.65e-10 * numdens_vector(sj) * std::sqrt(meanT) * std::pow((1. - 0.083 * std::log10(meanT)), 2.);
                                        
                                        //ccase = "resonant H-p+ ";
                                        
                                    } else {

                                        alpha_local = 1*2.21 * 3.141592 * numdens_vector(sj) * mass_vector(sj) * inv_totmasses(si,sj); //   /(mass_vector(sj)+mass_vector(si));
                                        alpha_local *= std::sqrt(0.66 / mumass ) * 1e-12 * qn * elm_charge; //0.66 is the polarizability of neutral atomic H                                        
                                        alpha_local *= resonant_pair_matrix(si,sj); //Crude way to emulate resonant collision cross-sections: pairs are multiplied by 10, everything else by 1
                                        
                                        //ccase = "resonant or nonresonant i-n ";
                                    }
                                }
                                else { //i-i collision
                                    alpha_local = alpha_collision_ions * 1.27 * species[si].static_charge * species[si].static_charge * species[sj].static_charge * species[sj].static_charge * std::sqrt(mumass_amu) * amu * species[si].inv_mass; //  / mass_vector(si);
                                    alpha_local *= numdens_vector(sj) / std::sqrt(meanT*meanT*meanT);
                                    
                                    //ccase = " i-i ";
                                }
                                
                                if(si!=sj && debug>3)
                                    cout<<"cell "<<j<<" colliding "<<species[si].speciesname<<" q="<<species[si].static_charge<<" with "<<species[sj].speciesname<<" q="<<species[sj].static_charge<<" ccase = "<<" alpha/n = "<<alpha_local/ numdens_vector(sj)<<" n_sj = "<<numdens_vector(sj)<<endl;
                                
                                double alpha_temp = alpha_collision;
                                if(globalTime < coll_rampup_time) { //Ramp up the friction factors linearly, if so desired
                                    
                                        alpha_temp = alpha_collision * ( init_coll_factor + (1.-init_coll_factor) * globalTime/coll_rampup_time );
                                    
                                } else {
                                    alpha_temp = alpha_collision;
                                }
                                
                                //alpha_local *= alpha_collision;
                                alpha_local *= alpha_temp;
                            }
                        }
                    }
                    
                    if(si==0 && sj == 1) {
                        alphas_sample(j) = alpha_local; //alphas_sample is saving alpha values for later use outside of the friction function, e.g. to output it
                        //cout<<"    spec "<<species[si].name<<" j = "<<j<<" alpha_local = "<<alpha_local<<endl;
                        if(debug > 0 && steps == 1 && ((j==2) || (j==num_cells-2) || (j==num_cells/2)) )
                            cout<<"    spec "<<species[si].speciesname<<" j = "<<j<<" alpha_local = "<<alpha_local<<endl;
                    }
                    
                    // Fill alpha_ij and alpha_ji
                    //if(si==sj)
                   //    friction_coefficients(si,sj) = 0.;
                    //else
                    friction_coefficients(si,sj) = friction_coeff_mask(si,sj) * alpha_local;
                    //friction_coefficients(sj,si) = friction_coeff_mask(sj,si) * alpha_local * dens_vector(si) / dens_vector(sj);
                    //alpha_local *=  ;
                
            }
        }
        //char a;
        //cin>>a;
}


/**
 * Build collisional heat exchange matrix. Uses pre-computed values from momentum transfer.
 * Called from compute_collisional_heat_exchange.
 * 21 May 2023: This function has been disabled in radiation_simple, the coefficients are computed in place there now
 * 
 * * @param[in] j Cell number in which to build the heat exchange matrix
 */

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


/**
 * Solve for collisional heat exchange
 */
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


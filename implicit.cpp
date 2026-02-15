/**
 * implicit.cpp
 * 
 * This file contains routines solving the hydrodynamic euler system with source terms in an implicit fashion
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"

inline double minmod(double a, double b) {
    if (a*b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Incompressible
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * Computes the implicit hydro solution for incompressible fluids - applicable to gases e.g. at low mach numbers for electrons
 * 
 * @param sp: the species of which to compute the timestep advance
 */
void c_Species::implicit_incompressible(double dt) {

    //cout<<" Hi i am species "<<speciesname<<" and i am attempting to be solved implicitly."<<endl;

    //Boundaries
    this->apply_boundary_left(this->u) ;
    this->apply_boundary_right(this->u) ;
    compute_pressure(this->u);
    std::vector<double> u_mask = np_zeros(num_cells+2); 
    reconstruct_edge_states(u_mask, 2) ;

    std::vector<double> t_damp = np_zeros(num_cells+2); 
    std::vector<double> ff     = np_zeros(num_cells+2); 
    /*
    Matrix_t adv_mat      = Matrix_t::Zero(num_cells+2, num_cells+2);
    Matrix_t adv_id       = Matrix_t::Identity(num_cells+2, num_cells+2);
    Vector_t adv_b        = Vector_t(num_cells+2);
    Vector_t results      = Vector_t(num_cells+2); 
    */

    //Allocate matrices
    int num_vars = 3; // + num_species
    int stride = num_vars * num_vars ;
    int size_r = (base->num_cells + 2) * num_vars ;
    int size_M = (base->num_cells + 2) * stride ;
    const int wall = 2;
    const int pwall= 2;
    const int ewall= 2;

    std::vector<double> 
        ll(size_M, 0.), dd(size_M, 0.), uu(size_M, 0.), r(size_r, 0.) ;
    
    ////////////////////////////////////////////////////////////////////////
    // precompute MUSCL slopes
    ////////////////////////////////////////////////////////////////////////
    std::vector<double> slope_rho(num_cells+2, 0.0);
    std::vector<double> slope_mom(num_cells+2, 0.0);
    std::vector<double> slope_E(num_cells+2, 0.0);
    std::vector<double> slope_p(num_cells+2, 0.0);

    const std::vector<double>& 
            x_i = base->x_i, 
            x_iVC = base->x_iVC,
            dx = base->dx ;
    
    
    for (int j = 2; j < num_cells-1; j++) {

        double cF = (x_iVC[j+1] - x_iVC[j]) / (x_i[j] - x_iVC[j]) ;
        double cB = (x_iVC[j] - x_iVC[j-1]) / (x_iVC[j] - x_i[j-1]) ;

        double dxF = (x_iVC[j+1] - x_iVC[j]) ;
        double dxB = (x_iVC[j] - x_iVC[j-1]) ;

        slope_rho[j] = 1. * reconstruct_pointer(u[j-1].u1, u[j].u1, u[j+1].u1, cF, cB, dxF, dxB) ;
        slope_mom[j] = 1. * reconstruct_pointer(u[j-1].u2, u[j].u2, u[j+1].u2, cF, cB, dxF, dxB) ;
        slope_E[j]   = 1. * reconstruct_pointer(u[j-1].u3, u[j].u3, u[j+1].u3, cF, cB, dxF, dxB) ;
        slope_p[j]   = 1. * reconstruct_pointer(prim[j-1].pres, prim[j].pres, prim[j+1].pres, cF, cB, dxF, dxB) ;

        if( std::isnan(slope_rho[j]) || std::isnan(slope_mom[j]) || std::isnan(slope_E[j]) || std::isnan(slope_p[j])) {
            slope_rho[j] = 0;
            slope_mom[j] = 0;
            slope_E[j] = 0;
            slope_p[j] = 0;
        }
    }

   
    /* if(base->steps == 49) {
        for (int j = 1; j < num_cells-1; j++) {
            cout<<" j / slope/ drho/ rho "<<j<<" "<<slope_rho[j]<<" "<<slope_rho[j]*(x_i[ j ] - x_iVC[j])<<" "<<prim[j].density<<endl;
        }
    }*/
    
    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////
    for (int j=0; j <= num_cells+1; j++) {

        double V   = base->vol[j];
        //cout<<" writing matrix j "<<j<<endl;
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Time dependent terms:
        //rho
        dd[idx]      += V / dt ;
        r[idx_r]     += V / dt * u[j].u1;

        //momentum
        dd[idx + 4]  += V / dt ;
        r[idx_r+ 1]  += V / dt * u[j].u2;

        //energy
        dd[idx + 8]  += V / dt ;
        r[idx_r+ 2]  += V / dt * u[j].u3;
    }
    //cout<<" steps "<<base->steps<<" v_l's = ";

    for (int j=1; j < num_cells+1; j++) {

        double V   = base->vol[j];
        double S_l = base->surf[j-1];
        double S_r = base->surf[j];
        //cout<<" writing matrix j "<<j<<endl;
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Face velocities (from non-advanced timestep)
        //double v_l = 0.5 * (prim[j-1].speed + prim[j].speed); //(u[j-1].u2 / u[j-1].u1 + u[j].u2 / u[j].u1);
        //double v_r = 0.5 * (prim[j].speed + prim[j+1].speed); // ;(u[j].u2 / u[j].u1     + u[j+1].u2 / u[j+1].u1);
        
        //
        //
        double v_l   = (std::sqrt(u[j-1].u1) * prim[j-1].speed + std::sqrt(u[j].u1) * prim[j].speed )/( std::sqrt(u[j-1].u1)  + std::sqrt(u[j].u1) );
        double v_r   = (std::sqrt(u[j+1].u1) * prim[j+1].speed + std::sqrt(u[j].u1) * prim[j].speed )/( std::sqrt(u[j+1].u1)  + std::sqrt(u[j].u1) );


        //
        // Momentum damping at low abundances to avoid numerical noise
        //
        double f      =  prim[j].pres/base->total_press[j];
        double f_lim  = base->edamp_lim;
        ff[j] = f;

        t_damp[j] = std::max(f*f_lim, 1e-10);
        if(j>= num_cells-1)
            t_damp[j] = 1e-2;
        
        dd[idx + 4]      += V / t_damp[j] ;

        //Boundaries 1
        if(0==1) {
            if(j<wall) {
                v_l = 0;
                v_r = 0;
            }
            if(j==wall) {
                v_l = 0.;
            }
        }
        if(j==num_cells-1)
            v_r = 0;
        if(j>num_cells-1) {
            v_l = v_r = 0;
        }
        
        double gsrc    = source_grav(u[j], j).u2;
        double gsrc3   = source_grav(u[j], j).u3;
        double gsrc_no = source_grav_noconserved(u[j], j).u2;
        double psrc    =  -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
                                          base->source_pressure_prefactor_right[j] * prim_r[j].pres);
        
        //if(j<6)
        //    cout<<v_l<<" ";
        //if(base->steps==430)                                          
        //    cout<<" j = "<<j<<" pl  = "<<prim_l[j].pres<<" pr "<<prim_r[j].pres<<" pc "<<prim[j].pres<<" src "<<psrc<<endl;                               
        //General left and right fluxes, flux = lam * density
        AOS lam_l = AOS(v_l * S_l, 0., 0.);
        AOS lam_r = AOS(v_r * S_r, 0., 0.);
        AOS lam_l_ex;
        AOS lam_r_ex;

        double theta = base->implicit_theta; //Regulates the balance between implicit and explicit terms in Crank-Nicolson scheme
        double drho  = 0, dmom = 0, dp = 0, dE=0, a=0, b=0, ap=0, bp=0;
        double tmp_r = 0, tmp_l=0;
        ////////////////////////////////////////////////////////////////////////
        // Right face
        ////////////////////////////////////////////////////////////////////////
        int sw=v_r>0? 0 : 1;
        int sgn = v_r>0 ? +1 : -1;
        ///////////////////////////////////
        //Jacobian entries
        double dg_dm   = 2*v_r;
        double dg_drho = -v_r*v_r;

        ///////////////////////////////////
        if(j<num_cells) {
        drho = -sgn * slope_rho[j+sw]  / (u[j+1].u1  - u[j].u1  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
        dmom = -sgn * slope_mom[j+sw]  / (u[j+1].u2  - u[j].u2  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); 
        dE   = -sgn * slope_E[j+sw]    / (u[j+1].u3  - u[j].u3  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); 
        dp   = -sgn * slope_p[j+sw]    / (prim[j+1].pres  - prim[j].pres  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); 
        }
        
        //rho
        a = v_r<0? 1.0 + 0.5  * drho : -0.5 * drho;
        b = v_r<0? -0.5 * drho : 1.0 + 0.5  * drho;
        //dd[idx + 0]    += theta     * lam_r.u1 * b;
        //uu[idx + 0]    += theta     * lam_r.u1 * a;

        //r[idx_r]       -= (1-theta) * lam_r.u1 * ( a * u[j+1].u1 + b * u[j].u1 );

        //momentum,
        //       /*
        a = v_r<0? 1.0 + 0.5  * dmom : -0.5 * dmom;
        b = v_r<0? -0.5 * dmom : 1.0 + 0.5  * dmom;

        dd[idx + 1]    +=       S_r * b;
        uu[idx + 1]    +=       S_r * a;

        dd[idx + 4]    += theta     * v_r * S_r * b; //mom_j+1/2
        uu[idx + 4]    += theta     * v_r * S_r * a;
        
        r[idx_r+1]       -= (1-theta) * v_r * S_r * ( a * u[j+1].u2 + b * u[j].u2);


        a = v_r<0? 1.0 + 0.5  * dp : -0.5 * dp;
        b = v_r<0? -0.5 * dp : 1.0 + 0.5  * dp;


        if(j>pwall && j<num_cells) {
            r[idx_r+1]       -=  0.5*S_r * (prim[j+1].pres + prim[j].pres )   ; //momentum,
        }
        //Jacobian entries
                                                     
        //        */
        //Energy
        a  = v_r<0? 1.0 + 0.5  * dE : -0.5 * dE;
        b  = v_r<0? -0.5 * dE : 1.0 + 0.5  * dE;
        
        dd[idx + 8]    += theta     * v_r * S_r * b; //E_j+1/2
        uu[idx + 8]    += theta     * v_r * S_r * a;
        r[idx_r+ 2]    -= (1-theta) * v_r * S_r * ( a * u[j+1].u3        + b * u[j].u3);
        if(j>=ewall && j<num_cells) {
            r[idx_r+ 2]    +=   v_r * S_r * ( a* prim_l[j+1].pres + b*prim_r[j].pres ); //p div v
        }
 
        /////////////////////////////////////////////////////////////
        // New terms
        /* a=b=0.5;
        std::vector<double> dPdu  = get_hydro_jacobian_P(j, j+1, a, b);
        std::vector<double> avgus = {b * u[j+1].u1 + a * u[j].u1, b * u[j+1].u2 + a * u[j].u2, b * u[j+1].u3 + a * u[j].u3};
        
        dd[idx + 8]    +=          v_r * S_r * b; //E_j+1/2
        uu[idx + 8]    +=          v_r * S_r * a;
        dd[idx + 8]    +=          v_r * S_r * b * dPdu[2]; //v dP_dE * E
        uu[idx + 8]    +=          v_r * S_r * a * dPdu[2];
 
        dd[idx + 7]    +=          v_r * S_r * b * dPdu[1]; //v dP_dmom * mom
        uu[idx + 7]    +=          v_r * S_r * a * dPdu[1];
        dd[idx + 6]    +=          v_r * S_r * b * dPdu[0]; //v dP_drho * rhp
        uu[idx + 6]    +=          v_r * S_r * a * dPdu[0];
 
        //r[idx_r+ 2]    -= (1-theta) * v_r * S_r * ( a * u[j+1].u3        + b * u[j].u3);
        if(j>=ewall && j<num_cells) {
            r[idx_r+ 2]    -=           v_r * S_r * (a* prim_l[j+1].pres + b*prim_r[j].pres ); //p div v
            r[idx_r+ 2]    +=           v_r * S_r * (dPdu[0]*avgus[0] + dPdu[1]*avgus[1] + dPdu[2]*avgus[2]); //p div v
        }
  */
        // End new terms
        /////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////
        // Left face
        //////////////////////////////////////////////////////////////////////
        sw  = v_l>0? -1 : 0;
        sgn = v_l>0 ? +1 : -1;
        ///////////////////////////////////

        drho = -sgn * slope_rho[j+sw]  / (u[j].u1  - u[j-1].u1  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
        dmom = -sgn * slope_mom[j+sw]  / (u[j].u2  - u[j-1].u2  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]);
        dE   = -sgn * slope_E[j+sw]    / (u[j].u3  - u[j-1].u3  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]);
        dp   = -sgn * slope_p[j+sw]    / (prim[j].pres  - prim[j-1].pres  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]);

        //rho
        a = v_l<0? 1.0 + 0.5  * drho : -0.5 * drho;
        b = v_l<0? -0.5 * drho : 1.0 + 0.5  * drho;
        //dd[idx + 0]    -= theta     * lam_l.u1 * a;
        //ll[idx + 0]    -= theta     * lam_l.u1 * b;

        //r[idx_r]       += (1-theta) * lam_l.u1 * ( a * u[j].u1 + b * u[j-1].u1 );
        
        //momentum,
        //       
        a = v_l<0? 1.0 + 0.5  * dmom : -0.5 * dmom;
        b = v_l<0? -0.5 * dmom : 1.0 + 0.5  * dmom;

        dd[idx + 1]    -=  S_l * a;
        ll[idx + 1]    -=  S_l * b;
       
        dd[idx + 4]    -= theta     * v_l * S_l * a; //mom_j-1/2
        ll[idx + 4]    -= theta     * v_l * S_l * b;
        r[idx_r+1]     += (1-theta) * v_l * S_l * ( a * u[j].u2 + b * u[j-1].u2) ;

        a = v_l<0? 1.0 + 0.5  * dp : -0.5 * dp;
        b = v_l<0? -0.5 * dp : 1.0 + 0.5  * dp;

        if(j>pwall && j<=num_cells) {
                r[idx_r+1]     += 0.5*S_l * ( prim[j].pres + prim[j-1].pres ) ; 
        }
        // 
        //energy
         
        a  = v_l<0? 1.0 + 0.5  * dE : -0.5 * dE;
        b  = v_l<0? -0.5 * dE : 1.0 + 0.5  * dE;                                                
        dd[idx + 8]    -= theta     * v_l * S_l * a; //E_j-1/2
        ll[idx + 8]    -= theta     * v_l * S_l * b;
        r[idx_r+2]     += (1-theta) * v_l * S_l * ( a * u[j].u3        + b * u[j-1].u3);         //u grad E
        if(j>ewall && j<=num_cells) {
            r[idx_r+2]     -=   v_l * S_l * ( a *prim_l[j].pres + b*prim_r[j-1].pres ); //p div v
        }

 
        /////////////////////////////////////////////////////////////
        // New terms
        /* a=b=0.5;
        dPdu  = get_hydro_jacobian_P(j, j-1, a, b);
        avgus = {b * u[j-1].u1 + a * u[j].u1, b * u[j-1].u2 + a * u[j].u2, b * u[j-1].u3 + a * u[j].u3};
        
        dd[idx + 8]    -= theta     * v_l * S_l * a; //E_j-1/2
        ll[idx + 8]    -= theta     * v_l * S_l * b;
        dd[idx + 8]    -= theta     * v_l * S_l * a * dPdu[2]; //v dP_dE * E
        ll[idx + 8]    -= theta     * v_l * S_l * b * dPdu[2];
 
        dd[idx + 7]    -= theta     * v_l * S_l * a * dPdu[1]; //v dP_dmom * mom
        ll[idx + 7]    -= theta     * v_l * S_l * b * dPdu[1];
        dd[idx + 6]    -= theta     * v_l * S_l * a * dPdu[0]; //v dP_drho * rho
        ll[idx + 6]    -= theta     * v_l * S_l * b * dPdu[0];
 
        r[idx_r+2]     += (1-theta) * v_l * S_l * ( a * u[j].u3        + b * u[j-1].u3);         //u grad E
        if(j>ewall && j<=num_cells) {
            r[idx_r+2]     +=           v_l * S_l * (a* prim_l[j].pres + b * prim_r[j-1].pres ); //p div v
            r[idx_r+2]     -=           v_l * S_l * (dPdu[0]*avgus[0] + dPdu[1]*avgus[1] + dPdu[2]*avgus[2]); //p div v
        }
  */
        // End new terms
        ///////////////////////////////////////////////////////////////
        
        //////////////////////////////////////////////////////////////////////
        // Sources
        //////////////////////////////////////////////////////////////////////
        double srccoeff_l = 0, srccoeff_c = 1, srccoeff_r = 0; //TODO: Determine those so that they reproduce the internal slope for the P source term

        //r[idx_r+1]       -= (1-theta) * V * ( gsrc + 2 / x_iVC[j] * prim[j].pres);   //Grav and geometric source
        if(j>pwall && j<num_cells) {
            r[idx_r+1]       += 1.0 * V * (gsrc + 1. * psrc);   //Grav and geometric source
            r[idx_r+2]       += 1.0 * V * (gsrc3);              //Grav and geometric source
        }
        
        //////////////////////////////////////////////////////////////////////
        // End Sources
        //////////////////////////////////////////////////////////////////////
    }
    //cout<<endl;
    
    //cout<<"Solving"<<endl;
    //
    // Solve!
    //

    base->implicit_tridiag.factor_matrix(&ll[0], &dd[0], &uu[0]) ;
    base->implicit_tridiag.solve(&r[0], &r[0]) ; // Solve in place

    //
    // End Solve
    //

    //
    // Write solution back into variables
    //
    if(base->steps >= 462e99) {
        cout<<" steps "<<base->steps;
    }
    for (int j=0; j < base->num_cells+1; j++) {
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        double rhotmp = std::max(r[idx_r + 0], 1e-60 );
        double momtmp = r[j*num_vars + 1];
        double vimpl_old = u[j].u2/u[j].u1;
        double vimpl_new = momtmp/rhotmp;
        double T_old = prim[j].temperature = std::min(std::max(prim[j].temperature, base->temperature_floor), base->max_temperature);
        double E_implied = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * T_old;      //This is the old temperature
        double E_floor   = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * base->temperature_floor;
        double T_pred = (u[j].u3 - 0.5*momtmp*momtmp/rhotmp ) / (rhotmp * cv);
        if(base->steps == -2 || base->steps == -2700) {
            cout<<" j = "<<j<<" old/new rho  = "<<u[j].u1<<" / "<<rhotmp<<" t_damp = "<<t_damp[j]<<" corresponding f = "<<ff[j]<<endl;
            cout<<" j = "<<j<<" old/new mom  = "<<u[j].u2<<" / "<<r[j*num_vars + 1]<<" v_implied old/new = "<<vimpl_old<<"/"<<vimpl_new<<endl;
            cout<<" j = "<<j<<" old/new E    = "<<u[j].u3<<" / "<<r[j*num_vars + 2]<<endl;
            cout<<" j = "<<j<<" predicted T  = "<<T_pred<<" current T "<<prim[j].temperature<<" u1/u2/u3: "<<rhotmp<<" / "<<momtmp<<" / "<<r[j*num_vars + 2]<<" E_implied "<<E_implied<<endl;
        }

        if(base->steps>2) {
            double rhonew = std::max(r[idx_r + 0], 1e-60 );//std::max(results(j), 1e-20 ); //std::max(r[idx_r + 0], 1e-20 );//r[idx_r + 0];
            double momnew = r[idx_r + 1];
            double Enew   = r[j*num_vars + 2];
            
            if(j>num_cells-1) { 
                //cout<<" boundary in cell "<<j<<" E_implied "<<E_implied<<" T_pred "<<T_pred<<" T_previous "<<prim[j].temperature<<endl; 
                Enew = E_implied;
            }

            if( (Enew < 0.) || (T_pred < base->temperature_floor) ) {
                //cout<<" negative enew in cell "<<j<<" Enew "<<Enew<<" T_pred "<<T_pred<<endl; 
                Enew = E_floor;
            }
            
                

            u[j].u1 = rhonew;//r[j*num_vars + 0];
            u[j].u2 = momnew;
            //u[j].u3 = 0.5*momnew*momnew/rhonew + rhonew * enew;//r[j*num_vars + 2] ;
            u[j].u3 = Enew ;
            //u[j].u3 = E_implied ; //This amounts to the approximation of local isothermality over the timestep. can lead to -p otherwise.
        }
    }

    // Update primitives. 
    eos->compute_primitive(&(u[0]), &(prim[0]), base->num_cells+2) ;   
    eos->compute_auxillary(&(prim[0]), base->num_cells+2);

}

/**
 * Computes the implicit hydro solution for incompressible fluids - applicable to gases e.g. at low mach numbers for electrons
 * This "version 2" acts on the internal energy, whereas the previous, above solver solves the total energy
 * 
 * @param sp: the species of which to compute the timestep advance
 */
void c_Species::implicit_incompressible2(double dt) {

    //cout<<" Hi i am species "<<speciesname<<" and i am attempting to be solved implicitly."<<endl;

    //Boundaries
    this->apply_boundary_left(this->u) ;
    this->apply_boundary_right(this->u) ;
    compute_pressure(this->u);
    std::vector<double> u_mask = np_zeros(num_cells+2); 
    reconstruct_edge_states(u_mask, 2) ;

    std::vector<double> t_damp = np_zeros(num_cells+2); 
    std::vector<double> ff     = np_zeros(num_cells+2); 
    /*
    Matrix_t adv_mat      = Matrix_t::Zero(num_cells+2, num_cells+2);
    Matrix_t adv_id       = Matrix_t::Identity(num_cells+2, num_cells+2);
    Vector_t adv_b        = Vector_t(num_cells+2);
    Vector_t results      = Vector_t(num_cells+2); 
    */

    //Allocate matrices
    int num_vars = 3; // + num_species
    int stride = num_vars * num_vars ;
    int size_r = (base->num_cells + 2) * num_vars ;
    int size_M = (base->num_cells + 2) * stride ;
    const int wall = 2;
    const int pwall= 2;
    const int ewall= 2;

    std::vector<double> 
        ll(size_M, 0.), dd(size_M, 0.), uu(size_M, 0.), r(size_r, 0.) ;
    
    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////
    for (int j=0; j <= num_cells+1; j++) {

        double V   = base->vol[j];
        //cout<<" writing matrix j "<<j<<endl;
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Time dependent terms:
        //rho
        dd[idx]      += V / dt ;
        r[idx_r]     += V / dt * u[j].u1;

        //momentum
        dd[idx + 4]  += V / dt ;
        r[idx_r+ 1]  += V / dt * u[j].u2;

        //energy
        dd[idx + 8]  += V / dt ;
        r[idx_r+ 2]  += V / dt * u[j].u3;
    }

    //////////////////////////////////////////////////////////////////////////
    // Precompute Jabobians interface by interface
    //////////////////////////////////////////////////////////////////////////

    for (int j=1; j <= num_cells; j++) { //Going interface by interface - interface j is between cell j and j+1. Ignore interface 0, that's just excess memory to get the index shifted by 1
        //write_roe_jacobians(base->roe_differentials_left[j], base->roe_differentials_right[j], j);  //interface j sits at position j-1 in memory
        write_hlle_jacobians(base->roe_differentials_left[j], base->roe_differentials_right[j], j);
    }

    int impl_debug = 0;

    //////////////////////////////////////////////////////////////////////////
    // 
    //////////////////////////////////////////////////////////////////////////

    double geometry_f =0.;
    if(base->geometry==Geometry::spherical) {
        geometry_f = 2;
    } else if(base->geometry!=Geometry::cartesian) {
        geometry_f = 1;
    }

    //cout<<" starting timestep "<<base->steps<<" with dt = "<<dt<<endl;

    for (int j=1; j <= num_cells; j++) {  //Reminder: cell 1 is boundaed by interface 0 (left) and 1 (right), cell 2 by 1 and 2

        double V   = base->vol[j];
        double S_l = base->surf[j-1];
        double S_r = base->surf[j];
        //cout<<" writing matrix j "<<j<<endl;
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Face velocities (from non-advanced timestep)
        //
        // Momentum damping at low abundances to avoid numerical noise
        //
        double f      =  prim[j].pres/base->total_press[j];
        double f_lim  = base->edamp_lim;
        ff[j] = f;

        t_damp[j] = std::max(f*f_lim, 1e-10);
        if(j>= num_cells-1)
            t_damp[j] = 1e-2;
        
        //Boundaries 1
        if(0==1) {
            dd[idx + 4]      += V / t_damp[j] ; //Momentum damping
        }
        
        double gsrc    = source_grav(u[j], j).u2;
        double gsrc3   = source_grav(u[j], j).u3;
        double gsrc_no = source_grav_noconserved(u[j], j).u2;
        double psrc    =  -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
                                          base->source_pressure_prefactor_right[j] * prim_r[j].pres);
        
        if(impl_debug)                                
            cout<<"starting cell j = "<<j<<" steps "<<base->steps<<endl;
        ////////////////////////////////////////////////////////////////////////
        // Right face
        ////////////////////////////////////////////////////////////////////////
        
        Vector3d ul             = u_to_vec(j);
        Vector3d ur             = u_to_vec(j+1);
        Matrix3d Jl = base->roe_differentials_left[j];
        Matrix3d Jr = base->roe_differentials_right[j];
        Vector3d dul            = Jl * ul;
        Vector3d dur            = Jr * ur;
        Vector3d f_last         = get_hlle_flux(j);//roe_flux_vec(j);

        if(impl_debug)                                {
            cout<<" in implicit roe solver Jl and Jr ="<<endl<<Jl<<endl<<endl<<Jr<<endl;
            cout<<" left and right states: "<<endl<<ul<<endl<<ur<<endl;
        }
        //Write R.H.S components
        for(int ii=0; ii<3; ii++) {
            r[idx_r + ii] +=  S_r * (- f_last(ii) + dul(ii) + dur(ii));

            if(impl_debug)                                
                cout<<" was writing into ii = "<<ii<<" r elms, f_roe = "<<f_last(ii) <<" dul+dur = "<<dul(ii) + dur(ii)<<endl;

            for(int jj=0; jj<3; jj++) {        
                int ci = ii*3 + jj;

                dd[idx + ci] += S_r * Jl(ii,jj);
                uu[idx + ci] += S_r * Jr(ii,jj);

                if(impl_debug)                                
                    cout<<" was writing into ci = "<<ci<<" matrix elms, dd = "<<Jl(ii,jj)<<" uu = "<<Jr(ii,jj)<<endl;
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Left face
        //////////////////////////////////////////////////////////////////////
        
        ul             = u_to_vec(j-1);
        ur             = u_to_vec(j);
        Jl = base->roe_differentials_left[j-1];
        Jr = base->roe_differentials_right[j-1];
        dul            = Jl * ul;
        dur            = Jr * ur;
        f_last         = get_hlle_flux(j-1);//roe_flux_vec(j-1);

        //Write L.H.S components
        for(int ii=0; ii<3; ii++) {
            r[idx_r + ii] -=  S_l * (- f_last(ii) + dul(ii) + dur(ii));

            for(int jj=0; jj<3; jj++) {        
                int ci = ii*3 + jj;

                dd[idx + ci] -= S_l * Jr(ii,jj);
                ll[idx + ci] -= S_l * Jl(ii,jj);
                
            }
        }
        
        //////////////////////////////////////////////////////////////////////
        // Sources
        //////////////////////////////////////////////////////////////////////

        //Source terms //2P/r    2 * S_l * base->x_i[j];
        if(base->geometry!=Geometry::cartesian) {
            std::vector<double> diff = get_hydro_jacobian_P(j);
            double prsc_total = (diff[0]*u[j].u1 + diff[1]*u[j].u2 + diff[2]*u[j].u3);

            double rm1 = (S_r-S_l)/V;
                    //rm1 = 2./base->x_iVC[j];
            
            double v2r = + 1.0 * V * rm1;
            dd[idx + 3]    -= v2r * diff[0]; 
            dd[idx + 4]    -= v2r * diff[1];  
            dd[idx + 5]    -= v2r * diff[2]; 
            //r[idx_r+ 1]    += v2r * (prim[j].pres - prsc_total) ;
            //r[idx_r+ 1]    += v2r * ( - prsc_total) + V * psrc;
            r[idx_r+ 1]    += v2r * (prim[j].pres - prsc_total);
        }
        if(j>pwall && j<num_cells) {
            r[idx_r+1]       += 1.0 * V * (gsrc );   //Grav and geometric source
            r[idx_r+2]       += 1.0 * V * (gsrc3);              //Grav and geometric source
        }
        ////////////////////////////////////////////psrc//////////////////////////
        // End Sources
        //////////////////////////////////////////////////////////////////////
    }
    
    //cout<<"Solving"<<endl;
    //
    // Solve!
    //
/* 
    cout<<"r: ";
    for (auto i : r)
        cout<<i<<endl;
    cout<<"dd: ";
    for (auto i : dd)
        cout<<i<<endl;
    cout<<"uu: ";
    for (auto i : uu)
        cout<<i<<endl;
    cout<<"ll: ";
    for (auto i : ll)
        cout<<i<<endl;

    char a;
    cin>>a;
 */
    //cout<<"before solve"<<endl;

    base->implicit_tridiag.factor_matrix(&ll[0], &dd[0], &uu[0]) ;
    base->implicit_tridiag.solve(&r[0], &r[0]) ; // Solve in place

    //cout<<"after solve"<<endl;
    //
    // End Solve
    //

    //
    // Write solution back into variables
    //
    if(base->steps >= 462e99) {
        cout<<" steps "<<base->steps;
    }
    for (int j=0; j <= base->num_cells; j++) {
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        double rhotmp = std::max(r[idx_r + 0], 1e-60 );
        double momtmp = r[j*num_vars + 1];
        double Etmp   = r[j*num_vars + 2];
        double Ekin = 0.5*momtmp*momtmp/rhotmp;

        double T_tmp = prim[j].temperature = std::min(std::max(prim[j].temperature, base->temperature_floor), base->max_temperature);
        
        double rhoe_implied = Etmp-Ekin;
        double E_implied = Ekin + rhotmp * cv * T_tmp;      //This is the old temperature
        double E_floor   = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * base->temperature_floor;
        
        if(base->steps == -2 || base->steps == -2700) {
            //cout<<" j = "<<j<<" old/new rho  = "<<u[j].u1<<" / "<<rhotmp<<" t_damp = "<<t_damp[j]<<" corresponding f = "<<ff[j]<<endl;
            //cout<<" j = "<<j<<" old/new mom  = "<<u[j].u2<<" / "<<r[j*num_vars + 1];
            //cout<<" j = "<<j<<" old/new E    = "<<u[j].u3<<" / "<<r[j*num_vars + 2]<<endl;
            cout<<" j = "<<j<<" predicted T  = "<<rhoe_implied/cv<<" current T "<<prim[j].temperature<<" u1/u2/u3: "<<rhotmp<<" / "<<momtmp<<" / "<<r[j*num_vars + 2]<<" E_implied "<<E_implied<<endl;
        }

        if(base->steps>2) {
            double rhonew = std::max(r[idx_r + 0], 1e-60 );//std::max(results(j), 1e-20 ); //std::max(r[idx_r + 0], 1e-20 );//r[idx_r + 0];
            double momnew = r[idx_r + 1];
            double Enew   = Etmp; //std::min(r[j*num_vars + 2], E_floor ); // we solved for rho e, not e
            
            if( rhoe_implied < 0.) {
                //cout<<" negative enew in cell "<<j<<" "<<endl; 
               Enew = E_floor;
            }

            u[j].u1 = rhonew;
            u[j].u2 = momnew;
            u[j].u3 = Enew;
        }
    }
    // Update primitives. 
    eos->compute_primitive(&(u[0]), &(prim[0]), base->num_cells+1) ;   
    eos->compute_auxillary(&(prim[0]), base->num_cells+1);

}




//std::vector<double> c_Species::get_hydro_jacobian(int jleft, int jright)

//
// Computes the derivatives of the E flux (i.e. third flux component) w.r.t the other conservative hydro variables at an interface
//
std::vector<double> c_Species::get_hydro_jacobian_df3(int jleft, int jright) {
    double g = gamma_adiabat;
    double m1 = 0.5*(u[jleft].u1+u[jright].u1);
    double m2 = 0.5*(u[jleft].u2+u[jright].u2);
    double m3 = 0.5*(u[jleft].u3+u[jright].u3);
    double mv = m2/m1;

    return  {-g*m2*m3/m1+(g-1)*mv*mv*mv,  g*m3/m1-1.5*(g-1)*mv*mv,  g*m2/m1};
}


//
// Computes the derivatives of E (i.e. third hydro vector component) w.r.t the other conservative hydro variables at an interface
//
std::vector<double> c_Species::get_hydro_jacobian_E(int jleft, int jright) {
    double g = gamma_adiabat;
    double m1 = 0.5*(u[jleft].u1+u[jright].u1);
    double m2 = 0.5*(u[jleft].u2+u[jright].u2);
    double m3 = 0.5*(u[jleft].u3+u[jright].u3);
    double mv = m2/m1;

    return  {1, 1,  1};
}

//
// Computes the derivatives of P w.r.t the other conservative hydro variables at an interface
//

std::vector<double> c_Species::get_hydro_jacobian_P(int j) {
    double g = gamma_adiabat;
    double m1 = u[j].u1;
    double m2 = u[j].u2;
    double m3 = u[j].u3;
    double mv = m2/m1;

    return  {(g-1)/2*(mv*mv), - (g-1)*mv,  (g-1)};
}

std::vector<double> c_Species::get_hydro_jacobian_P(int jleft, int jright) {
    double g = gamma_adiabat;
    double m1 = 0.5*(u[jleft].u1+u[jright].u1);
    double m2 = 0.5*(u[jleft].u2+u[jright].u2);
    double m3 = 0.5*(u[jleft].u3+u[jright].u3);
    double mv = m2/m1;

    return  {(g-1)/2*(mv*mv), - (g-1)*mv,  (g-1)};
}

std::vector<double> c_Species::get_hydro_jacobian_P(int jleft, int jright, double a, double b) {
    double g = gamma_adiabat;
    double m1 = 0.5*(a*u[jleft].u1+b*u[jright].u1);
    double m2 = 0.5*(a*u[jleft].u2+b*u[jright].u2);
    double m3 = 0.5*(a*u[jleft].u3+b*u[jright].u3);
    double mv = m2/m1;

    return  {(g-1)/2*(mv*mv), -(g-1)*mv, (g-1)};
}

Matrix3d c_Species::get_exact_Jacobian(AOS u) {
    double g = gamma_adiabat;
    double v = u.u2/u.u1;

    Matrix_t m(3,3);
    //m<< 0,1,0,   0, 1*v,  0.,  \
                0.,  0.,   1*v;

    m<< 0,1,0,   -0.5*(g-3)*v*v, (3-g)*v,  (g-1),  \
    -g*u.u2*u.u3/(u.u1*u.u1)+(g-1)*v*v*v,  g*u.u3/u.u1-1.5*(g-1)*v*v,  g*v;

    return m;
}

void c_Species::write_roe_jacobians(Matrix3d &left_m, Matrix3d &right_m, int interface) {

    AOS ravg = get_roe_averages(interface);

    //Matrix3d A_tilda = get_roe_matrix_abs(ravg);
    Matrix3d A_tilda2 = get_roe_matrix_abs(interface);
    Matrix3d  A_left = get_exact_Jacobian(u[interface]);
    Matrix3d A_right = get_exact_Jacobian(u[interface+1]);

    //cout<<" returnting Atilda and Atilda2: "<<endl<<A_tilda<<endl<<endl<<A_tilda2<<endl;
    left_m  = A_left  * 0.5 + A_tilda2 * 0.5;
    right_m = A_right * 0.5 - A_tilda2 * 0.5;
}

//
// Takes an AOS object, assuming the Roe-averaged states q_tilde, u_tilde, H_tilde are written into them
// Returns the absolute value of the Roe matrix for this interface.
//
Matrix3d c_Species::get_roe_matrix_abs(int j) 
{
    int jleft = j, jright = j+1;
    AOS_prim prim_l  = this->prim[jleft];
    AOS_prim prim_r  = this->prim[jright];
    
    AOS state_l      = AOS(prim_l.density, prim_l.speed*prim_l.density, prim_l.density*prim_l.internal_energy + 0.5*prim_l.density*prim_l.speed*prim_l.speed); //u[jleft];
    AOS state_r      = AOS(prim_r.density, prim_r.speed*prim_r.density, prim_r.density*prim_r.internal_energy + 0.5*prim_r.density*prim_r.speed*prim_r.speed); // = u[jright];
    
    AOS flux_l       = exact_flux(state_l);
    AOS flux_r       = exact_flux(state_r);
    
    double drho = state_r.u1 - state_l.u1;
    double dp   = prim_r.pres - prim_l.pres;
    double du   = prim_r.speed - prim_l.speed;
    
    //Roe averages
    double rho = std::sqrt(state_l.u1*state_r.u1);
    double u   = (std::sqrt(state_l.u1) * prim_l.speed + std::sqrt(state_r.u1) * prim_r.speed )/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double hl  = (prim_l.pres + state_l.u3)/state_l.u1;
    double hr  = (prim_r.pres + state_r.u3)/state_r.u1;
    double h   = (std::sqrt(state_l.u1) * hl + std::sqrt(state_r.u1) * hr)/(std::sqrt(state_l.u1) + std::sqrt(state_r.u1));
    double c   = std::sqrt((gamma_adiabat-1.) * (h-0.5*u*u));
 
    double lambdas[3] = {std::fabs(u-c), std::fabs(u), std::fabs(u+c) };

    double delta = 0.8*c;
    //HH-Entropy fix
    if(lambdas[1] < delta)
        lambdas[1] = (lambdas[1]*lambdas[1] + delta*delta)/(2.*delta); 

    double d[3] = {drho, state_r.u2 - state_l.u2, state_r.u3 - state_l.u3};
    double alphas[3];
    
    alphas[0] = (dp - rho * c * du)/(2*c*c);
    alphas[1] = drho - dp/(c*c);
    alphas[2] = (dp + rho * c * du)/(2*c*c);
    
    Vector3d ev0(1, u-c, h - u*c);
    Vector3d ev1(1, u, 0.5*u*u);
    Vector3d ev2(1, u+c, h + u*c);

    Eigen::DiagonalMatrix<double, 3> m(lambdas[0], lambdas[1], lambdas[2]);

    Matrix3d R;
    R << ev0,ev1,ev2;
    Matrix3d Rm1 = R.inverse();

    return R*(m*Rm1);
}
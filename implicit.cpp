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
        double S_l = base->surf[j-1];
        double S_r = base->surf[j];
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

        // Face velocities (from non-advanced timestep)
        double v_l = 0.5 * (prim[j-1].speed + prim[j].speed); //(u[j-1].u2 / u[j-1].u1 + u[j].u2 / u[j].u1);
        double v_r = 0.5 * (prim[j].speed + prim[j+1].speed); // ;(u[j].u2 / u[j].u1     + u[j+1].u2 / u[j+1].u1);

        //cout<<" j = "<<j<<" old rho  = "<<u[j].u1<<" / "<<" vl / vr = "<<v_l<<" / "<<v_r<<endl;
        //cout<<" j = "<<j<<" old mom  = "<<u[j].u2<<" / "<<" vl / vr = "<<v_l<<" / "<<v_r<<endl;

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
        if(0==0) {
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
        drho = -sgn * slope_rho[j+sw]  / (u[j+1].u1  - u[j].u1  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
        dmom = -sgn * slope_mom[j+sw]  / (u[j+1].u2  - u[j].u2  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); 
        dE   = -sgn * slope_E[j+sw]    / (u[j+1].u3  - u[j].u3  + 1e-50) * (x_i[ j ] - x_iVC[j+sw]); 

        //rho
        a = v_r<0? 1.0 + 0.5  * drho : -0.5 * drho;
        b = v_r<0? -0.5 * drho : 1.0 + 0.5  * drho;
        dd[idx + 0]    += theta     * lam_r.u1 * b;
        uu[idx + 0]    += theta     * lam_r.u1 * a;
        r[idx_r]       -= (1-theta) * lam_r.u1 * ( a * u[j+1].u1 + b * u[j].u1 );

        //momentum,
        //       /*
        a = v_r<0? 1.0 + 0.5  * dmom : -0.5 * dmom;
        b = v_r<0? -0.5 * dmom : 1.0 + 0.5  * dmom;
        dd[idx + 4]    += theta     * v_r * S_r * b; //mom_j+1/2
        uu[idx + 4]    += theta     * v_r * S_r * a;
        r[idx_r+1]       -= (1-theta) * v_r * S_r * ( a * u[j+1].u2 + b * u[j].u2);
        if(j>pwall && j<num_cells) {
            r[idx_r+1]       -=  0.5 * S_r * (prim_l[j+1].pres + prim_r[j].pres )   ; //momentum,
        }
                                                     
        //        */
        //Energy
        a  = v_r<0? 1.0 + 0.5  * dE : -0.5 * dE;
        b  = v_r<0? -0.5 * dE : 1.0 + 0.5  * dE;
        dd[idx + 8]    += theta     * v_r * S_r * b; //E_j+1/2
        uu[idx + 8]    += theta     * v_r * S_r * a;
        r[idx_r+ 2]    -= (1-theta) * v_r * S_r * ( a * u[j+1].u3        + b * u[j].u3);
        if(j>=ewall && j<num_cells)
            r[idx_r+ 2]    -=   0.5 * v_r * S_r * ( prim_l[j+1].pres + prim_r[j].pres ); //p div v

        //////////////////////////////////////////////////////////////////////
        // Left face
        //////////////////////////////////////////////////////////////////////
        sw  = v_l>0? -1 : 0;
        sgn = v_l>0 ? +1 : -1;
        ///////////////////////////////////

        drho = -sgn * slope_rho[j+sw]  / (u[j].u1  - u[j-1].u1  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
        dmom = -sgn * slope_mom[j+sw]  / (u[j].u2  - u[j-1].u2  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]);
        dE   = -sgn * slope_E[j+sw]    / (u[j].u3  - u[j-1].u3  + 1e-50) * (x_i[ j-1 ] - x_iVC[j+sw]);

        //rho
        a = v_l<0? 1.0 + 0.5  * drho : -0.5 * drho;
        b = v_l<0? -0.5 * drho : 1.0 + 0.5  * drho;
        dd[idx + 0]    -= theta     * lam_l.u1 * a;
        ll[idx + 0]    -= theta     * lam_l.u1 * b;
        r[idx_r]       += (1-theta) * lam_l.u1 * ( a * u[j].u1 + b * u[j-1].u1 );
        
        //momentum,
        //       /*
        a = v_l<0? 1.0 + 0.5  * dmom : -0.5 * dmom;
        b = v_l<0? -0.5 * dmom : 1.0 + 0.5  * dmom;
        dd[idx + 4]    -= theta     * v_l * S_l * a; //mom_j-1/2
        ll[idx + 4]    -= theta     * v_l * S_l * b;
        r[idx_r+1]     += (1-theta) * v_l * S_l * ( a * u[j].u2 + b * u[j-1].u2) ;
        if(j>pwall && j<=num_cells) {
                r[idx_r+1]     +=  0.5 * S_l * ( prim_l[j].pres + prim_r[j-1].pres ) ; 
        }
        // */
        //energy
        a  = v_l<0? 1.0 + 0.5  * dE : -0.5 * dE;
        b  = v_l<0? -0.5 * dE : 1.0 + 0.5  * dE;                                                
        dd[idx + 8]    -= theta     * v_l * S_l * a; //E_j-1/2
        ll[idx + 8]    -= theta     * v_l * S_l * b;
        r[idx_r+2]     += (1-theta) * v_l * S_l * ( a * u[j].u3        + b * u[j-1].u3);         //u grad E
        if(j>ewall && j<=num_cells)
            r[idx_r+2]     +=   0.5 * v_l * S_l * ( prim_l[j].pres + prim_r[j-1].pres ); //p div v

        //cout<<" j = "<<j<<" tmp_l = "<<tmp_l<<" tmp_r "<<tmp_r<<" sum "<<tmp_l+tmp_r<<endl;
        //cout<<" j = "<<j<<" ap = "<<ap<<" bp "<<bp<<" sum "<<ap+bp<<endl;

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

        //////////////////////////////////////////////////////////////////////
        // Proof of concept in big, slow matrix for only advection
        //////////////////////////////////////////////////////////////////////
        /*
        adv_mat(j,j) += V / dt;
        adv_b(j)     += V / dt * u[j].u1;
        adv_mat(j,j)   += 0.5 * lam_r.u1;
        adv_mat(j,j-1) -= 0.5 * lam_l.u1;
        adv_b(j)       -= 0.5 * ( lam_r.u1 * u[j].u1 - lam_l.u1 * u[j-1].u1 ) ;*/
    }

    //cout<<"Boundaries"<<endl;
    // More Boundaries:
    // Left boundary:
    //    Reflecting / no flux or planetary temperature
    for (int j=0; j < base->num_ghosts; j++) {
        //int idx   = j*stride  ;
        //int idx_r = j*num_vars ;
        //ll[idx] = 0 ;
        //uu[idx] = -dd[idx] ;
        //r[idx_r] = 0 ; 
    }
    
    /*
    int cnt = 0;
    cout<<" dd :"<<endl;
    for (auto i: dd){
                std::cout << i << ' ';
                cnt++;
                if(cnt%9==0) cout<<endl;
        }
    cout<<endl;

    cout<<" ll :"<<endl;
    for (auto i: ll) {
                std::cout << i << ' ';
                cnt++;
                if(cnt%9==0) cout<<endl;
    }
    cout<<endl;

    cout<<" uu :"<<endl;
    for (auto i: uu) {
        std::cout << i << ' ';
        cnt++;
        if(cnt%9==0) cout<<endl;        
    }
    cout<<endl;

   cout<<"r :"<<endl;
    for (auto i: r) {
                std::cout << i << ' ';
                cnt++;
                if(cnt%3==0) cout<<endl;
        }
    cout<<endl;
    */
    
    //cout<<"Solving"<<endl;
    //
    // Solve!
    //

    base->implicit_tridiag.factor_matrix(&ll[0], &dd[0], &uu[0]) ;
    base->implicit_tridiag.solve(&r[0], &r[0]) ; // Solve in place

    //base->LUchem_ptr[0].compute(adv_id + adv_mat.transpose()) ;
    //base->LUchem_ptr[0].compute(adv_id + adv_mat) ;
    //results.noalias() = base->LUchem_ptr[0].solve(adv_b);

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
        double vimpl_old = u[j].u2/u[j].u1;
        double vimpl_new = momtmp/rhotmp;
        double E_implied = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * prim[j].temperature;      //This is the old temperature
        double E_floor   = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * base->temperature_floor;
        double T_pred = (u[j].u3 - 0.5*momtmp*momtmp/rhotmp ) / (rhotmp * cv);
        if(base->steps == 440 || base->steps == 2700) {
            cout<<" j = "<<j<<" old/new rho  = "<<u[j].u1<<" / "<<rhotmp<<" t_damp = "<<t_damp[j]<<" corresponding f = "<<ff[j]<<endl;
            cout<<" j = "<<j<<" old/new mom  = "<<u[j].u2<<" / "<<r[j*num_vars + 1]<<" v_implied old/new = "<<vimpl_old<<"/"<<vimpl_new<<endl;
            cout<<" j = "<<j<<" old/new E    = "<<u[j].u3<<" / "<<r[j*num_vars + 2]<<endl;
            cout<<" j = "<<j<<" predicted T  = "<<T_pred<<" current T "<<prim[j].temperature<<" u1/u2/u3: "<<rhotmp<<" / "<<momtmp<<" / "<<r[j*num_vars + 2]<<" E_implied "<<E_implied<<endl;
        }

        if(base->steps>2) {
            double rhonew = std::max(r[idx_r + 0], 1e-60 );//std::max(results(j), 1e-20 ); //std::max(r[idx_r + 0], 1e-20 );//r[idx_r + 0];
            double momnew = r[idx_r + 1];
            double Enew   = r[j*num_vars + 2];
            if(Enew < 0.)
                Enew = E_floor;
            
            u[j].u1 = rhonew;//r[j*num_vars + 0];
            u[j].u2 = momnew;
            //u[j].u3 = 0.5*momnew*momnew/rhonew + rhonew * enew;//r[j*num_vars + 2] ;
            //u[j].u3 = Enew ;
            u[j].u3 = E_implied ; //This amounts to the approximation of local isothermality over the timestep. can lead to -p otherwise.
        }else {
            u[j].u3 = E_implied;
        }
    }

    // Update energies. 
    eos->compute_primitive(&(u[0]), &(prim[0]), base->num_cells+2) ;   
    eos->compute_auxillary(&(prim[0]), base->num_cells+2);

    if(base->steps > 524)
        cout<<"";
}
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
    
    Matrix_t adv_mat      = Matrix_t::Zero(num_cells+2, num_cells+2);
    Matrix_t adv_id       = Matrix_t::Identity(num_cells+2, num_cells+2);
    Vector_t adv_b        = Vector_t(num_cells+2);
    Vector_t results      = Vector_t(num_cells+2);

    //Allocate matrices
    int num_vars = 3; // + num_species
    int stride = num_vars * num_vars ;
    int size_r = (base->num_cells + 2) * num_vars ;
    int size_M = (base->num_cells + 2) * stride ;
    const int wall = 5;

    std::vector<double> 
        ll(size_M, 0.), dd(size_M, 0.), uu(size_M, 0.), r(size_r, 0.) ;
    
    ////////////////////////////////////////////////////////////////////////
    // precompute MUSCL slopes
    ////////////////////////////////////////////////////////////////////////
    std::vector<double> slope_rho(num_cells+2, 0.0);

    const std::vector<double>& 
            x_i = base->x_i, 
            x_iVC = base->x_iVC,
            dx = base->dx ;
    
    for (int j = 2; j < num_cells-1; j++) {

        double cF = (x_iVC[j+1] - x_iVC[j]) / (x_i[j] - x_iVC[j]) ;
        double cB = (x_iVC[j] - x_iVC[j-1]) / (x_iVC[j] - x_i[j-1]) ;

        double dxF = (x_iVC[j+1] - x_iVC[j]) ;
        double dxB = (x_iVC[j] - x_iVC[j-1]) ;

        slope_rho[j] = reconstruct_pointer(
                 prim[j-1].density, prim[j].density, prim[j+1].density, cF, cB, dxF, dxB) ;
    }

    if(base->steps == 49) {

        for (int j = 1; j < num_cells-1; j++) {
            cout<<" j / slope/ drho/ rho "<<j<<" "<<slope_rho[j]<<" "<<slope_rho[j]*(x_i[ j ] - x_iVC[j])<<" "<<prim[j].density<<endl;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////
    for (int j=0; j < num_cells+1; j++) {

        double V   = base->vol[j];
        double S_l = base->surf[j-1];
        double S_r = base->surf[j];
        //cout<<" writing matrix j "<<j<<endl;
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Time dependent terms:
        //rho
        adv_mat(j,j) += V / dt;
        adv_b(j)     += V / dt * u[j].u1;

        dd[idx]      += V / dt ;
        r[idx_r]     += V / dt * u[j].u1;

        //momentum
        dd[idx + 4]  += V / dt ;
        r[idx_r+ 1]  += V / dt * u[j].u2;

        //energy
        dd[idx + 8]  += V / dt ;
        r[idx_r+ 2]  += V / dt * prim[j].internal_energy;

        // Face velocities (lagged)
        double v_l = -3e7; //0.5 * (u[j-1].u2 / u[j-1].u1 + u[j].u2 / u[j].u1);
        double v_r = -3e7; //0.5 * (u[j].u2 / u[j].u1     + u[j+1].u2 / u[j+1].u1);

        if(0==1) {
            if(j<wall) {
                v_l = 0;
                v_r = 0;
            }
            if(j==wall) {
                v_l = 0.;
            }
        }
        if(j==num_cells)
            v_r = 0;
        if(j>num_cells) {
            v_l = v_r = 0;
        }
        
        //General left and right fluxes, flux = lam * density
        AOS lam_l = AOS(v_l * S_l, 0., 0.);
        AOS lam_r = AOS(v_r * S_r, 0., 0.);
        AOS lam_l_ex;
        AOS lam_r_ex;
        
        //dd[idx + 0] += v_l * S_l;

        ////////////////////////////////////////////////////////////////////////
        // Right face
        ////////////////////////////////////////////////////////////////////////
        double rho_r = 0;
        double drho  = 0, drho2  = 0;
        if(v_r > 0) {

            drho = - slope_rho[j]  / (u[j+1].u1  - u[j].u1  + 1e-50) * (x_i[ j ] - x_iVC[j]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson

            double a = 1.0 + 0.5  * drho;
            double b = -0.5 * drho;

            uu[idx + 0] += 0.5 * lam_r.u1  * b;
            dd[idx + 0] += 0.5 * lam_r.u1  * a;

            r[idx_r] -= 0.5 * lam_r.u1 * ( b * u[j+1].u1 + a * u[j].u1 );
            
            // prim_l[i].density +=  slope * (x_i[j-1] - x_iVC[j]) ; 
            // prim_r[i].density +=  slope * (x_i[ j ] - x_iVC[j]) ;

            //dd[idx + 0] += 0.5 * (lam_r.u1 );
            //rho_r = u[j].u1;

        }
        else {

            drho = +slope_rho[j+1]  / (u[j+1].u1  - u[j].u1  + 1e-50) * (x_i[ j ] - x_iVC[j+1]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson

            double a = 1.0 + 0.5  * drho;
            double b = -0.5 * drho;

            uu[idx + 0]    += 0.5 * lam_r.u1 * a;
            dd[idx + 0]    += 0.5 * lam_r.u1 * b;

            r[idx_r] -= 0.5 * lam_r.u1 * ( a * u[j+1].u1 + b * u[j].u1 );
            //rho_r = u[j+1].u1;
            
        }
        //lam_r_ex = AOS(lam_r.u1 * rho_r ,0.,0.);

        //////////////////////////////////////////////////////////////////////
        // Left face
        //////////////////////////////////////////////////////////////////////
        double rho_l = 0;
        drho  = 0;
        drho2 = 0;
        if(v_l > 0) {
            drho = - slope_rho[j-1]  / (u[j].u1  - u[j-1].u1  + 1e-50) * (x_i[ j-1 ] - x_iVC[j-1]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
            
            double a = 1.0 + 0.5  * drho;
            double b = -0.5 * drho;

            dd[idx + 0] -= 0.5 * lam_l.u1  * b;
            ll[idx + 0] -= 0.5 * lam_l.u1  * a;

            r[idx_r] += 0.5 * lam_l.u1 * ( b * u[j].u1 + a * u[j-1].u1 );

            //orig:
            //ll[idx + 0]    -= 0.5 * (lam_l.u1 );
            //rho_l = u[j-1].u1;
        }
        else{
            drho = +slope_rho[j]  / (u[j].u1  - u[j-1].u1  + 1e-50) * (x_i[ j-1 ] - x_iVC[j]); //Note: Slope 0 reduces this to the old, first order Crank-Nicolson
            
            double a = 1.0 + 0.5  * drho;
            double b = -0.5 * drho;

            dd[idx + 0]    -= 0.5 * lam_l.u1 * a;
            ll[idx + 0]    -= 0.5 * lam_l.u1 * b;

            r[idx_r] += 0.5 * lam_l.u1 * ( a * u[j].u1 + b * u[j-1].u1 );

           if(base->steps == -849) {
                cout<<"left  j, a, b, rho "<<j<<" "<<a<<" "<<b<<" "<<u[j].u1<<" drho, slope = "<<drho<<" "<<slope_rho[j]<<endl;
            }

            //rho_l = u[j].u1;
        }
        //lam_l_ex = AOS(lam_l.u1 * rho_l, 0.,0.);

        //////////////////////////////////////////////////////////////////////
        // Crank-Nicolson terms
        //////////////////////////////////////////////////////////////////////
        //r[idx_r]         -= 0.5 * ( lam_r_ex.u1 - lam_l_ex.u1 ); 

        //////////////////////////////////////////////////////////////////////
        // Proof of concept in big, slow matrix for only advection
        //////////////////////////////////////////////////////////////////////
        adv_mat(j,j)   += 0.5 * lam_r.u1;
        adv_mat(j,j-1) -= 0.5 * lam_l.u1;
        adv_b(j)       -= 0.5 * ( lam_r.u1 * u[j].u1 - lam_l.u1 * u[j-1].u1 ) ;
    }

    //cout<<"Boundaries"<<endl;
    // More Boundaries:
    // Left boundary:
    //    Reflecting / no flux or planetary temperature
    for (int j=0; j < base->num_ghosts; j++) {
        int idx   = j*stride  ;
        int idx_r = j*num_vars ;
                
        //ll[idx] = 0 ;
        //uu[idx] = -dd[idx] ;
        //r[idx_r] = 0 ; 
    }
    
    //cout<<" dd and r "<<endl;
    
    /*
    for (auto i: dd)
        std::cout << i << ' ';
    cout<<endl;

    for (auto i: r)
        std::cout << i << ' ';
    cout<<endl;

    */

    //cout<<"Solving"<<endl;
    //
    // Solve!
    //
    base->implicit_tridiag.factor_matrix(&ll[0], &dd[0], &uu[0]) ;
    base->implicit_tridiag.solve(&r[0], &r[0]) ; // Solve in place

    //base->LUchem_ptr[0].compute(adv_id + adv_mat.transpose()) ;
    base->LUchem_ptr[0].compute(adv_id + adv_mat) ;
    results.noalias() = base->LUchem_ptr[0].solve(adv_b);

    if(base->steps==100) {
        cout<<" mat results: "<<endl;
        for (auto i: results)
            std::cout << i << ' ';
        cout<<endl;
    }

    //
    // End Solve
    //

    //
    // Write solution back into variables
    //
    for (int j=0; j <= base->num_cells; j++) {
        int idx   = j*stride  ;
        int idx_r = j*num_vars;


        if(base->steps>2) {
            double rhonew = std::max(r[idx_r + 0], 1e-20 );//std::max(results(j), 1e-20 ); //std::max(r[idx_r + 0], 1e-20 );//r[idx_r + 0];
            double momnew = r[idx_r + 1];
            double enew   = cv * 500;//r[j*num_vars + 2]
            
            
            //cout<<" j = "<<j<<" old/new rho  = "<<u[j].u1<<" / "<<r[idx_r + 0]<<endl;
            //cout<<" j = "<<j<<" old/new mom  = "<<u[j].u2<<" / "<<r[j*num_vars + 1]<<endl;
            //cout<<" j = "<<j<<" old/new E    = "<<u[j].u3<<" / "<<r[j*num_vars + 2]<<endl;

            u[j].u1 = rhonew;//r[j*num_vars + 0];
            u[j].u2 = 0.; //momnew;
            u[j].u3 = 0.5*momnew*momnew/rhonew + rhonew * enew;//r[j*num_vars + 2] ;
        }
    }

    // Update energies. 
    eos->compute_primitive(&(u[0]), &(prim[0]), base->num_cells+2) ;   
    eos->compute_auxillary(&(prim[0]), base->num_cells+2);
    //    species[si].eos->compute_conserved(&(species[si].prim[0]), &(species[si].u[0]), num_cells+2);        

}
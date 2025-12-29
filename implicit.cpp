/**
 * implicit.cpp
 * 
 * This file contains routines solving the hydrodynamic euler system with source terms in an implicit fashion
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"

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
    
    //Write into matrices
    for (int j=0; j < num_cells; j++) {

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

        adv_mat(j,j) += V / dt;
        adv_b(j)     += V / dt * u[j].u1;

        //momentum
        dd[idx + 4]  += V / dt ;
        r[idx_r+ 1]  += V / dt * u[j].u2;

        //energy
        dd[idx + 8]  += V / dt ;
        r[idx_r+ 2]  += V / dt * prim[j].internal_energy;

        // Face velocities (lagged)
        double v_l = +1e7;//0.5 * (u[j-1].u2 / u[j-1].u1 + u[j].u2 / u[j].u1);
        double v_r = +1e7;//0.5 * (u[j].u2 / u[j].u1     + u[j+1].u2 / u[j+1].u1);

        if(0==1) {
            if(j<wall) {
                v_l = 0;
                v_r = 0;
            }
            if(j==wall) {
                v_l = 0.;
            }

        }
        
        int offset = +9;

        //General left and right boundary fluxes
        AOS lam_l = AOS(v_l * S_l, 0., 0.);
        AOS lam_r = AOS(v_r * S_r, 0., 0.);
        AOS lam_l_ex;
        AOS lam_r_ex;
        
        //dd[idx + 0] += v_l * S_l;
        if(v_r > 0) {
            dd[idx + 0]    += 0.5 * lam_r.u1;
            
            lam_r_ex = AOS(lam_r.u1 * u[j].u1 ,0.,0.);
        }
        else {
            uu[idx + 0]    += 0.5 * lam_r.u1;

            lam_r_ex = AOS(lam_r.u1 * u[j+1].u1 ,0.,0.);
        }

        if(v_l > 0) {
            ll[idx + 0]    -= 0.5 * lam_l.u1;

            lam_l_ex = AOS(lam_l.u1 * u[j-1].u1, 0.,0.);
        }
        else{
            dd[idx + 0]    -= 0.5 * lam_l.u1;

            lam_l_ex = AOS(lam_l.u1 * u[j].u1, 0.,0.);
        }
        
        r[idx_r]         -= 0.5 * ( lam_r_ex.u1 - lam_l_ex.u1 ); //Crank-Nicholson

        adv_mat(j,j)   += 0.5 * lam_r.u1;
        adv_mat(j,j-1) -= 0.5 * lam_l.u1;
        adv_b(j)       -= 0.5 * ( lam_r.u1 * u[j].u1 - lam_l.u1 * u[j-1].u1 ) ;
        
        //adv_b(j)       += v_l * S_l * 0.5*(u[j].u1+u[j-1].u1);
        //adv_b(j+1)     -= v_l * S_l * 0.5*(u[j].u1+u[j-1].u1);

        // RIGHT FACE (outgoing from j)
        //dd[idx + 0] += v_r * S_r;
        //uu[idx + 0] -= v_r * S_r;
        /*
        if (v_l > 0.0) {
                // Flux uses U_j
                dd[idx + 0] +=  v_l * S_l;
                dd[idx + 4] +=  v_l * S_l;
                dd[idx + 8] +=  v_l * S_l;
        } else {
                // Flux uses U_{j-1}
                ll[idx + 0] +=  v_l * S_l;
                ll[idx + 4] +=  v_l * S_l;
                ll[idx + 8] +=  v_l * S_l;
        }

            // -----------------------
            // RIGHT FACE (j+1/2)
            // -----------------------
        if (v_r > 0.0) {
                // Flux uses U_j
                dd[idx + 0] -=  v_r * S_r;
                dd[idx + 4] -=  v_r * S_r;
                dd[idx + 8] -=  v_r * S_r;
        } else {
                // Flux uses U_{j+1}
                uu[idx + 0] +=  v_r * S_r;
                uu[idx + 4] -=  v_r * S_r;
                uu[idx + 8] -=  v_r * S_r;
        }*/

       

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


    //Eigen::PartialPivLU<Matrix_t> LUadv;
    //
    //LUadv          = Eigen::PartialPivLU<Matrix_t>;
    //LUchem_ptr[i]          = Eigen::PartialPivLU<Matrix_t>;
    //LUadv.compute(adv_id + adv_mat.transpose()) ;
    //results.noalias() = LUadv.solve(adv_b);

    //base->LUchem_ptr[0].compute(adv_id + adv_mat.transpose()) ;
    base->LUchem_ptr[0].compute(adv_id + adv_mat) ;
    results.noalias() = base->LUchem_ptr[0].solve(adv_b);

    if(base->steps==100) {
        cout<<" mat results: "<<endl;
        for (auto i: results)
            std::cout << i << ' ';
        cout<<endl;
    }

    //    LUchem_ptr[loc_thr].compute(identity_matrix + reaction_matrix_ptr[loc_thr].transpose()) ;
    //n_news.noalias() = LUchem_ptr[loc_thr].solve(reaction_b_ptr[loc_thr]);
    
    //
    // End Solve
    //

    //cout<<"Writing"<<endl;
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
            
            //cout<<" j = "<<j<<" old/new rho  = "<<u[j].u1<<" / "<<results(j)<<endl;
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
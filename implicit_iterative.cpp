/**
 * implicit_iterative.cpp
 * 
 * This file contains the iterative hydrodynamic solver
 */

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "aiolos.h"


/**
 * Computes the implicit hydro solution for incompressible fluids - applicable to gases e.g. at low mach numbers for electrons
 * This "version 2" acts on the internal energy, whereas the previous, above solver solves the total energy
 * 
 * @param sp: the species of which to compute the timestep advance
 */
void c_Species::implicit_iterative(double dt) {
/* 
    cout<<" Hi i am species "<<speciesname<<" and i am attempting to be solved implicitly."<<endl;

    char a;
    cin>>a; */
    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////

    int num_iter_max = 10;
    double convergence_measure = 0;
    double epsilon_desired = 0;

 //cout<<" Hi i am species "<<speciesname<<" and i am attempting to be solved implicitly."<<endl;

    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////

    //Boundaries
    this->apply_boundary_left(this->u) ;
    this->apply_boundary_right(this->u) ;
    compute_pressure(this->u);
    std::vector<double> u_tmp = np_zeros(num_cells+2); 
    std::vector<AOS>    u_n   = std::vector<AOS>(num_cells+2); 
    std::vector<AOS>    u_k   = std::vector<AOS>(num_cells+2);
    reconstruct_edge_states(u_mask, 2) ;

    std::vector<double> t_damp = np_zeros(num_cells+2); 
    std::vector<double> ff     = np_zeros(num_cells+2); 

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
    // Write initial data
    ////////////////////////////////////////////////////////////////////////

    for(int i=0; i<num_cells; i++) {
        u_k[i] = u[i];
        u_n[i] = u[i];
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Construct implicit Matrix
    ////////////////////////////////////////////////////////////////////////
    for (int j=0; j <= num_cells+1; j++) {

        double V   = base->vol[j];
        int idx   = j*stride  ;
        int idx_r = j*num_vars;

        // Time dependent terms:
        //rho
        dd[idx]      += V / dt ;
        //momentum
        dd[idx + 4]  += V / dt ;
        //energy
        dd[idx + 8]  += V / dt ;

    }
    /////////////////////////////////////////////////////////////////////////
    // Precompute Jabobians interface by interface, define flux function pointer
    //////////////////////////////////////////////////////////////////////////

    for (int j=1; j <= num_cells; j++) { //Going interface by interface - interface j is between cell j and j+1. Ignore interface 0, that's just excess memory to get the index shifted by 1
        (this->*write_jacobians)(base->impl_jacobian_left[j], base->impl_jacobian_right[j], prim_to_u(prim_r[j]),prim_to_u(prim_l[j+1])); //u[j], u[j+1]
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

        int imin = base->ignore_electron_cfl_cell; // default is 1, default for ignore_electron_cfl_cell is also 1
        //determine imin based on minimum mixing ratio, to avoid instabilities of electrons at low densities
        double f =0;
        for (int j=1; j <= num_cells; j++) {  
            f = prim[j].number_density/base->total_numdens[j];
            imin = j;
            if(f > base->mix_p3)
                break;
        }
        if(base->steps%1000==0)
            cout<<" limiting cell in implicit solver found as "<<imin<<endl;

        for (int j=imin; j <= num_cells; j++) {  //Going cell by cell. Reminder: cell 1 is boundaed by interface 0 (left) and 1 (right), cell 2 by 1 and 2
            double V   = base->vol[j];
            double S_l = base->surf[j-1];
            double S_r = base->surf[j];
            //cout<<" writing matrix j "<<j<<endl;
            int idx   = j*stride  ;
            int idx_r = j*num_vars;

            double gsrc    = source_grav(u[j], j).u2;
            double gsrc3   = source_grav(u[j], j).u3;
            double gsrc_no = source_grav_noconserved(j).u2;
            //double psrc    =  -(base->source_pressure_prefactor_left[j] * prim_l[j].pres - 
            //                                  base->source_pressure_prefactor_right[j] * prim_r[j].pres);
            
            if(impl_debug)                                
                cout<<"starting cell j = "<<j<<" steps "<<base->steps<<endl;
            ////////////////////////////////////////////////////////////////////////
            // Right face
            ////////////////////////////////////////////////////////////////////////
            
            Vector3d ul             = u_to_vec(prim_to_u(prim_r[j]));
            Vector3d ur             = u_to_vec(prim_to_u(prim_l[j+1]));
            Matrix3d Jl             = base->impl_jacobian_left[j];
            Matrix3d Jr             = base->impl_jacobian_right[j];
            Vector3d dul            = Jl * ul;
            Vector3d dur            = Jr * ur;
            Vector3d f_last         = (this->*flux_pointer)(prim_to_u(prim_r[j]),prim_to_u(prim_l[j+1]));

            if(impl_debug) {
                cout<<" in implicit roe solver Jl and Jr ="<<endl<<Jl<<endl<<endl<<Jr<<endl;
                cout<<" left and right states: "<<endl<<ul<<endl<<ur<<endl;
                cout<<" f_last "<<f_last<<endl;
            }
            //Write R.H.S components
            for(int ii=0; ii<3; ii++) {
                //r[idx_r + ii] +=  S_r * (- f_last(ii) + dul(ii) + dur(ii));

                //if(impl_debug)                                
                //    cout<<" was writing into ii = "<<ii<<" r elms, f_roe = "<<f_last(ii) <<" dul+dur = "<<dul(ii) + dur(ii)<<endl;

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
            
            ul             = u_to_vec(prim_to_u(prim_r[j-1]));
            ur             = u_to_vec(prim_to_u(prim_l[j]));
            Jl             = base->impl_jacobian_left[j-1];
            Jr             = base->impl_jacobian_right[j-1];
            dul            = Jl * ul;
            dur            = Jr * ur;
            f_last         = (this->*flux_pointer)(prim_to_u(prim_r[j-1]),prim_to_u(prim_l[j]));

            //Write L.H.S components
            for(int ii=0; ii<3; ii++) {
                //r[idx_r + ii] -=  S_l * (- f_last(ii) + dul(ii) + dur(ii));

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
                std::vector<double> diff = get_hydro_jacobian_P(u[j]);
                double prsc_total = (diff[0]*u[j].u1 + diff[1]*u[j].u2 + diff[2]*u[j].u3);
                    
                double v2r = + V * 2./base->x_iVC[j] ;
                dd[idx + 3]    -= v2r * diff[0]; 
                dd[idx + 4]    -= v2r * diff[1];  
                dd[idx + 5]    -= v2r * diff[2]; 

                //Gravity
                dd[idx + 3]    -= V * gsrc_no;             //  dS_g/drho * rho in momentum equation
                dd[idx + 7]    -= V * gsrc_no;             // 
            }
        //////////////////////////////////////////////////////////////////////
        // End Sources
        //////////////////////////////////////////////////////////////////////
    }



    //
    // Jacobian finished, factor matrix for use and reuse
    //
    base->implicit_tridiag.factor_matrix(&ll[0], &dd[0], &uu[0]) ;

    //
    // Compute Residual R(U^k, U^n)
    //
    Vector3d ul  ;
    Vector3d ur       ;
    Matrix3d Jl      ;
    Matrix3d Jr       ;
    Vector3d dul      ;
    Vector3d dur     ;
    Vector3d f_last   ;
    for(int k=0; k<num_iter_max; k++) {
        //
        // Compute Residual R(U^k, U^n)
        //
        for (int j=0; j <= num_cells+1; j++) { 
            double V   = base->vol[j];
            double S_l = base->surf[j-1];
            double S_r = base->surf[j];
            int idx   = j*stride  ;
            int idx_r = j*num_vars;
            //Evaluate residual at current U^k and write it into r
            
            r[idx_r]     = V / dt * (u_n[j].u1 );
            r[idx_r+ 1]  = V / dt * (u_n[j].u2 );
            r[idx_r+ 2]  = V / dt * (u_n[j].u3 );
        }

        for (int j=2; j <= num_cells; j++) {
            double V   = base->vol[j];
            double S_l = base->surf[j-1];
            double S_r = base->surf[j];
            int idx   = j*stride  ;
            int idx_r = j*num_vars;
            //Evaluate residual at current U^k and write it into r
            r[idx_r]     -= V / dt * (u_k[j].u1);
            r[idx_r+ 1]  -= V / dt * (u_k[j].u2);
            r[idx_r+ 2]  -= V / dt * (u_k[j].u3);
            ////////////////////////////////////Fluxes at left and right faces//////////////////////////////////
            //Vector3d f_last_r         = (this->*flux_pointer)(prim_to_u(prim_r[j]),prim_to_u(prim_l[j+1]));
            //Vector3d f_last_l         = (this->*flux_pointer)(prim_to_u(prim_r[j-1]),prim_to_u(prim_l[j]));
            //for(int ii=0; ii<3; ii++) { //Regular implicit solver has the order Left - Right
            //    r[idx_r + ii] += S_l * f_last_l(ii) - S_r * f_last_r(ii);
            //}


            //@@@@@@@@@@@@@@@@@@@@
            ul             = u_to_vec(prim_to_u(prim_r[j]));
            ur             = u_to_vec(prim_to_u(prim_l[j+1]));
            Jl             = base->impl_jacobian_left[j];
            Jr             = base->impl_jacobian_right[j];
            dul            = Jl * ul;
            dur            = Jr * ur;
            f_last         = (this->*flux_pointer)(prim_to_u(prim_r[j]),prim_to_u(prim_l[j+1]));
            //Write R.H.S components
            for(int ii=0; ii<3; ii++) {
                r[idx_r + ii] +=  S_r * (- f_last(ii));
            }

            //@@@@@@@@@@@@@@@@@@@@
            ul             = u_to_vec(prim_to_u(prim_r[j-1]));
            ur             = u_to_vec(prim_to_u(prim_l[j]));
            Jl             = base->impl_jacobian_left[j-1];
            Jr             = base->impl_jacobian_right[j-1];
            dul            = Jl * ul;
            dur            = Jr * ur;
            f_last         = (this->*flux_pointer)(prim_to_u(prim_r[j-1]),prim_to_u(prim_l[j]));
            //Write L.H.S components
            for(int ii=0; ii<3; ii++) {
                r[idx_r + ii] -=  S_l * (- f_last(ii));
            }
            //@@@@@@@@@@@@@@@@@@@@

            //////////////////////////////////////Source terms/////////////////////////////
            if(base->geometry!=Geometry::cartesian) {
                std::vector<double> diff = get_hydro_jacobian_P(u[j]);
                double prsc_total = (diff[0]*u[j].u1 + diff[1]*u[j].u2 + diff[2]*u[j].u3);

                double p_tmp    = (u[j].u3-0.5*u[j].u2*u[j].u2/u[j].u1)/(gamma_adiabat-1);
                //r[idx_r+ 1]    += V * 2./base->x_iVC[j] * (p_tmp - prsc_total ); //2P/r
            }
            if(j>pwall && j<num_cells) {
                 AOS gsrc    = source_grav(u[j], j);
                //r[idx_r+1]       += 1.0 * V * (gsrc.u2 );   //Grav and geometric source
                //r[idx_r+2]       += 1.0 * V * (gsrc.u3);              //Grav and geometric source
            }
        }
        //End Residual

        //Inverse sign
        for (int j=0; j <= num_cells+1; j++) { 
            int idx_r = j*num_vars;
            /* r[idx_r]     *= -1;  //No need because we already have accounted for the -1 by keeping our signs from the original implicit code
            r[idx_r+ 1]  *= -1;
            r[idx_r+ 2]  *= -1;   */
        }

        //Calculate residual metric before solving
        //
        double abs_res = 0;
        double dens_res = 0;
        double mom_res  = 0;
        double energy_res = 0;
        double rmax=0, mmax=0, emax=0;
        for (int j=0; j <= num_cells+1; j++) {
            double V   = base->vol[j];
            double S_l = base->surf[j-1];
            double S_r = base->surf[j];
            int idx   = j*stride  ;
            int idx_r = j*num_vars;

            if((j>2) &&  (j<num_cells) ) {

                dens_res    += std::fabs(r[idx_r + 0]) / u0[j].u1         ;  
                if(std::fabs(u0[j].u2)>0)
                    mom_res += std::fabs(r[idx_r + 1]) / std::fabs(u0[j].u2)     ;  
                energy_res  += std::fabs(r[idx_r + 2]) /  u0[j].u3        ;
                abs_res      = dens_res + energy_res;

                //AOS ttmp(r[idx_r + 0], r[idx_r + 1], r[idx_r + 2]);
                if(eint(u_k[j])<0)
                    energy_res += 1e10;

            }
            
        }
        
        if(impl_debug) {
            cout<<"steps "<<base->steps<<" iter k == "<<k<<" residual metric = "<<abs_res<<" "<<endl;
            
            cout<<"r ="<<endl; 
            int cnt=0;
            for(auto rr: r){
                cout<<rr<<" ";
                cnt++;
                if(cnt%3==0)
                    cout<<endl;
            }

            cout<<"dd ="<<endl; 
            cnt=0;
            for(auto d: dd){
                cout<<d<<" ";
                cnt++;
                if(cnt%9==0)
                    cout<<endl;
            }
            
            cout<<" just before solve "<<endl;
        }
        //char a;
        //cin>>a; 
        //cout<<"Solving"<<endl;
        //
        // Solve!
        //
        
        base->implicit_tridiag.solve(&r[0], &r[0]) ; //r and dx 
        
        //Determine some residual metric to gauge convergence
        //Break condition if convergence metric is good

        //
        // Write new solution into Uk and recompute primitive slopes
        //

        //cout<< "After solve variables "<<endl;
        for (int j=0; j <= base->num_cells; j++) {
            int idx   = j*stride  ;
            int idx_r = j*num_vars;

            double drhotmp = r[idx_r + 0];
            double dmomtmp = r[idx_r + 1];
            double dEtmp   = r[idx_r + 2];

            //Do checks on solution tmps
            
            u_k[j] += AOS(drhotmp, dmomtmp, dEtmp);

            if(j>num_cells-2) {
                u[k] = u[k-1];
            }

            if(impl_debug) 
                cout<<u_k[j].u1<<" "<<u_k[j].u2<<" "<<u_k[j].u3<<" cell "<<j<<endl;

            if(j<base->num_cells-1)
                if(eint(u_k[j])<0)
                        energy_res += 1e10;
        }

        //cout<<"steps "<<base->steps<<" old convergence metric was "<<abs_res<<" dens "<<dens_res<<" en "<<energy_res<<" it took num_iters = "<<k<<" energies[30-33] = "<<eint(u_k[30])<<" "<<eint(u_k[31])<<" "<<eint(u_k[32])<<" "<<eint(u_k[62])<<endl;
        cout<<"steps "<<base->steps<<" convergence: r "<<dens_res<<" m "<<mom_res<<" e "<<energy_res<<" it took num_iters = "<<k<<" applied dt = "<<dt<<" ens:"<<(u_k[300].u3)<<" "<<(u_k[301].u3)<<" "<<(u_k[302]).u3<<" "<<(u_k[303]).u3 <<endl;
        if( (std::fabs(abs_res) < 1e-13 ) || (k==num_iter_max) ) {
            /* cout<<" inspect solution "<<endl;
            for (int j=0; j <= base->num_cells; j++) {
                cout<<u_k[j].u1<<" "<<u_k[j].u2<<" "<<u_k[j].u3<<" cell "<<j<<" eint "<<eint(u_k[j])<<endl;
            } */

            break;
        }
        

            

        eos->compute_primitive(&(u_k[0]), &(prim[0]), base->num_cells) ;   
        eos->compute_auxillary(&(prim[0]), base->num_cells); 
        reconstruct_edge_states(u_mask, 2) ;
        //reconstruct_edge_states(u_mask, 2) ;

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

            double rhotmp = std::max(u_k[j].u1, 1e-60 ); 
            double momtmp = u_k[j].u2;
            double Etmp   = u_k[j].u3;
            double Ekin = 0.5*momtmp*momtmp/rhotmp;

            double T_tmp = prim[j].temperature = std::min(std::max(prim[j].temperature, base->temperature_floor), base->max_temperature);
            
            double rhoe_implied = Etmp-Ekin;
            double E_implied = Ekin + rhotmp * cv * T_tmp;      //This is the old temperature
            double E_floor   = 0.5*momtmp*momtmp/rhotmp + rhotmp * cv * base->temperature_floor;
            
            if(base->steps>=0) {
                double rhonew = std::max(u_k[j-1].u1, 1e-60 );//std::max(results(j), 1e-20 ); //std::max(r[idx_r + 0], 1e-20 );//r[idx_r + 0];
                double momnew = u_k[j-1].u2;
                double Enew   = Etmp; //std::min(r[j*num_vars + 2], E_floor ); // we solved for rho e, not e
                
                
                if( rhoe_implied < 0.) {
                    cout<<" negative enew in cell "<<j<<" "<<endl; 
                    if(j>num_cells-1) {
                        rhonew = u_k[j-1].u1;   
                        momnew = u_k[j-1].u2;
                        Enew = u_k[j-1].u3;
                        cout<<"fixed cell "<<j<<endl;
                    }
                        
                }
                if( Enew < 0.) {
                    cout<<" negative Enew in cell "<<j<<" "<<endl; 
                }

                u[j].u1 = rhonew;
                u[j].u2 = momnew;
                u[j].u3 = Enew;
            }

            

    }
    //cout<<" Hi im the famous boundary cell density "<<u_n[num_cells-1].u1<<endl;
    //Fix boundaries
    u[0] = u_n[0]; 
    u[1] = u_n[1]; 
    u[num_cells] = u_n[num_cells-1]; 
    u[num_cells+1] = u_n[num_cells-1]; 
    
        
    // Update primitives. 
    eos->compute_primitive(&(u[0]), &(prim[0]), base->num_cells+1) ;   
    eos->compute_auxillary(&(prim[0]), base->num_cells+1); 
/* 
    char b;
    cin>>b; */
}

double c_Species::compute_hydro_residual(std::vector<AOS>& u_tmp, std::vector<double>& res ) {

    //Allocate matrices
    int num_vars = 3; // + num_species
    int stride = num_vars * num_vars ;
    int size_r = (base->num_cells + 2) * num_vars ;
    int size_M = (base->num_cells + 2) * stride ;

    for (int j=1; j <= num_cells; j++) {

        double V   = base->vol[j];
        double S_l = base->surf[j-1];
        double S_r = base->surf[j];

        int idx   = j*stride  ;
        int idx_r = j*num_vars;


    }
}
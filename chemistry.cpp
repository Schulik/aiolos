///////////////////////////////////////////////////////////
//
//
//  chemistry.cpp
//
//  A module to solve for time-dependent chemistry and photochemistry
//
//
//
//
///////////////////////////////////////////////////////////

#define EIGEN_RUNTIME_NO_MALLOC

#include <array>
#include <cassert>

#include "aiolos.h"

double thermo_get_r(double T, double r1, double delta_stoch);

//
// Wrapper function to make it possible to call different chemistry solvers (such as specialized one-equation solvers, that might be much faster than the general solver)
//
//
//
//
//
void c_Sim::do_chemistry() {
    
    for (int j = num_cells + 1; j > 0; j--) {
        solver_cchem_implicit_general(dt, j, 0);
        
        //update_tau_in_chem();
        update_dS_jb(int j, int b);
    }
    
    //if(chemistry == 2)
    //     solver_cchem_implicit_something_specialized();
}


c_reaction::c_reaction(int reacnum, int num_species, std::vector<int> e_indices,std::vector<int> p_indices,std::vector<double> e_stoch, std::vector<double> p_stoch, double reaction_rate) 
{
        for(int e=0; e<e_indices.size(); e++) {
            if(e_indices[e] > num_species)
                throw 0;
        }
        for(int p=0; p<p_indices.size(); p++) {
            if(p_indices[p] > num_species)
                throw 1;
        }
        
        this->reaction_number    = reacnum;
        this->educts   = e_indices;
        this->products = p_indices;
        
        this->e_num   =     e_indices.size();
        this->p_num   = p_indices.size();
        
        this->e_stoch = std::vector<double>(num_species); //e_stoch;
        this->p_stoch = std::vector<double>(num_species); //p_stoch;
        this->r       = reaction_rate   ;
        
        for(int s=0; s<num_species; s++) {
            this->e_stoch[s] = 0;
            this->p_stoch[s] = 0;
        }
        
        delta_stoch = 0.;
        for(int e=0; e<e_num; e++) {
            this->e_stoch[e_indices[e]] = e_stoch[e];
            delta_stoch += e_stoch[e];
        }
            
        for(int p=0; p<p_num; p++) {
            this->p_stoch[p_indices[p]] = p_stoch[p];
            delta_stoch -= p_stoch[p];
        }
            
    }

c_photochem_reaction::c_photochem_reaction(int reacnum, int num_species, int num_bands, std::vector<int> e_indices, std::vector<int> p_indices, std::vector<double> e_stoch, std::vector<double> p_stoch, int band)
{
        //cout<<"in init photochem reaction "<<endl;
        for(int e=0; e<e_indices.size(); e++) {
            //cout<<" e = "<<e_indices[e]<<endl;
            if(e_indices[e] > num_species)
                throw 2;
        }
        for(int p=0; p<p_indices.size(); p++) {
            //cout<<" p ="<<p_indices[p]<<endl;
            if(p_indices[p] > num_species)
                throw 3;
        }
        
        //cout<<"e_indices = "<<e_indices<<" p_indices = "<<p_indices<<endl;
        //cout<<"pos1"<<endl;
        this->reaction_number = reacnum;
        this->educts   = e_indices;
        this->products = p_indices;
        //cout<<"pos2"<<endl;
        this->e_num   = e_indices.size();
        this->p_num   = p_indices.size();
        //cout<<"pos3, enum/pnum = "<<this->e_num<<" "<<this->p_num<<endl;
        this->e_stoch = std::vector<double>(num_species);//  e_stoch;
        this->p_stoch = std::vector<double>(num_species); // p_stoch;
        //cout<<"pos4"<<endl;
        
        for(int s=0; s<num_species; s++) {
            this->e_stoch[s] = 0.;
            this->p_stoch[s] = 0.;
        }
        for(int s=0; s<num_species; s++) {
            //cout<<" s = "<<s<<" "<<e_stoch[s]<<"/"<<p_stoch[s]<<endl;
        }
        //cout<<"pos5"<<endl;
        for(int e=0; e<e_num; e++) {
            this->e_stoch[e_indices[e]] = e_stoch[e];
            //cout<<" set e_stoch = "<<e_stoch[e]<<endl;
            //cout<<" set e_indices = "<<e_indices[e]<<endl;
        }
            
        for(int p=0; p<p_num; p++) {
            this->p_stoch[p_indices[p]] = p_stoch[p];
            //cout<<" set p_stoch = "<<p_stoch[p]<<endl;
            //cout<<" set p_indices = "<<p_indices[p]<<endl;
        }
            
        
        //cout<<"pos6"<<endl;
        if(band+1 > num_bands)
            throw 4;
        
        this->band = band; //TODO: Check that band <= num_bands
        
        //cout<<"pos7"<<endl;
    }


void c_Sim::init_reactions(int cdebug) {
    
    //std::vector<c_reaction> reactions;
    //std::vector<c_photochem_reaction> photoreactions;
    cout<<"starting photoreactions"<<endl;
    
    double kf = 1e0;
    double ns = num_species;
    
    //photoreactions.push_back(c_photochem_reaction( 4, ns, num_bands_in, {1}, {2,3}, {1.}, {1.,1.}, 0 ));
    
    cout<<"starting reactions"<<endl;
    
    //reactions.push_back(c_reaction(1, ns, {2,3}, {1}, {1.,1.}, {1.}, 2.7e-13 )); //Electron-proton recombination
    //reactions.push_back(c_reaction(2, ns, {0}, {1}, {1.}, {2.}, kf));            //H2 Chemistry there
    //reactions.push_back(c_reaction(3, ns, {1}, {0}, {2.}, {1.}, kf*1e-2));             //H2 Chemistry back
    reactions.push_back(c_reaction(0, ns, {0,1}, {2,3}, {3.,1.}, {1.,1.}, 1e-40 )); //CO + 3 H2 -> CH4 + H2O
    reactions.push_back(c_reaction(0, ns, {2,3}, {0,1}, {1.,1.}, {3.,1.}, 1e-30 )); //CO + 3 H2 <- CH4 + H2O
    
    //if(cdebug > 0) {
        
        cout<<" Reporting reactions in init chemistry... "<<endl;
        
        for(c_photochem_reaction& reaction : photoreactions) {
            cout<<" Photoreaction "<<reaction.reaction_number<<": ";
            
            for(int& ei : reaction.educts) {
                cout<<int(reaction.e_stoch[ei]*1.01)<<" "<<species[ei].speciesname<<" + ";
            }
            cout<<" hv -> ";
            for(int& pi : reaction.products) {
                cout<<int(reaction.p_stoch[pi]*1.01)<<" "<<species[pi].speciesname<<" + ";
            }
        }
        cout<<endl;
        
        for(c_reaction& reaction : reactions) {
            cout<<" Reaction "<<reaction.reaction_number<<": ";
            
            for(int& ei : reaction.educts) {
                cout<<int(reaction.e_stoch[ei]*1.01)<<" "<<species[ei].speciesname<<" + ";
            }
            cout<<" -> ";
            for(int& pi : reaction.products) {
                cout<<int(reaction.p_stoch[pi]*1.01)<<" "<<species[pi].speciesname<<" + ";
            }
            cout<<" #### dStoch = "<<reaction.delta_stoch<<endl;
        }
        
    //char a;
    //cin>>a;
        
    //}
    
    
}


//
//
// Here we split all reactions into forward, backward and individual summation steps.
// Each term in a dn/dt equation is approximate via the semi-implicit Euler approximation. 
// k denotes the timestep, k+1 being the advanced one, A and B are number densities for species A and B
//
// 1/s * dn/dt = f(x^k+1) = k_r * (  (A^k+1)^a * (B^k+1)^b \approx f(x^k) + sum_educts dfdn_e * (e^k+1 - e^k )  )
// where we define labels 
//                                                      t1     + t2_e1 +  t3_e1                         +t2_e2 + t3_e2 + ....
// 
// It is t1 = f(x^k), and f1 and all f3 terms need to appear in the b-vector element i, as b_i = stoch[i] (t1 + sum_e t3_e)
//                      whereas all t2-terms go on the left hand side of the equation, as they are multiplied by the advanced time educt e^k+1
//                      they will find their ways into the matrix.
//
//
int c_Sim::solver_cchem_implicit_general(double dt, int cell, int cdebug) {
    
    //cout<<"Starting solver_chem_implicit_general."<<endl;
    
    reaction_matrix = Matrix_t::Zero(num_species, num_species);
    reaction_b      = Vector_t(num_species);
    n_olds          = Vector_t(num_species);
    double Tcell           = species[0].prim[cell].temperature;
    double n_tot = 0.;
    
    for(int s=0;s<num_species; s++) {
        n_tot += species[s].prim[cell].number_density;
    }
    
    for(int s=0;s<num_species; s++) {
        n_olds[s]     = species[s].prim[cell].number_density / n_tot;
        reaction_b[s] = n_olds[s];
        //cout<<n_olds[s]<<" ";
    }
    //cout<<endl;
    //cout<<"pos2."<<endl;    
    //
    // Photochemical reactions
    //
    //for ri,reac in enumerate(photochem_reactions) {
    
    for(c_photochem_reaction& reaction : photoreactions) {
        //cout<<"pc pos3. reaction = "<<reaction.reaction_number<<" b = "<<reaction.band<<" cell = "<<cell<<endl;
    
        int b = reaction.band;
        double temptau = 0.;
        
        double dlognu = 1.;
        double ds = dx[cell];
        //cout<<" dx "<<ds<<endl;
        //cout<<" S "<<solar_heating(b)<<endl;
        //cout<<" E "<<photon_energies[b]<<endl;
        
        if(cell < num_cells)
            temptau = radial_optical_depth_twotemp(cell+1,b);
        
        double F = 0.25 * solar_heating(b) / photon_energies[b] * std::exp(-temptau) / ds * dlognu;
        
        //cout<<"tau "<<temptau<<endl;
        //cout<<" F "<<F<<endl;
        //cout<<" kappa "<<kappa<<endl;
        //cout<<"pc pos4 "<<endl;
        //for i,educt in enumerate(reac.educt) {
        for(int& ei : reaction.educts) { //We assue that photochem reactions look likes reaction.educts[ei] + hv -> products, but leave the more general case open here
            
            double dtau = 0.;
            double dfdn = 0.;
            double dfdn_po = 0.;
            //cout<<"pc pos5, ei = "<<ei<<endl;
            double kappa = species[reaction.educts[ei]].opacity_twotemp(cell, b) * species[reaction.educts[ei]].mass_amu*amu; 
            
            dtau    = ds*kappa*n_olds[ei];
            dfdn    = dt*F*kappa*std::exp(-dtau) / n_tot;  //TODO: Replaces F, kappa, dx with the appropriate values for species and flux in cell 
            
            dfdn_po = -dfdn*n_olds[ei];  //df/dn|k * n^k
            
            //print("i/reac_educt[i]/dfdn/dfdn_po = " + repr([i,reac.educt[i], dfdn, dfdn_po]))
            //
            // This distributes t2 = dfdn*stoch throughout the matrix
            // and adds dfdn*n to b for this educt
            //
            //for j,edu in enumerate(reac.educt) {
            for(int& ej : reaction.educts ) {
                //cout<<"pc pos5, ej = "<<ej<<endl;
                //print(" in i/j = " + repr([i,j]) + " indices = " + repr([reac.educt[i],reac.educt[j]]) + " adding " + repr([reac.e_stoch[j], dfdn]))
                
                reaction_matrix(ei, ej) += reaction.e_stoch[ej] * dfdn;           //Change e_stoch j --> ej!
                reaction_b(ej)          -= reaction.e_stoch[ej] * dfdn_po;
            }
            
            //Now add alpha in column for products
            //for j,pro in enumerate(reac.product) {
            for(int& pj : reaction.products ) { 
                //cout<<"pc pos6, pj = "<<pj<<endl;
                
                reaction_matrix(ei, pj) -= reaction.p_stoch[pj] * dfdn;        //Same with p_stoch --> change j-> pj for memory access!
                reaction_b(pj)          += reaction.p_stoch[pj] * dfdn_po;
            }
        }
            
        //
        // Term t1 = F/dx (1-exp(-dtau))
        //
        double t1    = 0.;
        //for i,educt in enumerate(reac.educt) { //For multibands different kappa, F
        for(int& ei : reaction.educts) {
            double dtau    = ds*kappa*n_olds[ei];
            t1     += dt*F/ds*(-std::expm1(-dtau))/n_tot;
        }
        
        for(int& ei : reaction.educts) {
        //for i,educt in enumerate(reac.educt) {
            reaction_b(ei)              -= reaction.e_stoch[ei] * t1;           
        }
        
        for(int& pi : reaction.products) {
        //for i,product in enumerate(reac.product) {
            reaction_b(pi)            += reaction.p_stoch[pi] * t1;
        }
    }
    //cout<<"Intermediate results: "<<endl;
    //cout<<"matrix = "<<endl<<reaction_matrix<<endl<<" b ="<<reaction_b<<endl;
    //
    // Regular reactions
    //
    //for ri,reac in enumerate(reactions) {
    for(c_reaction& reaction : reactions) {
        //cout<<"pos4. reaction = "<<reaction.reaction_number<<endl;
        
        double reac_r = thermo_get_r(Tcell, reaction.r, reaction.delta_stoch);
        double eStochsum = 0.;
        for(int& ei : reaction.educts) {
            eStochsum += reaction.e_stoch[ei];
        }
        
        //for i,educt in enumerate(reac.educt) {
        for(int& ei : reaction.educts) {
            //cout<<" educt i "<<ei<<endl;
            double dfdn    = dt * reac_r * std::pow(n_tot, eStochsum-1.);
            double dfdn_po = 0.;
            //print("dt / reac.r / dfdn = " + repr([dt, reac.r, dfdn]))
            //
            // This loop constructs df/dn_d as (mathematical) product of all educts (A^k)^a * (B^k)^b * .... (d-1) * (D^k)^(d-1) ..... (Z^k)^z
            //
            for(int& ej : reaction.educts) {
                //cout<<" educt j "<<ej<<endl;
                //for j,eparter in enumerate(reac.educt):
                if(ei==ej) {
                    //print(" in i/j = " + repr([i,j]) + " e[j]/n[e[j]] = " + repr([reac.educt[j],n_olds[int(reac.educt[j])]] ) + " adding " + repr([(float(reac.e_stoch[j]))*n_olds[  reac.educt[j]  ]**(reac.e_stoch[ j ]-1.)]) + " g-1 = " +repr(reac.e_stoch[ j ]-1.))
                    dfdn *= (double(reaction.e_stoch[ej]))* std::pow(n_olds[  ej  ], (reaction.e_stoch[ ej ]-1.));
                }
                    
                else {
                    //print(" in i/j = " + repr([i,j]) + " e[j]/n[e[j]] = " + repr([reac.educt[j],n_olds[int(reac.educt[j])]] ) + " adding " + repr([n_olds[reac.educt[j]]**reac.e_stoch[j]]) + " n,g-1 = " + repr([n_olds[reac.educt[j]],reac.e_stoch[ j ]]))
                    dfdn *= std::pow(n_olds[ej], reaction.e_stoch[ej]);
                }
                    
            }
            //
            // Term t3
            //
            dfdn_po = -dfdn * n_olds[ei];  //df/dn|k * n^k
            //cout<<" dfdn_po "<<dfdn_po<<endl;
            //print("i/reac_educt[i]/dfdn/dfdn_po = " + repr([i,reac.educt[i], dfdn, dfdn_po]))
            //
            // This distributes t2 = dfdn*stoch throughout the matrix
            // and adds dfdn*n to b for this educt
            //
            //Now add alpha in column for educts
            
            for(int& ej : reaction.educts) {
                //cout<<"matrix educt j "<<ej<<endl;
            //for j,edu in enumerate(reac.educt):
                //print(" in i/j = " + repr([i,j]) + " indices = " + repr([reac.educt[i],reac.educt[j]]) + " adding " + repr([reac.e_stoch[j], dfdn]))
                reaction_matrix(ei, ej) += reaction.e_stoch[ej] * dfdn;
                reaction_b(ej)          -= reaction.e_stoch[ej] * dfdn_po;
            }
                
            //Now add alpha in column for products
            for(int& pj : reaction.products) {
                //cout<<"matrix product j "<<pj<<endl;
            //for j,pro in enumerate(reac.product):
                reaction_matrix(ei,pj) -= reaction.p_stoch[pj] * dfdn;
                reaction_b(pj)         += reaction.p_stoch[pj] * dfdn_po;
            }
            
        }
            
        //
        // Term t1 = k * (A^k)^a * (B^k)^b ...  as final act
        //
        double t1    = dt * reac_r * std::pow(n_tot, eStochsum-1.);
        //for i,educt in enumerate(reac.educt) {
        for(int& ei : reaction.educts) {
            t1 *= std::pow(n_olds[ei], reaction.e_stoch[ei]);
        }
        
        //for i,educt in enumerate(reac.educt)
        for(int& ei : reaction.educts) {
            reaction_b(ei)              -= reaction.e_stoch[ei] * t1;           
        }
        for(int& pi : reaction.products) {
        //for i,product in enumerate(reac.product)
            reaction_b(pi)            += reaction.p_stoch[pi] * t1 ;
        }
        
    }
    
    if(cdebug > 0) {
        cout<<"In solver_implicit_cchem3 matrix, = "<<endl;
        cout<<reaction_matrix<<endl;
        cout<<"b="<<endl;
        cout<<reaction_b<<endl;
    
    }
    
    LUchem.compute(identity_matrix + reaction_matrix.transpose()) ;
    n_news.noalias() = LUchem.solve(reaction_b);
    
    //reaction_matrix.diagonal().noalias() += dt * (friction_coefficients * unity_vector);
     //
     //   chem_vec_output.noalias() = LUchem.solve(b);
    
    for(int s=0;s<num_species; s++) {
        if(n_news(s) < 0) {
            cout<<"In solver_implicit_cchem3 cell, = "<<cell<<" matrix = "<<endl;
            cout<<reaction_matrix.transpose()<<endl;
            cout<<"b="<<endl;
            cout<<reaction_b<<endl;
            cout<<"nolds = "<<endl<<n_olds<<endl;
            cout<<"nnews = "<<endl<<n_news* n_tot<<endl;
            cout<<"NEGATIVE DENSITY IN SPECIES["<<s<<"] = " <<n_news(s)<<" at time = "<<globalTime<<endl;
            char aa;
            cin>>aa;
            
            return 1;
        }
        else {
            
            species[s].prim[cell].number_density = n_news(s) * n_tot;
            species[s].prim[cell].density            = species[s].prim[cell].number_density * species[s].mass_amu*amu;
            
        }
            
    }
    
    for (int s = 0; s < num_species; s++) {
                    //species[s].prim[cell].number_density = n_news(s);
        
                    //species[s].prim[j].speed = vX[s] ;
                    //species[s].prim[j].internal_energy = uX[s];
                    //species[s].prim[j].temperature = TX[s];
                }
    
    for (int s = 0; s < num_species; s++) {
                species[s].eos->update_eint_from_T(&species[s].prim[0], num_cells + 2);
                species[s].eos->update_p_from_eint(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_auxillary(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_conserved(&species[s].prim[0], &species[s].u[0], num_cells + 2);
            }
    
    /*
    cout<<"Intermediate results with dt= "<<dt<<endl;
    cout<<"matrix = "<<endl<<reaction_matrix.transpose()<<endl<<" b ="<<reaction_b<<endl;
    
    for(int si=0; si<num_species; si++)
        cout<<n_news(si)<<" ";
    
    char aa;
    cin>>aa;
    */
    //
    // Add conversions
    //
            
    return 0; //No need for returning x anymore, all things are done in place here
}


//
// Example specialized fast solver for one specific reaction (in this case: co + 3 h2 <-> h2o + ch4) with reaction coefficients kf (k_forward) and kr (k_reverse)
// Edit: too lazy to translate python code, do if necessary
//
int solver_cchem_implicit_specialized_cochem(double dt, int num_spec, int cdebug) {
    
    /*
    n_olds = np.array([H2_old, CO_old, H2O_old, CH4_old])
    double kf = 1.;
    double kr = kf*1e-2;
    
    a = get_total_atoms(n_olds)
    nh = a[0]
    nc = a[1]
    no = a[2]
    
    fk  = - kf * n_olds[1]*n_olds[0]**3. + kr * n_olds[2]*n_olds[3]
    d   = nh - 4.*nc - 2.*no + 6.* n_olds[1]
    dak = - kf * 1./8. * (d**3 + 18. * n_olds[1] * d**2)
    dbk = + kr * (2.*n_olds[1] - no - nc)
    
    denom = (1. - dt *( dak  + dbk))
    nom   = n_olds[1] + dt * ( fk - dak * n_olds[1] - dbk * n_olds[1])
    
    // Use atom conservation to get the other molecule number densities
    CO_new  = nom/denom
    CH4_new = nc - CO_new
    H2O_new = no - CO_new
    H2_new  = 0.5*nh - H2O_new - 2*CH4_new
    
    return np.array([H2_new, CO_new, H2O_new, CH4_new])*/
    
    return 0.;
}

//
// Computes the reaction rate for a reaction,  considers also the case that this might be a backward reaction
//
double thermo_get_r(double T, double r1, double delta_stoch) {
    
    return r1;
}


//
// Returns H0/RT, s0/R and cp/R according to the NASA polynomials
//
std::vector<double> thermo_poly(double T,double a1,double a2,double a3,double a4,double a5,double a6,double a7,double a8,double a9){
    
    double H0 = -a1*std::pow(T,-2) + a2 *std::log(T)/T + a3 + 0.5*a4*T + a5/3*T*T + a6/4*T*T*T + a7/5*T*T*T*T + a8/T;
    double s0 = -0.5*a1*std::pow(T,-2) - a2/T + a3 * std::log(T) + a4*T + 0.5*a5*T*T + a6/3.*T*T*T + a7/4.*T*T*T*T + a9;
    double cp = a1*std::pow(T,-2) + a2/T + a3 + a4*T + a5*T*T + a6*T*T*T + a7*T*T*T*T;

    std::vector<double> temp = {H0, s0, cp};
    return temp;
    
}
    

std::vector<double> get_thermo_variables(double T,string species_string) {
    //double h = 0., s = 0., c = 0.;
    std::vector<double> temp = {0., 0., 0.};
    //
    // Numbering of a-coefficients deliberately omits a[0] as nonexistent value, due to confusion. Data from Tsai et al. (2017) The VULCAN paper
    //
    if(species_string.compare("H")==0)
        if(T < 1000)
            temp = thermo_poly(T, 0.00000000E+00, 0.00000000E+00, 2.50000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 2.54737080E+04, -4.46682853E-01);
        else
            temp = thermo_poly(T, 6.07877425E+01, -1.81935442E-01, 2.50021182E+00, -1.22651286E-07, 3.73287633E-11, -5.68774456E-15, 3.41021020E-19, 2.54748640E+04, -4.48191777E-01);
    
   else if(species_string.compare("H2")==0)
        if(T < 1000)
            temp = thermo_poly(T, 4.07832321E+04, -8.00918604E+02, 8.21470201E+00, -1.26971446E-02, 1.75360508E-05, -1.20286027E-08, 3.36809349E-12, 2.68248466E+03, -3.04378884E+01);
        else
            temp = thermo_poly(T, 5.60812801E+05, -8.37150474E+02, 2.97536453E+00, 1.25224912E-03, -3.74071619E-07, 5.93662520E-11, -3.60699410E-15, 5.33982441E+03, -2.20277477E+00);
    
    else if(species_string.compare("O")==0)
        if(T < 1000)
            temp = thermo_poly(T, -7.95361130E+03, 1.60717779E+02, 1.96622644E+00, 1.01367031E-03, -1.11041542E-06, 6.51750750E-10, -1.58477925E-13, 2.84036244E+04, 8.40424182E+00);
        else
            temp = thermo_poly(T, 2.61902026E+05, -7.29872203E+02, 3.31717727E+00, -4.28133436E-04, 1.03610459E-07, -9.43830433E-12, 2.72503830E-16, 3.39242806E+04, -6.67958535E-01);
    
    else if(species_string.compare("OH")==0)
        if(T < 1000)
            temp = thermo_poly(T, -1.99885899E+03, 9.30013616E+01, 3.05085423E+00, 1.52952929E-03, -3.15789100E-06, 3.31544618E-09, -1.13876268E-12, 2.99121423E+03, 4.67411079E+00);
        else
            temp = thermo_poly(T, 1.01739338E+06, -2.50995728E+03, 5.11654786E+00, 1.30529993E-04, -8.28432226E-08, 2.00647594E-11, -1.55699366E-15, 2.01964021E+04, -1.10128234E+01);
    
    else if(species_string.compare("H2O")==0)
        if(T < 1000)
            temp = thermo_poly(T, -3.94796083E+04, 5.75573102E+02, 9.31782653E-01, 7.22271286E-03, -7.34255737E-06, 4.95504349E-09, -1.33693325E-12, -3.30397431E+04, 1.72420578E+01);
        else
            temp = thermo_poly(T, 1.03497210E+06, -2.41269856E+03, 4.64611078E+00, 2.29199831E-03, -6.83683048E-07, 9.42646893E-11, -4.82238053E-15, -1.38428651E+04, -7.97814851E+00);
    
    else if(species_string.compare("CH")==0)
        if(T < 1000)
            temp = thermo_poly(T, 2.22059013E+04, -3.40541153E+02, 5.53145229E+00, -5.79496426E-03, 7.96955488E-06, -4.46591159E-09, 9.59633832E-13, 7.24078327E+04, -9.10767305E+00);
        else
            temp = thermo_poly(T, 2.06076344E+06, -5.39620666E+03, 7.85629385E+00, -7.96590745E-04, 1.76430830E-07, -1.97638627E-11, 5.03042951E-16, 1.06223659E+05, -3.15475744E+01);
    
    else if(species_string.compare("C")==0)
        if(T < 1000)
            temp = thermo_poly(T, 6.49503147E+02, -9.64901086E-01, 2.50467548E+00, -1.28144803E-05, 1.98013365E-08, -1.60614403E-11, 5.31448341E-15, 8.54576311E+04, 4.74792429E+00);
        else
            temp = thermo_poly(T,-1.28913647E+05, 1.71952857E+02, 2.64604439E+00, -3.35306895E-04, 1.74209274E-07, -2.90281783E-11, 1.64218238E-15, 8.41059785E+04, 4.13004742E+00 );
            
    else if(species_string.compare("CH2")==0)
        if(T < 1000)
            temp = thermo_poly(T, 3.21892173E+04, -2.87760181E+02, 4.20358382E+00, 3.45540596E-03, -6.74619334E-06, 7.65457164E-09, -2.87032842E-12, 4.73362471E+04, -2.14362860E+00);
        else
            temp = thermo_poly(T, 2.55041803E+06, -7.97162539E+03, 1.22892449E+01, -1.69912292E-03, 2.99172860E-07, -2.76700749E-11, 1.05134174E-15, 9.64221689E+04, -6.09473991E+01);
            
    else if(species_string.compare("CH3")==0)
       if(T < 1000)
            temp = thermo_poly(T, -2.87618881E+04, 5.09326866E+02, 2.00214395E-01, 1.36360583E-02, -1.43398935E-05, 1.01355673E-08, -3.02733194E-12, 1.40827182E+04, 2.02277279E+01);
        else
            temp = thermo_poly(T, 2.76080266E+06, -9.33653117E+03, 1.48772961E+01, -1.43942977E-03, 2.44447795E-07, -2.22455578E-11, 8.39506576E-16, 7.48180948E+04, -7.91968240E+01);
            
    else if(species_string.compare("CH4")==0)
        if(T < 1000)
            temp = thermo_poly(T, -1.76685100E+05, 2.78618102E+03, -1.20257785E+01, 3.91761929E-02, -3.61905443E-05, 2.02685304E-08, -4.97670549E-12, -2.33131436E+04, 8.90432275E+01);
        else
            temp = thermo_poly(T, 3.73004276E+06, -1.38350148E+04, 2.04910709E+01, -1.96197476E-03, 4.72731304E-07, -3.72881469E-11, 1.62373721E-15, 7.53206691E+04, -1.21912489E+02);
        
    else if(species_string.compare("C2")==0)
       if(T < 1000)
            temp = thermo_poly(T, 5.55963451E+05, -9.98012644E+03, 6.68162037E+01, -1.74343272E-01, 2.44852305E-04, -1.70346758E-07, 4.68452773E-11, 1.44586963E+05, -3.44822970E+02);
        else
            temp = thermo_poly(T, -9.68926793E+05, 3.56109299E+03, -5.06413893E-01, 2.94515488E-03, -7.13944119E-07, 8.67065725E-11, -4.07690681E-15, 7.68179683E+04, 3.33998524E+01);
            
    else if(species_string.compare("C2H2")==0)
       if(T < 1000)
            temp = thermo_poly(T, 1.59811209E+05, -2.21664412E+03, 1.26570781E+01, -7.97965108E-03, 8.05499275E-06, -2.43330767E-09, -7.52923318E-14, 3.71261906E+04, -5.24433890E+01);
        else
            temp = thermo_poly(T, 1.71384741E+06, -5.92910666E+03, 1.23612794E+01, 1.31418699E-04, -1.36276443E-07, 2.71265579E-11, -1.30206620E-15, 6.26657897E+04, -5.81896059E+01);
            
    else if(species_string.compare("C2H3")==0)
       if(T < 1000)
            temp = thermo_poly(T, -3.34789687E+04, 1.06410410E+03, -6.40385706E+00, 3.93451548E-02, -4.76004609E-05, 3.17007135E-08, -8.63340643E-12, 3.03912265E+04, 5.80922618E+01);
        else
            temp = thermo_poly(T, 2.71808009E+06, -1.03095683E+04, 1.83657981E+01, -1.58013115E-03, 2.68059494E-07, -2.43900400E-11, 9.20909639E-16, 9.76505559E+04, -9.76008686E+01);
            
    else if(species_string.compare("C2H")==0)
        if(T < 1000)
            temp = thermo_poly(T, 1.34366949E+04, -5.06797072E+02, 7.77210741E+00, -6.51233982E-03, 1.03011785E-05, -5.88014767E-09, 1.22690186E-12, 6.89226999E+04, -1.87188163E+01);
        else
            temp = thermo_poly(T, 3.92233457E+06, -1.20475170E+04, 1.75617292E+01, -3.65544294E-03, 6.98768543E-07, -6.82516201E-11, 2.71926279E-15, 1.43326663E+05, -9.56163438E+01);
            
    else if(species_string.compare("C2H4")==0)
     if(T < 1000)
            temp = thermo_poly(T, -1.16360584E+05, 2.55485151E+03, -1.60974643E+01, 6.62577932E-02, -7.88508186E-05, 5.12522482E-08, -1.37034003E-11, -6.17619107E+03, 1.09333834E+02);
        else
            temp = thermo_poly(T, 3.40876367E+06, -1.37484790E+04, 2.36589807E+01, -2.42380442E-03, 4.43139566E-07, -4.35268339E-11, 1.77541063E-15, 8.82042938E+04, -1.37127811E+02 );
            
    else if(species_string.compare("C2H5")==0)
      if(T < 1000)
            temp = thermo_poly(T, -1.41131255E+05, 2.71428509E+03, -1.53497773E+01, 6.45167258E-02, -7.25914396E-05, 4.59911601E-08, -1.21836754E-11, 5.98141884E+02, 1.09096652E+02);
        else
            temp = thermo_poly(T, 4.16922040E+06, -1.66298214E+04, 2.79544213E+01, -3.05171576E-03, 5.68516004E-07, -5.68286360E-11, 2.35564856E-15, 1.13701009E+05, -1.63935800E+02);
            
    else if(species_string.compare("C2H6")==0)
        if(T < 1000)
            temp = thermo_poly(T, -1.86204416E+05, 3.40619186E+03, -1.95170509E+01, 7.56583559E-02, -8.20417322E-05, 5.06113580E-08, -1.31928199E-11, -2.70293289E+04, 1.29814050E+02);
        else
            temp = thermo_poly(T, 5.02578213E+06, -2.03302240E+04, 3.32255293E+01, -3.83670341E-03, 7.23840586E-07, -7.31918250E-11, 3.06546870E-15, 1.11596395E+05, -2.03941058E+02);
    
    else if(species_string.compare("C4H2")==0)
      if(T < 1000)
            temp = thermo_poly(T, 2.46754257E+05, -3.89785564E+03, 2.36608046E+01, -2.20807780E-02, 2.78110114E-05, -1.57734001E-08, 3.42316546E-12, 7.08690782E+04, -1.10917356E+02);
        else
            temp = thermo_poly(T, 2.32817991E+06, -8.92518609E+03, 2.11432688E+01, -1.36887128E-03, 2.32750316E-07, -2.12451762E-11, 8.05331302E-16, 1.05778842E+05, -1.08831357E+02);
    
    else if(species_string.compare("CO")==0)
        if(T < 1000)
            temp = thermo_poly(T, 1.48904533E+04, -2.92228594E+02, 5.72452717E+00, -8.17623503E-03, 1.45690347E-05, -1.08774630E-08, 3.02794183E-12, -1.30313188E+04, -7.85924135E+00);
        else
            temp = thermo_poly(T, 4.61919725E+05, -1.94470486E+03, 5.91671418E+00, -5.66428283E-04, 1.39881454E-07, -1.78768036E-11, 9.62093557E-16, -2.46626108E+03, -1.38741311E+01);
    
    else if(species_string.compare("CO2")==0)
      if(T < 1000)
            temp = thermo_poly(T, 4.94365054E+04, -6.26411601E+02, 5.30172524E+00, 2.50381382E-03, -2.12730873E-07, -7.68998878E-10, 2.84967780E-13, -4.52819846E+04, -7.04827944E+00);
        else
            temp = thermo_poly(T, 1.17696242E+05, -1.78879148E+03, 8.29152319E+00, -9.22315678E-05, 4.86367688E-09, -1.89105331E-12, 6.33003659E-16, -3.90835059E+04, -2.65266928E+01);
    
    
    else if(species_string.compare("CH2OH")==0)
        if(T < 1000)
            temp = thermo_poly(T, -1.56007624E+05, 2.68544628E+03, -1.34202242E+01, 5.75713947E-02, -7.28444999E-05, 4.83664886E-08, -1.29349260E-11, -1.59682041E+04, 9.96303370E+01);
        else
            temp = thermo_poly(T, 2.25034951E+06, -8.17318606E+03, 1.59963918E+01, -8.70413372E-04, 6.06918395E-08, 4.40834946E-12, -5.70230950E-16, 4.64531343E+04, -7.83515845E+01);
    
    
    else if(species_string.compare("H2CO")==0)
      if(T < 1000)
            temp = thermo_poly(T, -1.17391634E+05, 1.87362885E+03, -6.89028857E+00, 2.64156167E-02, -2.18638930E-05, 1.00569301E-08, -2.02347695E-12, -2.30735177E+04, 6.42042055E+01);
        else
            temp = thermo_poly(T, 1.70082541E+06, -7.62085384E+03, 1.47244755E+01, -1.64911175E-03, 3.29214472E-07, -3.49504977E-11, 1.52613500E-15, 3.14681295E+04, -7.38647850E+01);
    
    
    else if(species_string.compare("HCO")==0)
        if(T < 1000)
            temp = thermo_poly(T, -1.18985189E+04, 2.15153611E+02, 2.73022403E+00, 1.80651611E-03, 4.98430057E-06, -5.81456792E-09, 1.86968989E-12, 2.90575564E+03, 1.13677254E+01);
        else
            temp = thermo_poly(T, 6.94960612E+05, -3.65622338E+03, 9.60473117E+00, -1.11712928E-03, 2.87532802E-07, -3.62624774E-11, 1.80832960E-15, 2.54370444E+04, -3.58247372E+01);
    
    
    else if(species_string.compare("CH3O")==0)
        if(T < 1000)
            temp = thermo_poly(T, 8.65711766E+04, -6.63168525E+02, 2.25745567E+00, 2.26628379E-02, -2.97056640E-05, 2.19934135E-08, -6.58804338E-12, 4.17410213E+03, 8.17477790E+00);
        else
            temp = thermo_poly(T, 2.10118824E+06, -8.84196880E+03, 1.82264573E+01, -1.74348503E-03, 3.34043427E-07, -3.43067316E-11, 1.47389777E-15, 5.30958206E+04, -9.42250059E+01);
    
    else if(species_string.compare("CH3OH")==0)
       if(T < 1000)
            temp = thermo_poly(T, -2.41664289E+05, 4.03214719E+03, -2.04641544E+01, 6.90369807E-02, -7.59893269E-05, 4.59820836E-08, -1.15870674E-11, -4.43326117E+04, 1.40014219E+02);
        else
            temp = thermo_poly(T, 3.41157076E+06, -1.34550020E+04, 2.26140762E+01, -2.14102918E-03, 3.73005054E-07, -3.49884639E-11, 1.36607344E-15, 5.63608156E+04, -1.27781428E+02);
    
    else if(species_string.compare("CH3CO")==0)
        if(T < 1000)
            temp = thermo_poly(T, -7.19389413E+04, 1.46446517E+03, -6.63227613E+00, 4.10846838E-02, -4.22625664E-05, 2.48576682E-08, -6.29255848E-12, -9.30937081E+03, 6.42289762E+01);
        else
            temp = thermo_poly(T, 2.48538815E+06, -1.12071420E+04, 2.27752544E+01, -2.31426055E-03, 4.53618917E-07, -4.74263555E-11, 2.04466390E-15, 6.38008841E+04, -1.21535093E+02);
    
    else if(species_string.compare("O2")==0)
       if(T < 1000)
            temp = thermo_poly(T, -3.42556342E+04, 4.84700097E+02, 1.11901096E+00, 4.29388924E-03, -6.83630052E-07, -2.02337270E-09, 1.03904002E-12, -3.39145487E+03, 1.84969947E+01);
        else
            temp = thermo_poly(T, -1.03793902E+06, 2.34483028E+03, 1.81973204E+00, 1.26784758E-03, -2.18806799E-07, 2.05371957E-11, -8.19346705E-16, -1.68901093E+04, 1.73871651E+01);
    
    else if(species_string.compare("H2CCO")==0)
        if(T < 1000)
            temp = thermo_poly(T, 3.54959809E+04, -4.06306283E+02, 3.71892192E+00, 1.58350182E-02, -1.72619569E-05, 1.15737696E-08, -3.30584263E-12, -5.20999258E+03, 3.83960422E+00);
        else
            temp = thermo_poly(T, 2.01356492E+06, -8.20088746E+03, 1.75969407E+01, -1.46454452E-03, 2.69588697E-07, -2.66567484E-11, 1.09420452E-15, 4.17777688E+04, -8.72580358E+01);
    
    else if(species_string.compare("HCCO")==0)
     if(T < 1000)
            temp = thermo_poly(T, 6.95961270E+04, -1.16459440E+03, 9.45661626E+00, -2.33124063E-03, 5.16187360E-06, -3.52616997E-09, 8.59914323E-13, 2.53500399E+04, -2.72635535E+01);
        else
            temp = thermo_poly(T, 1.09392200E+06, -4.49822821E+03, 1.24644643E+01, -6.34331740E-04, 1.10854902E-07, -1.12548868E-11, 5.68915194E-16, 4.65228030E+04, -5.09907043E+01);
            
    return temp;
    
}

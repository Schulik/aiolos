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

extern double HOnly_cooling(const std::array<double, 3> nX, double Te);

double thermo_get_r(double T, double r1, double delta_stoch);
std::vector<double> get_thermo_variables(double T,string species_string);

//
//
// Define and push the wanted (photo)chemical reactions into the reaction arrays.
// This is for now how the user 
//
void c_Sim::init_reactions(int cdebug) {
    
    //std::vector<c_reaction> reactions;
    //std::vector<c_photochem_reaction> photoreactions;
    cout<<"Init photoreactions..."<<endl;
    
    double ns = num_species;
    
    // H+ H2+ ion system
    //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {0}, {1,4}, {1.}, {1.,1.}, 1., 13.6 )); //H + gamma -> H+ + e-
    photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {0}, {1,2}, {1.}, {1.,1.}, 1., 13.6 )); //H + gamma -> H+ + e-
    
    //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {2}, {3,4}, {1.}, {1.,1.}, 1./3., 13.6)); //H2 + gamma-> H2+ + e-
    //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {2}, {0}, {1.}, {2.}, 1./3., 13.6)); //H2 + gamma-> 2H
    //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {2}, {0,1,4}, {1.}, {1.,1.,1.}, 1./3., 13.6)); //H2 + gamma-> H + H+ + e-
    
    //photoreactions.push_back(c_photochem_reaction( ns, num_bands_in, 0, {3}, {0,1}, {1.}, {1.,1.}, 1., 13.6)); //H2+ + gamma-> H + H+
    
    cout<<"Init reactions..."<<endl;
    
    //reactions.push_back(c_reaction(0, ns, {0,2}, {2}, {2.,1.}, {2.}, 2.7e-31, -0.600, 0.));       //H+H+H2 -> H2+H2  From Gail&Sedlmayer Chapter 10.5
    //reactions.push_back(c_reaction(1, ns, {0,2}, {2}, {2.,1.}, {2.}, 2.7e-31, -0.600, 0.));       //H+H+H2 <- H2+H2
    
    //reactions.push_back(c_reaction(0, ns, {0}, {0,2}, {3.}, {1.,1.}, 2.3e-35)); //8.82e-33      //H+H+H -> H2 + H
    //reactions.push_back(c_reaction(1, ns, {0}, {0,2}, {3.}, {1.,1.}, 2.3e-35));      //H+H+H <- H2 + H
    
    reactions.push_back(c_reaction(0, ns, {1,2}, {0}, {1.,1.}, {1.}, 2.7e-13 )); //Electron-proton recombination
    //reactions.push_back(c_reaction(0, ns, {1,4}, {0}, {1.,1.}, {1.}, 2.7e-13 )); //Electron-proton recombination
    //reactions.push_back(c_reaction(0, ns, {3,4}, {2}, {1.,1.}, {1.}, 2.7e-13 )); //Electron-H2+ recombination
    
    
    int cnt = 0;
    for(c_reaction& reaction : reactions) {
            reaction.set_reac_number(cnt);
            reaction.set_base_pointer(this);
            cnt++;
    }
    for(c_photochem_reaction& reaction : photoreactions) {
            reaction.set_reac_number(cnt);
            reaction.set_base_pointer(this);
            cnt++;
    }
    
    //if(cdebug > 0) {
        
    cout<<" Reporting reactions in init chemistry... "<<endl;
    for(c_reaction& reaction : reactions) {
            cout<<" Reaction #"<<reaction.reaction_number<<": ";
            
            if(reaction.is_reverse_reac) {
                
                for(int& pi : reaction.products) {
                        cout<<int(reaction.p_stoch[pi]*1.01)<<" "<<species[pi].speciesname<<" ";
                        if(reaction.products.back() != pi) 
                            cout<<"+ ";
                    
                        
                }
                cout<<" <- ";
                for(int& ei : reaction.educts) {
                    cout<<int(reaction.e_stoch[ei]*1.01)<<" "<<species[ei].speciesname<<" ";
                    if(reaction.educts.back() != ei) 
                        cout<<"+ ";
                }
                
                cout<<" ...... dStoch = "<<reaction.delta_stoch<<" initial rrate(100.K, 3000.K) = "<<reaction.get_reaction_rate(100.)<<"/"<<reaction.get_reaction_rate(3000.)<<endl;
                
                
            }
            else {
                
                for(int& ei : reaction.educts) {
                    cout<<int(reaction.e_stoch[ei]*1.01)<<" "<<species[ei].speciesname<<" ";
                    
                    if(reaction.educts.back() != ei) 
                        cout<<"+ ";
                }
                cout<<" -> ";
                for(int& pi : reaction.products) {
                    cout<<int(reaction.p_stoch[pi]*1.01)<<" "<<species[pi].speciesname<<" ";
                    if(reaction.products.back() != pi) 
                        cout<<"+ ";
                }
                cout<<" ...... dStoch = "<<reaction.delta_stoch<<" initial rrate(100.K, 3000.K) = "<<reaction.get_reaction_rate(100.)<<"/"<<reaction.get_reaction_rate(3000.)<<endl;
                
            }
            
            
    }
    cout<<"Photoreaction: "<<endl;
    for(c_photochem_reaction& reaction : photoreactions) {
            cout<<" Reaction #"<<reaction.reaction_number<<": ";
            
            for(int& ei : reaction.educts) {
                cout<<int(reaction.e_stoch[ei]*1.01)<<" "<<species[ei].speciesname;
                
                        cout<<" + ";
            }
            cout<<photon_energies[reaction.band]/(ev_to_K * kb)<<" eV -> ";
            for(int& pi : reaction.products) {
                cout<<int(reaction.p_stoch[pi]*1.01)<<" "<<species[pi].speciesname;
                if(reaction.products.back() != pi) 
                        cout<<" + ";
            }
            
            if(photon_energies[reaction.band]/(ev_to_K * kb) < 1.) {
                cout<<" WARNING! Ionizing with photons of <1 eV energy! Current band photon energy = "<<photon_energies[reaction.band]/(ev_to_K * kb)<<endl;
                char stop;
                cin>>stop;
            }
            cout<<endl;
    }
    cout<<endl;
        
    //char a;
    //cin>>a;
        
    //}
    
    
}


//
// Wrapper function to make it possible to call different chemistry solvers (such as specialized one-equation solvers, that might be much faster than the general solver)
// and do sub-cycling in case negative densitites occur.
//
//
//
//
void c_Sim::do_chemistry() {
    
    for (int j = num_cells+1; j > 0; j--) {
        
        double temp_dt = dt;
        double t_div = 1.;
        
        //
        // Repeat levels
        //
        for(int i=1; i<=32; i*=2) {
            int check_chem;
            //
            // Individual repeats
            //
            for(int k = 0; k < i; k++) {
                
                double dt_eff = dt/((double)i);
                
                check_chem = solver_cchem_implicit_general(dt_eff, j, 0);
                //if(i>4)
                //    cout<<" chem cell "<<j<<" repeat level "<<i<<" repeat number "<<k<<" dt_div = "<<dt_eff<<" dt = "<<dt<<" species[0].n = "<<species[0].prim[j].number_density<<endl; 
                
            }
            if(check_chem == 0) break;
            
        }
        
        /*
        cout<<" finished chemistry in cell "<<j<<endl;
        char a;
        cin>>a;
        */
        
        for(int b=0; b<num_bands_in; b++) {
            update_tau_s_jb(j, b);
            //update_dS_jb(j, b);
            
            update_dS_jb_photochem(j, b);
        }
        
        //update_dS_jb
        
        
    }
    
    //if(chemistry == 2)
    //     solver_cchem_implicit_something_specialized();
}


void c_reaction::update_reaction_rate(double T) {
    
    this->r = reac_a * std::pow(T, reac_b) * std::exp(-reac_c/T);
    
    //TODO: If (is_reverse_reaction) r *= gibbs energy stuff;
    if(is_reverse_reac) {
        
        std::vector<double> educt = {0.,0.,0.};
        std::vector<double> product = {0.,0.,0.};
        std::vector<double> temp = {0.,0.,0.};
        for(int& ej : educts ) {
                
                temp = get_thermo_variables( T, base->species[ej].speciesname); //{H0, s0, cp};
                educt[0] += e_stoch[ej] * temp[0];
                educt[1] += e_stoch[ej] * temp[1];
                //cout<<" came back with dH = "<<temp[0]<<" dS = "<<temp[1]<<" stoch = "<<e_stoch[ej]<<" for educt species "<<base->species[ej].speciesname<<endl;
        }
        for(int& pj : products ) {
                temp = get_thermo_variables( T, base->species[pj].speciesname); //{H0, s0, cp};
                product[0] += p_stoch[pj] * temp[0];
                product[1] += p_stoch[pj] * temp[1];
                
                //cout<<" came back with dH = "<<temp[0]<<" dS = "<<temp[1]<<" stoch = "<<p_stoch[pj]<<" for product species "<<base->species[pj].speciesname<<endl;
        }
        
        double dH = product[0]-educt[0]; //products - reactants
        double dS = product[1]-educt[1]; //products - reactants
        current_dG = dH - T*dS;
        
        this->r *= pow(kb*T/1e6, delta_stoch) * std::exp(current_dG/(Rgas*T));
    } 
    
    //cout<<" Reaction "<<reaction_number<<" r = "<<r<<endl;
    //char a;
    //cin>>a;
}

double c_reaction::get_reaction_rate(double T) {
    
    update_reaction_rate(T);
    
    return r;
}

//
//
// Here we split all reactions into forward, backward and individual summation steps.
// Each term in a dn/dt equation is approximated via the semi-implicit Euler method. 
// k denotes the timestep, k+1 being the advanced one, A and B are number densities for species A and B (e.g. n_A, n_B etc.)
// where the chemical rate of change is the product of power-laws f(x^k) = f(A^k, B^k, ...) = = k_r * (A^k)^a * (B^k)^b * ....
// and the number density equation becomes with S-I-Euler:
//
//        1/s * dn/dt = f(x^k+1)  \approx f(x^k) + sum_reactants dfdn_e * (e^k+1 - e^k )  )
// 
// where we define labels 
//                                       t1     + t2_e1 +  t3_e1                         +t2_e2 + t3_e2 + ....
//
// and the stochiometric coefficient is s.
// 
// It is t1 = f(x^k), and f1 and all f3 terms need to appear in the b-vector element i, as b_i = stoch[i] (t1 + sum_e t3_e)
//                      whereas all t2-terms go on the left hand side of the equation, as they are multiplied by the advanced time educt e^k+1
//                      they will find their ways into the matrix.
// 
// The same is done for the photoionizing reactions, only f(x^k) changes, it is related to the ionisation rate per particle n_H as Gamma G = F_b*exp(-\tau) /(hnu * n_H * dx) * (1 - exp(-d\tau))
// and the photoionising rate equation, for example for the atomic hydrogen reaction H + hnu -> p+ + e- is then 
//
//           dn_H/dt = - n * G
//
// Recombination reactions are regular chemical reactions which don't conserve particle number (by using the appropriate stochiometric coefficients).
//
int c_Sim::solver_cchem_implicit_general(double dtt, int cell, int cdebug) {
    
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
    }
    //
    // Photochemical reactions
    //
    //for ri,reac in enumerate(photochem_reactions) {
    
    //
    // For all bands, check participant photoreactions, construct ionization rates and generate chemical matrix entries
    //
    for(int b=0; b<num_bands_in; b++) {
        
        //int b = reaction.band;
        double temptau = 0.;
            
        double dlognu = 1.;
        double ds = dx[cell];
            
        if(cell < num_cells)
            temptau = radial_optical_depth_twotemp(cell+1,b);
            
        double F = 0.25 * solar_heating(b) / photon_energies[b] * std::exp(-temptau) * dlognu;
        
        double ntot_b = 0.;
        double tau_tot_b = 0.;
    
        //
        // First run through photoreactions: Compute n_tot and dtau_tot at retarded time k. Take branching ratios into account for not-double counting.
        // CAUTION: This requires the user to correctly specify branching ratios, so that they add up to one 
        //
        for(c_photochem_reaction& reaction : photoreactions) {
            
            if(reaction.band >= b) {
                
                ntot_b    +=  reaction.branching_ratio * n_olds[reaction.educts[0]];
                tau_tot_b +=  reaction.branching_ratio * n_olds[reaction.educts[0]] * species[reaction.educts[0]].opacity_twotemp(cell, b) * species[reaction.educts[0]].mass_amu*amu ;
                
            }
            
        }
        ntot_b    *= n_tot;
        tau_tot_b *= n_tot*ds; 
        
        for(c_photochem_reaction& reaction : photoreactions) {
            
            if(reaction.band >= b) {
            
                //cout<<"pc pos3. Starting photoreaction = "<<reaction.reaction_number<<" b = "<<reaction.band<<" cell = "<<cell<<endl;
                
                //cout<<" F = "<<F<<" kappa = "<<kappa<<" tau = "<<temptau<<endl;
                //cout<<" dx = "<<ds<<" S = "<<solar_heating(b)<<" E = "<<photon_energies[b]<<endl;
                
                //
                // Term t1 = F/dx (1-exp(-dtau))
                //
                double t1    = 0.;
                
                for(int& ei : reaction.educts) { //We assue that photochem reactions look likes reaction.educts[ei] + hv -> products, but leave the more general case open here
                    
                    double dtau = 0.;
                    double dfdn = 0.;
                    double dfdn_po = 0.;
                    //cout<<"pc pos5, ei = "<<ei<<endl;
                    double kappa = species[ei].opacity_twotemp(cell, b) * species[ei].mass_amu*amu; 
                    
                    dtau    = ds*kappa*n_olds[ei]*n_tot;
                    //dfdn    = dtt*F*kappa / n_tot * reaction.branching_ratio * std::exp(-dtau);  //TODO: Replaces F, kappa, dx with the appropriate values for species and flux in cell 
                    dfdn    = dtt*F*kappa/ n_tot * reaction.branching_ratio;
                    dfdn   *= (dtau * std::exp(-tau_tot_b) - std::expm1(-tau_tot_b) * (1. - dtau/tau_tot_b)  ) / tau_tot_b;
                    
                    if(globalTime > 2.9e10 && cell == 280) {
                        double nh = species[ei].prim[cell].number_density;
                        double np = species[reaction.products[0]].prim[cell].number_density;
                        double ne = species[reaction.products[1]].prim[cell].number_density;
                        cout<<cell<<" F "<<F<<" dx = "<<ds<<" kappa "<<kappa<<" hnu/eV = "<<photon_energies[b]/(kb * ev_to_K)<<" F*kappa ="<<F*kappa<<endl;
                        cout<<" F/hnu = "<<0.25 * solar_heating(b) / photon_energies[b]<<" F/hnu/dx = "<<0.25 * solar_heating(b) / photon_energies[b]/ds<<endl;
                        cout<<" ntot = "<<(nh + np);
                        cout<<" other ntot = "<<(nh + ne);
                        cout<<" steady-steate fzero/fplus = "<<nh/(nh+np)<<"/"<<np/(nh+np)<<endl;
                        
                        double lhs = F*kappa;
                        double a = 2.7e-13/lhs;
                        double b = 1.;
                        double c = -(nh+np);
                        double nplus = (-1 + std::sqrt(1.-4.*a*c))/2./a;
                        double fplus = -nplus / c;
                        //double ac = 1./(a*c);
                        //double fplus = 0.5*(+1.*ac - std::sqrt(ac*(ac-4.)));
                        long double temp2 = 1.-4.*a*c;
                        //double fzero = 1.-(+1 - std::sqrt(1.-4.*a*c))/(2.*a*c); //(1.-fplus);
                        double fzero = 1.-(+1 - (double)sqrtl(temp2))/(2.*a*c); //(1.-fplus);
                        
                        double temp = -0.5 * (b + 1. * std::sqrt(b*b - 4*a*c));
                        double x1 = temp / a;
                        double x2 = - c / temp;
                        fzero = 1.+1./temp;
                        
                        cout<<" 1-4ac = "<<1.-4*a*c<<" sqrt(...) = "<<std::sqrt(1.-4.*a*c)<<" -b + sqrt(...) = "<<-1 + std::sqrt(1.-4.*a*c)<<" a = "<<a<<" a*c = "<<a*c<<endl;
                        cout<<" theory steady-steate fzero/fplus = "<<fzero<<"/"<<fplus<<" nplus = "<<nplus<<" other fzero = "<<(c+nplus)/c<<"/"<<1.+nplus/c<<endl;
                        char aaa;
                        cin>>aaa;
                    }
                    
                    dfdn_po = -dfdn*n_olds[ei];  //df/dn|k * n^k
                    
                    //
                    //
                    // This distributes t2 = dfdn*stoch throughout the matrix
                    // and adds dfdn*n to b for this educt
                    //
                    //
                    for(int& ej : reaction.educts ) {
                        
                        reaction_matrix(ei, ej) += reaction.e_stoch[ej] * dfdn;           //Change e_stoch j --> ej!
                        reaction_b(ej)          -= reaction.e_stoch[ej] * dfdn_po;
                    }
                    
                    //Now add alpha in column for products
                    //for j,pro in enumerate(reac.product) {
                    for(int& pj : reaction.products ) { 
                        
                        reaction_matrix(ei, pj) -= reaction.p_stoch[pj] * dfdn;        //Same with p_stoch --> change j-> pj for memory access!
                        reaction_b(pj)          += reaction.p_stoch[pj] * dfdn_po;
                    }
                    
                    //t1     += dtt*F/ds/n_tot * reaction.branching_ratio*(-std::expm1(-dtau));
                    t1     += dtt*F/ds/n_tot * reaction.branching_ratio * (-std::expm1(-tau_tot_b)) * dtau/tau_tot_b;
                }
                
                
                //
                // Write the RHS of the Matrix equation
                //
                for(int& ei : reaction.educts) {
                //for i,educt in enumerate(reac.educt) {
                    reaction_b(ei)              -= reaction.e_stoch[ei] * t1;           
                }
                
                for(int& pi : reaction.products) {
                //for i,product in enumerate(reac.product) {
                    reaction_b(pi)            += reaction.p_stoch[pi] * t1;
                }
            
            }
        }
    
    }
    
    //cout<<"matrix = "<<endl<<reaction_matrix<<endl<<" b ="<<reaction_b<<endl;
    //
    // Regular reactions
    //
    //for ri,reac in enumerate(reactions) {
    double total_chem_dG = 0.;
    
    for(c_reaction& reaction : reactions) {
        double meanT = 0.;
        double denom = 0.;
        for(int& ei : reaction.educts) {
            meanT += species[ei].prim[cell].temperature * species[ei].prim[cell].density;
            denom += species[ei].prim[cell].density;
        }
        meanT = std::sqrt(meanT/denom);
        
        double reac_r = reaction.get_reaction_rate(meanT);
        total_chem_dG += reaction.current_dG;
        
        /*if(globalTime > 1e3) {
                
                cout<<"reaction "<<reaction.reaction_number<<" reac_r "<<reac_r<<" meanT "<<meanT<<endl;
                
                char a;
                cin>>a;
            }*/
        
        //cout<<" reaction number "<<reaction.reaction_number<<" r = "<<reaction.r<<" a/b/c = "<<reaction.reac_a<<" / "<<std::pow(meanT, reaction.reac_b)<<" / "<<std::exp(-meanT/reaction.reac_c)<<" c = "<<reaction.reac_c<<endl;
        
        double eStochsum = 0.;
        for(int& ei : reaction.educts) {
            eStochsum += reaction.e_stoch[ei];
        }
        
        for(int& ei : reaction.educts) {
            
            double dfdn    = dtt * reac_r * std::pow(n_tot, eStochsum-1.);
            double dfdn_po = 0.;
            
            //
            // This loop constructs df/dn_d as product of all educts (A^k)^a * (B^k)^b * .... (d-1) * (D^k)^(d-1) ..... (Z^k)^z
            //
            for(int& ej : reaction.educts) {
                if(ei==ej) {
                    dfdn *= (double(reaction.e_stoch[ej]))* std::pow(n_olds[  ej  ], (reaction.e_stoch[ ej ]-1.));
                }
                    
                else {
                    dfdn *= std::pow(n_olds[ej], reaction.e_stoch[ej]);
                }
                    
            }
            //
            // Term t3
            //
            dfdn_po = -dfdn * n_olds[ei];  //df/dn|k * n^k
            //cout<<" dfdn_po "<<dfdn_po<<endl;
            //
            // This distributes t2 = dfdn*stoch throughout the matrix
            // and adds dfdn*n to b for this educt
            //
            //Now add alpha in column for educts
            
            for(int& ej : reaction.educts) {
                //cout<<"matrix educt j "<<ej<<endl;
                reaction_matrix(ei, ej) += reaction.e_stoch[ej] * dfdn;
                reaction_b(ej)          -= reaction.e_stoch[ej] * dfdn_po;
            }
                
            //Now add alpha in column for products
            for(int& pj : reaction.products) {
                reaction_matrix(ei,pj) -= reaction.p_stoch[pj] * dfdn;
                reaction_b(pj)         += reaction.p_stoch[pj] * dfdn_po;
            }
            
        }
            
        //
        // Term t1 = k * (A^k)^a * (B^k)^b ...  as final act
        //
        double t1    = dtt * reac_r * std::pow(n_tot, eStochsum-1.);
        for(int& ei : reaction.educts) {
            t1 *= std::pow(n_olds[ei], reaction.e_stoch[ei]);
        }
        
        for(int& ei : reaction.educts) {
            reaction_b(ei)              -= reaction.e_stoch[ei] * t1;           
        }
        for(int& pi : reaction.products) {
            reaction_b(pi)            += reaction.p_stoch[pi] * t1 ;
        }
        
    }
    
    if(cdebug > 0) {
        cout<<"In solver_implicit_cchem3 matrix, = "<<endl;
        cout<<reaction_matrix<<endl;
        cout<<"b="<<endl;
        cout<<reaction_b<<endl;
    
        char a;
        cin>>a;
    }
    
    LUchem.compute(identity_matrix + reaction_matrix.transpose()) ;
    n_news.noalias() = LUchem.solve(reaction_b);
    
    //reaction_matrix.diagonal().noalias() += dt * (friction_coefficients * unity_vector);
     //
     //   chem_vec_output.noalias() = LUchem.solve(b);
    
    for(int s=0;s<num_species; s++) {
        if(n_news(s) < 0) {
            
            return 1;
            
            cout<<"NEGATIVE DENSITY IN SPECIES["<<s<<"] = " <<n_news(s)<<" at time = "<<globalTime<< " steps/cell = "<<steps<<"/"<<cell<<" dt_chem = "<<dtt<<endl;
            
            cout<<"In solver_implicit_cchem3 cell, = "<<cell<<" matrix = "<<endl;
            cout<<reaction_matrix.transpose()<<endl;
            cout<<"b="<<endl;
            cout<<reaction_b<<endl;
            cout<<"nolds = "<<endl<<n_olds<<endl;
            cout<<"nnews = "<<endl<<n_news* n_tot<<endl;
            cout<<"NEGATIVE DENSITY IN SPECIES["<<s<<"] = " <<n_news(s)<<" at time = "<<globalTime<< " steps/cell = "<<steps<<"/"<<cell<<endl;
            cout<<" This species density is being now floored. If you don't think this is a good idea, please ctrl+x."<<endl;
            
            //species[s].prim[cell].number_density = n_olds(s); //density_floor / species[s].mass_amu / amu;
            //species[s].prim[cell].density        = species[s].prim[cell].number_density * species[s].mass_amu*amu;
            //species[s].prim[cell].density        = density_floor * species[s].mass_amu;
            
            cout<<" Replaced new relative number density (comapre to nnews just above) = "<<species[s].prim[cell].number_density/n_tot<<endl;
            
            
            
            //char aa;
            //cin>>aa;
            
        }
        else {
            
            if(n_news(s) < 1e-30)
                n_news(s) = 1e-30;
            
            species[s].prim[cell].number_density = n_news(s) * n_tot;
            species[s].prim[cell].density        = species[s].prim[cell].number_density * species[s].mass_amu*amu;
            //species[s].prim[cell].internal_energy *= n_news(s)/n_olds[s];
        }
            
    }
    
    //Correct species internal energy and speed due to change in number density
    //Correct species internal energy due to total change in Gibbs energy dU = dG - PdV + TdS
    for (int s = 0; s < num_species; s++) {
                    //species[s].prim[cell].number_density = n_news(s);
        
                    //species[s].prim[j].speed = vX[s] ;
                    //species[s].prim[j].internal_energy = uX[s];
                    //species[s].prim[j].temperature = TX[s];
    }
    
    for (int s = 0; s < num_species; s++) {
                //species[s].eos->update_eint_from_T(&species[s].prim[cell], 1);
                species[s].eos->update_p_from_eint(&species[s].prim[cell], 1);
                species[s].eos->compute_auxillary(&species[s].prim[cell], 1);
                species[s].eos->compute_conserved(&species[s].prim[cell], &species[s].u[cell], 1);
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
//
// Take the already solved chemical and photochemical equations and update dS(j,s) for each species according to the ionizing photons absorbed
//
//
void c_Sim::update_dS_jb_photochem(int cell, int b) {
    
    n_olds          = Vector_t(num_species);
    double n_tot = 0.;
    
    for(int s=0;s<num_species; s++) {
        n_tot += species[s].prim[cell].number_density;
    }
    
    for(int s=0;s<num_species; s++) {
        n_olds[s]     = species[s].prim[cell].number_density / n_tot;
        reaction_b[s] = n_olds[s];
    }
    
    for (int b=0; b<num_bands_in; b++) {
        
        cell_optical_depth_highenergy(cell, b) = 0.;
        
        double temptau = 0.;
                
        double dlognu = 1.;
        double ds = dx[cell];
                
        if(cell < num_cells)
            temptau = radial_optical_depth_twotemp(cell+1,b);
                
        double F = 0.25 * solar_heating(b) * std::exp(-temptau) * dlognu;
            
        double ntot_b = 0.;
        double tau_tot_b = 0.;
        
        //
        // First run through photoreactions: Compute n_tot and dtau_tot at retarded time k. Take branching ratios into account for not-double counting.
        // CAUTION: This requires the user to correctly specify branching ratios, so that they add up to one 
        //
        for(c_photochem_reaction& reaction : photoreactions) {
            
            if(reaction.band >= b) {
                
                ntot_b    +=  reaction.branching_ratio * n_olds[reaction.educts[0]];
                tau_tot_b +=  reaction.branching_ratio * n_olds[reaction.educts[0]] * species[reaction.educts[0]].opacity_twotemp(cell, b) * species[reaction.educts[0]].mass_amu*amu ;
                
            }
            
        }
        ntot_b    *= n_tot;
        tau_tot_b *= n_tot*ds;
        cell_optical_depth_highenergy(cell, b) = tau_tot_b;
        
        double dS = -std::expm1(-tau_tot_b) * F / ds;
    
        for(c_photochem_reaction& reaction : photoreactions) {
            
            if(reaction.band >= b) {
            
                //Does this work???
                double tau_i =  reaction.branching_ratio * n_olds[reaction.educts[0]] * species[reaction.educts[0]].opacity_twotemp(cell, b) * species[reaction.educts[0]].mass_amu*amu ;
                
                //Get fraction of total cell-heating
                double fractional_dS = tau_i/tau_tot_b * dS * (1.- reaction.threshold_energy/photon_energies[b]);
                
                //Distribute energy according to mass
                for(int& pj : reaction.products) {
                    if(reaction.count_p == 1)
                        species[pj].dS(cell) += fractional_dS;
                    else {
                        species[pj].dS(cell) += fractional_dS * reaction.energy_split_factor/species[pj].mass_amu;
                        //cout<<" In highenergy heating cell "<<cell<<": target = "<<species[pj].speciesname<<" recieving a fraction = "<< reaction.energy_split_factor/species[pj].mass_amu<<" of the fractional "<<fractional_dS<<" heating. tau_tot_b = "<<tau_tot_b<<" expm-dtau = "<<-std::expm1(-tau_tot_b)<<endl;
                        
                    }
                        
                }
                
            }
        }
    }
    
    
    //
    // Do highenergy-cooling
    //
    int e_idx     = get_species_index("e-");
    int hnull_idx = get_species_index("H0"); 
    int hplus_idx = get_species_index("H+");
    
    if( e_idx!=-1 && hnull_idx!=-1 && hplus_idx!=-1  ) {
        
        std::array<double, 3> nX = {
                    species[hnull_idx].prim[cell].number_density,
                    species[hplus_idx].prim[cell].number_density,
                    species[e_idx].prim[cell].number_density};
    
        double cooling = nX[2] * HOnly_cooling(nX, species[e_idx].prim[cell].temperature);
        species[e_idx].dG(cell) += cooling;
        
        cout<<" We indeed found all cooling species, and cooling is = "<<cooling;
    }
    
    
    
    
    //
    // Do momentum correction
    //
    
    
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
    if(species_string.compare("H0")==0)
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

///////////////////////////////////////////////////////////
//
//
//  photochem.cpp
//
//
//
//
//
//
///////////////////////////////////////////////////////////

#define EIGEN_RUNTIME_NO_MALLOC

#include <array>
#include <cassert>

#include "aiolos.h"
#include "brent.h"

/*

This requires 8-12 species that we track, besides neutral H/He.
The couplings between high-energy radiation, ionization and heating happens in
this routine.

For this we assume a fixed ordering and numbering of species, as follows:
0   1   2   3   4    5    6       7    8     9    10    11
H0  H+  e-  He  He+  H++  He 23S  H2   H2+   H3+  HeH+  H-

Heating and cooling rates according to Black (1981), the physical state of
primordial intergalactic clouds*/

// Atomic data for Hydrogen:
double H_radiative_recombination(double T_e) {
    using std::pow;
    if (T_e < 220) T_e = 220.;
    
    double x = 2 * 157807 / T_e;
    // Case B
    //return 2.753e-14 * pow(x, 1.500) / pow(1 + pow(x / 2.740, 0.407), 2.242);
    ////       2.753e-14 * pow(2 * 157807 / x, 1.500) / pow(1 + pow(2 * 157807 / x / 2.740, 0.407), 2.242)
    return 2.7e-13 * pow(T_e/1e4, -0.9); //MC2009 value
}
double H_threebody_recombination(double T_e) {
    using std::pow;
    
    if (T_e < 220) T_e = 220.;

    double x = 2 * 157807 / T_e;
    return 1. * (1.005e-14 / (T_e * T_e * T_e)) * pow(x, -1.089) /
           pow(1 + pow(x / 0.354, 0.874), 1.101);
}
double H_collisional_ionization(double T_e) {
    using std::exp;
    using std::pow;

    if (T_e < 220) T_e = 220.; //return 1e-50; //T_e = 220.;//return 0 ;

    double x = 2 * 157807. / T_e;
    double term = 21.11 * pow(T_e, -1.5) * pow(x, -1.089) /
                  pow(1 + pow(x / 0.354, 0.874), 1.101);
    return 1e-50* term * exp(-0.5 * x);
}

// Cooling rate per electron.
double HOnly_cooling(const std::array<double, 3> nX, double Te) {
    using std::exp;
    using std::pow;
    using std::sqrt;

    const double T_HI = 157807.;
    if (Te < 220) Te = 220.;
    
    double x = 2 * T_HI / Te;
    double cooling = 0;

    // Recombination Cooling:
    double term ;
    if (x < 1e5) 
        //term = 3.435e-30 * Te * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);
        term = 3.435e-30 * Te * x*x / pow(1 + pow(x / 2.250, 0.376), 3.720);
    else 
        term = 7.562-27 * std::pow(Te, 0.42872) ;
    cooling += nX[1] * term;
    

    // Collisional ionization cooling:
    term = kb * T_HI * H_collisional_ionization(Te);
    cooling += nX[0] * term;

    // HI Line cooling (Lyman alpha):
    //term = 7.5e-19 * exp(-0.75 * T_HI / Te) / (1 + sqrt(Te / 1e5));
    term = 7.5e-19 * exp(-0.75 * T_HI / Te); //MC2009
    cooling += nX[0] * term;

    // Free-Free:
    term = 1.426e-27 * 1.3 * sqrt(Te) ;
    cooling += nX[1] * term  ;

    return 1.*cooling;
}

/* class C2Ray_HOnly_ionization
 *
 * Helper class for updating ionization fraction in a Hydrogen only model.
 *
 * Uses the C2Ray approximation, i.e. it computes the new and average values,
 * x_e and <x_e>, of the electron fraction, n_e/n_H, by solving
 *    dx_p/dt = Gamma0*exp(-tau0(1-x_p)) - R*nH*x_p*<x_e>
 *        + C*nH*<x_e>*(1-x_p) - B*nH*nH*<x_e>*<x_e>*x_p,
 * and assuming that <x_e> = <x_p>.
 *
 * Note that this method can be extended to support multiple species by adding
 * their contribution to <x_e>.
 */ 
class C2Ray_HOnly_ionization {
   public:
    C2Ray_HOnly_ionization(double* Gamma0_,double* tau0_, double dt_,
                           const std::array<double, 3>& nX, double Te, int num_he_bands_)
     : Gamma0(Gamma0_),
       tau0(tau0_),
       dt(dt_),
       nH(nX[0] + nX[1]),
       ne0(nX[2]),
       x0(nX[1] / (nX[0] + nX[1])),
       R(H_radiative_recombination(Te)),
       C(H_collisional_ionization(Te)),
       B(H_threebody_recombination(Te)),
       num_he_bands(num_he_bands_)
    {
        C = std::max(C, 1e-15*R) ;
    };

    // Computes <x_e> - xe, the residual in the average electron abundance
    double operator()(double xe) const { return _xe_average(xe) - xe; }

   private:
    // Compute the average electron abundance, given a guess for the average
    // electron abundance.
    double _xe_average(double x) const {
        double ne = ne0 + (x - x0) * nH;

        // Compute equilibrium
        double ion = photoionization_rate(x) ;

        double x_eq = (ion + ne * C) / (ion + ne * (C + R + ne * B));
        double t_rat = dt * (ion + ne * (C + R + ne * B));

        if (t_rat > 1e-10)
            return x_eq - (x0 - x_eq) * std::expm1(-t_rat) / t_rat;
        else
            return x_eq + (x0 - x_eq) * (1 - 0.5 * t_rat * (1 - t_rat / 3));
    }

   public:
    // Compute the average abundance fractions over the time-step, given an
    // electron fraction
    std::array<double, 3> nX_average(double xe) const {
        double xbar = _xe_average(xe);
        return {nH * (1 - xbar), nH * xbar, nH * (xbar - x0) + ne0};
    }

    // Compute the new abundance abundance fractions, given an electron
    // fraction.
    std::array<double, 3> nX_new(double x_bar) const {
        double ne = ne0 + (x_bar - x0) * nH;

        // Compute equilibrium
        double ion = photoionization_rate(x_bar) ;
        
        double x_eq = (ion + ne * C) / (ion + ne * (C + R + ne * B));
        double t_rat = dt * (ion + ne * (C + R + ne * B));

        double x_new = x_eq + (x0 - x_eq) * std::exp(-t_rat);

        //cout<<" in nX_newm "<<Gamma0[0]<<" = Gamma0, taub = "<<tau0[0]<<", ion = "<<ion<<" x_eq "<<x_eq<<" t_rat = "<<t_rat<<endl;
        //cout<<" dt = "<<dt<<" ion = "<<ion<<"  ne * (C + R + ne * B) = "<< ne * (C + R + ne * B)<<" ne = "<<ne<<endl;
        
        //return {nH * (1 - x_new), nH * x_new, nH * (x_new - x0) + ne0};
        return {nH * (1 - x_new), nH * x_new, nH * x_new};
    }

    double photoionization_rate(double x_bar) const {
        double ion = 0;
        for(int b=0; b < num_he_bands; b++)
            ion += Gamma0[b]/nH * -std::expm1(-tau0[b] * (1 - x_bar));
        return ion ;
    }

    double *Gamma0, *tau0, dt, nH, ne0, x0, R, C, B, num_he_bands;
};

/* class C2Ray_HOnly_heating
 *
 * Helper class for updating the temperature after ionization in a H-only
 * model.
 *
 * The goal of this class is to provide the net heating rate (photoheating +
 * line cooling). We compute the cooling rate implicitly to help reach
 * equilibrium, taking into account the collisional heat exchange between the
 * electrons, ions and neutrals.
 */
class C2Ray_HOnly_heating {
    using Mat3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
    using Vec3 = Eigen::Matrix<double, 3, 1>;

   public:
    using Matrix_t =
        Eigen::Matrix<double, NUM_SPECIES, NUM_SPECIES, Eigen::RowMajor>;

    C2Ray_HOnly_heating(double GammaH_, double dt_,
                        const std::array<double, 3>& nX,
                        const std::array<double, 3>& TX,
                        double ion_rate, double recomb_rate,
                        const Matrix_t& collisions)
     : GammaH(GammaH_),
       dt(dt_),
       ion(ion_rate), recomb(recomb_rate),
       nX(nX),
       T(TX),
       coll_mat(Mat3::Identity() - collisions * dt)
       {
           // Add in heat exchange due to ionization / recombination
           //   Here I've assumed m_e / m_p ~ 0
           coll_mat(0,0) += ion*dt/nX[0] ;
           coll_mat(1,0) -= ion*dt/nX[1] ;
           
           coll_mat(1,1) += recomb*dt/nX[1] ;
           coll_mat(0,1) -= recomb*dt/nX[0] ;           
       }; 

    // Compute the temperature residual.
    double operator()(double Te) const { return compute_T(Te)(2) - Te; }

    // Compute the cooling rate give the electron temperature.
    std::array<double, 3> net_heating_rate(double Te) const {
        Vec3 Tx = compute_T(Te) ;
        double heat_exch = 1.5*kb * (recomb*Tx(1) - ion*Tx(0)) ;

        return {heat_exch, -heat_exch, GammaH - nX[2]*HOnly_cooling(nX, Tx(2))};
    }
    
    std::array<double, 3> heating_rate(double /*Te*/) const {
        return {0.0, 0.0, GammaH};
    }
    
    std::array<double, 3> cooling_rate(double Te) const {
        Vec3 Tx = compute_T(Te) ;
        double heat_exch = 1.5*kb * (recomb*Tx(1) - ion*Tx(0)) ;

        return {heat_exch, -heat_exch, - nX[2]*HOnly_cooling(nX, Tx(2))};
    }
    
    Vec3 compute_T(double Te) const {
        // Setup RHS
        Vec3 RHS(T.data());
        RHS(2) += (GammaH/nX[2] - HOnly_cooling(nX, Te)) * dt / (1.5 * kb);

        // Solve for temperature
        return coll_mat.partialPivLu().solve(RHS);
    }

   private:
    

   public:
    double GammaH, dt, ion, recomb ;

    std::array<double, 3> nX, T;
    Mat3 coll_mat;
};

void c_Sim::do_photochemistry() {
    // cout<<" Doing photochemistry "<<endl;
            
            //double tau0 = 0;
        
            for (int b = 0; b < num_he_bands; b++) {
                //radial_optical_depth(num_cells + 1,b)         = 0.;
                radial_optical_depth_twotemp(num_cells + 1,b) = 0.;
            }
            
            for (int j = num_cells + 1; j > 0; j--) {
                
                for (int b = 0; b < num_he_bands; b++) {
                    cell_optical_depth_twotemp(j,b)           = 0.;
                    //total_opacity(j,b)                = 0 ;
                    total_opacity_twotemp(j,b)        = 0 ;
                }
                
                std::array<double, 3> nX = {
                    species[0].prim[j].number_density,
                    species[1].prim[j].number_density,
                    species[2].prim[j].number_density};

                std::array<double, 3> uX = {
                    species[0].prim[j].internal_energy,
                    species[1].prim[j].internal_energy,
                    species[2].prim[j].internal_energy};
                    
                std::array<double, 3> pX = {
                    species[0].prim[j].pres,
                    species[1].prim[j].pres,
                    species[2].prim[j].pres};
                    
                std::array<double, 3> vX = {
                    species[0].prim[j].speed,
                    species[1].prim[j].speed,
                    species[2].prim[j].speed};

                std::array<double, 3> TX = {
                    species[0].prim[j].temperature,
                    species[1].prim[j].temperature,
                    species[2].prim[j].temperature};

                std::array<double, 3> mX = {
                    species[0].mass_amu * amu,
                    species[1].mass_amu * amu,
                    species[2].mass_amu * amu};
                
                std::vector<double> dtaus   = np_zeros(num_he_bands);
                std::vector<double> Gamma0  = np_zeros(num_he_bands);

                // Remove negative densities/temperatures
                for (int s=0; s < 3; s++) {
                    nX[s] = std::max(nX[s],0.0);
                    TX[s] = std::max(TX[s],1.0);
                }

                // The local optical depth to ionizing photons
                // under the "initially all neutral" assumption

                // Assume 20eV photons for now
                //double E_phot = 20 * ev_to_K * kb;
                
                //std::array<double, num_he_bands> Gamma0;
                double dlognu = 1.; //1.8833e1;
                
                for (int b = 0; b < num_he_bands; b++) {
                    Gamma0[b] = 0.25 * solar_heating(b) / photon_energies[b] * std::exp(-radial_optical_depth_twotemp(j+1,b)) / dx[j] * dlognu;
                    dtaus[b]  = species[0].opacity_twotemp(j, b) * (nX[0] + nX[1]) * mX[0] * dx[j]; //TODO: Replace this approx with the actual HE optd values, but get everything else to run properly first
                    
                    if(debug >= 0 && j >= 1053 && steps <= 2) {
                        cout<<" PH j= "<<j<<" b == "<<b<<" Gamma0[b] = "<<Gamma0[b]<<" dtau = "<<dtaus[b]<<" opa/nX0+nX1/dx[j] = "<<species[0].opacity_twotemp(j, b)<<"/"<<(nX[0] + nX[1])<<"/"<<dx[j]<<endl;
                    }
                }
                
                if(debug >= 0 && j >= 52 && j <= 53 && steps >= 8955 && steps <= 8965) {
                        //cout<<" j = "<<j<<" G0/G1 before "<<Gamma0[0]<<"/"<<Gamma0[1];
                }
                
                // Update the ionization fractions
                C2Ray_HOnly_ionization ion(&Gamma0[0], &dtaus[0], dt, nX, TX[2], num_he_bands);
                
                if(debug >= 0 && j >= debug && j <= debug_cell && steps >= debug_steps) {
                        //cout<<" j = "<<j<<" G0/G1 after "<<Gamma0[0]<<"/"<<Gamma0[1]<<endl;
                }
                
                Brent root_solver_ion(ion_precision, ion_maxiter);
                double xmin = std::max((nX[1] - nX[2]) / (nX[0] + nX[1]), 0.0);
                double x_bar;
                try {
                        x_bar = root_solver_ion.solve(xmin, 1, ion);
                }
                catch(...) {
                    
                    
                        double Tn = species[0].prim[j].temperature;
                        double Tp = species[1].prim[j].temperature;
                        double Te = species[2].prim[j].temperature;
                        double Gamma0final = 0.;
        
                        for(int b=1; b < num_he_bands; b++)
                            Gamma0final += Gamma0[b]/nX[0] * -std::expm1(-dtaus[b] * (1 - x_bar));
                        
                        double R = H_radiative_recombination(Te);
                        double C = H_collisional_ionization(Te);
                        double B = H_threebody_recombination(Te);
                        double ne = nX[2];
                        double x_eq = (Gamma0final + ne * C) / (Gamma0final + ne * (C + R + ne * B));
                        double t_rat = dt * (Gamma0final + ne * (C + R + ne * B));
                        double x1 = species[1].u[j].u1  /(species[0].u[j+1].u1+species[1].u[j+1].u1);
                        double x2 = species[1].u[j+1].u1/(species[0].u[j+1].u1+species[1].u[j+1].u1);
                        double x_grad = std::abs( x1/x2 );
                        
                        if(globalTime < -1e5) {
                            cout<<" ERROR IN ROOT steps = "<<steps<<" j= "<<j<<" b == "<<1<<" G0,G1 = "<<Gamma0[0]<<"/"<<Gamma0[1]<<" taus0/1 = "<<dtaus[0]<<"/"<<dtaus[1]<<"/"<<" opa/xbar/xeq = "<<species[0].opacity_twotemp(j, 1)<<"/"<<x_bar<<"/"<<x_eq<<" dS = "<<dS_band(j,1)<<" SBand = "<<S_band(j,1)<<" t_rat ="<<t_rat<<endl;
                            cout<<" ERROR: time/G0 = "<<globalTime<<"/"<<Gamma0final<<" R/C/B = "<<R<<"/"<<C<<"/"<<B<<" Tn/Tp/Te ="<<Tn<<"/"<<Tp<<"/"<<Te<<" x_grad "<<x_grad<<" dtau = "<<dtaus[1]<<endl;
                            std::cout << j << " " <<"\n"
                                << "\t" << TX[0] << " " << TX[1] << " " << TX[2] << "\n"
                                << "\t" << uX[0] << " " << uX[1] << " " << uX[2] << "\n" 
                                << "\t" << nX[0] << " " << nX[1] << " " << nX[2] << "\n" 
                                << "\t" << pX[0] << " " << pX[1] << " " << pX[2] << "\n" ;;
                        
                        
                            
                        }
                        
                }
                

                std::array<double, 3> nX_new = ion.nX_new(x_bar);
                std::array<double, 3> nX_bar = ion.nX_average(x_bar);

                // Store the optical depth through this cell:
                for (int b = 0; b < num_he_bands; b++)
                    dtaus[b] *= (1 - x_bar);
                
                //tau0 += dtau;

                // Next update the primitive quantities to ensure momentum / energy
                // conservation 

                // Momentum conservation:
                double ne = nX_bar[2], nH = nX_bar[0] ;
                double dn_R = (ion.R + ion.B*ne)*ne*ne*dt ;
                double dn_I = (ion.C*ne + ion.photoionization_rate(x_bar))*nH*dt ;

                std::array<double, 3> mom = { 
                    nX[0]*mX[0]*vX[0], nX[1]*mX[1]*vX[1], nX[2]*mX[2]*vX[2] } ;

                double fe = 1/(1 + mX[2]/mX[1]), fp = 1/(1 + mX[1]/mX[2]) ;
                double f = nX_new[0] + dn_I - dn_I*dn_R*(fp/(nX_new[1]+dn_R) + fe/(nX_new[2] + dn_R)) ;
                mom[0] = (mom[0] + dn_R*(mom[1]/(nX_new[1]+dn_R) + mom[2]/(nX_new[2]+dn_R))) / f ;
                mom[1] = (mom[1] + dn_I*mom[0]*fp)/(nX_new[1]+dn_R) ;
                mom[2] = (mom[2] + dn_I*mom[0]*fe)/(nX_new[2]+dn_R) ;

                // Compute the change in kinetic energy:
                double dEk = 0.5*(mom[0]*mom[0]*nX_new[0] / mX[0] - mX[0]*nX[0]*vX[0]*vX[0] +
                                  mom[1]*mom[1]*nX_new[1] / mX[1] - mX[1]*nX[1]*vX[1]*vX[1] +
                                  mom[2]*mom[2]*nX_new[2] / mX[2] - mX[2]*nX[2]*vX[2]*vX[2]) ;
                
                // Update the velocity
                vX[0] = mom[0] / mX[0] ;
                vX[1] = mom[1] / mX[1] ;
                vX[2] = mom[2] / mX[2] ;

                // Next, compute the heating and cooling
                //    It is useful to re-scale the temperatures first to conserve the
                //    internal energy.
                //TX[0] *= nX[0] / nX_new[0] ;
                //TX[1] *= nX[1] / nX_new[1] ;
                //TX[2] *= nX[2] / nX_new[2] ;

                for (int s = 0; s < 3; s++) {
                    species[s].prim[j].number_density = nX_new[s];
                    species[s].prim[j].speed = vX[s] ;
                    species[s].prim[j].density = nX_new[s] * mX[s];
                    species[s].prim[j].internal_energy = uX[s];
                    species[s].prim[j].temperature = TX[s];
                }
                
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                double GammaH = 0.;
                double Te;
                std::array<double, 3> heating;
                std::array<double, 3> cooling;
                
                for (int b = 0; b < num_he_bands; b++) {
                    GammaH += 0.25 * solar_heating(b)  * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) / dx[j] * dlognu;
                }
                
                if(debug >= 2) {
                    
                    double dtaufinal = dtaus[1];
                    double ipp = Gamma0[1]/(nX[0]+nX[1]) * -std::expm1(-dtaufinal); 
                    double hpp = GammaH/(nX[0]+nX[1]);
                    
                    cout<<" j= "<<num_cells-3<<" steps == "<<steps<<" ion per particle = "<<ipp<<" heating per particle = "<<hpp<<endl;
                    
                }
                
                if(use_rad_fluxes) {
                        
                    Te = TX[2];
                    
                    heating = {0.0, 0.0, GammaH - dEk/(dt+1e-300)}; //    heat.heating_rate(Te);
                    cooling = {0., 0., - nX[2]*HOnly_cooling(nX, Te) }; //
                    
                    //
                    // We cannot cool more than what we got from highenergy photons
                    //
                    if(std::abs(cooling[2]) > heating[2])
                        cooling[2] = -heating[2];
                    
                    if(j > 39 && j <= 42 && steps <= 2) {

                        double Gamma0final = 0.;
        
                        for(int b=0; b < num_he_bands; b++)
                            Gamma0final += Gamma0[b]/nX[0] * -std::expm1(-dtaus[b] * (1 - x_bar));
                        
                        double R = H_radiative_recombination(Te);
                        double C = H_collisional_ionization(Te);
                        double B = H_threebody_recombination(Te);
                        double ne = nX[2];
                        cout<<" Cooling peak detected in photochem. DEBUG data:"<<endl;
                        cout<<" Te = "<<Te<<endl;
                        cout<<GammaH<<" "<<Gamma0[0]<<" "<<Gamma0[1]<<" "<<Gamma0final<<endl;
                        cout<<R<<" "<<C<<" "<<B<<endl;

                        
                        double x_eq = (Gamma0final + ne * C) / (Gamma0final + ne * (C + R + ne * B));
                        double t_rat = dt * (Gamma0final + ne * (C + R + ne * B));
                        //double x_new = x_eq + (x0 - x_eq) * std::exp(-t_rat);
                        
                        cout<<x_eq<<" "<<t_rat<<" "<<x_bar<<endl;
                         
                    }
                    
                } else {
                    
                    // Collisional heat exchange:
                    //if(!use_rad_fluxes) {
                    fill_alpha_basis_arrays(j);
                    compute_collisional_heat_exchange_matrix(j);
                    //}
                    
                    // Solve for radiative cooling implicitly
                    C2Ray_HOnly_heating heat(GammaH - dEk/(dt+1e-300), dt, nX_new, TX,
                           (ion.R + ion.B*ne)*ne*ne, (ion.C*ne + ion.photoionization_rate(x_bar))*nH,
                           friction_coefficients);

                    // Bracket the temperature:
                    double Te1 = TX[2], Te2;
                    if (heat(Te1) < 0) {
                        Te2 = Te1;
                        Te1 = Te2 / 1.4;
                        while (heat(Te1) < 0) {
                            Te2 = Te1 ;
                            Te1 /= 1.4;
                        }
                    } else {
                        Te2 = Te1 * 1.4;
                        while (heat(Te2) > 0) {
                            Te1 = Te2 ;
                            Te2 *= 1.4;
                        }
                    }

                    // Compute the new temperature and save the net heating rates
                    if (!(heat(Te1)*heat(Te2)<=0)) {
                        cout<<"if (!(heat(Te1)*heat(Te2)<=0)) {"<<endl;
                        std::cout << j << " " << Te1 << " " << Te2 << " " << TX[2] << "\n"  
                            << "\t" << heat(Te1) << " " << heat(Te2) << "\n" 
                            << "\t" << nX_bar[0] << " " << nX_bar[1] << " " << nX_bar[2] << "\n" 
                            << "\t" << TX[0] << " " << TX[1] << " " << TX[2] << "\n" 
                            << "\t\t" << heat.coll_mat(0,0) << " " <<  heat.coll_mat(0,1) << " " << heat.coll_mat(0,2) << "\n"
                            << "\t\t" << heat.coll_mat(1,0) << " " <<  heat.coll_mat(1,1) << " " << heat.coll_mat(1,2) << "\n"
                            << "\t\t" << heat.coll_mat(2,0) << " " <<  heat.coll_mat(2,1) << " " << heat.coll_mat(2,2) << "\n"
                            << "\t" << uX[0] << " " << uX[1] << " " << uX[2] << "\n" 
                              << "\t" << pX[0] << " " << pX[1] << " " << pX[2] << "\n" ;;
                    }

                    Brent root_solver_heat(ion_heating_precision, ion_heating_maxiter);
                    
                    Te = root_solver_heat.solve(Te1, Te2, heat);
                        
                    heating = heat.heating_rate(Te);
                    cooling = heat.cooling_rate(Te);
                    
                    Eigen::Matrix<double, 3, 1> newT = heat.compute_T(Te);
                    for (int s = 0; s < 3; s++)
                        species[s].prim[j].temperature = newT(s);   
                }
                
                if(debug >= 1 && j==num_cells) {
                    cout<<" End of photochem. Found xbar ="<<x_bar<<" in j ="<<j<<endl;
                }
                if(debug >= 2 ) {
                    cout<<" End of photochem. Found xbar ="<<x_bar<<" Gamma_H = "<<GammaH<<" in j ="<<j<<endl;
                }
                
                //
                // This loop documents quantities consistent with the ionization found, computes the highenergy heating residual that was not used to ionize and dumps it into
                // the species that do not participate in the C2-ray scheme, proportional to their optical depth per cell
                //
                for (int s = 0; s < 3; s++) {
                    //Even in cases when radiative transport is not used, we want to document those quanitites in the output
                    species[s].dG(j) += cooling[s];
                    species[s].dS(j) += heating[s];
                }
                double adddS[3] = {0.,0.,0.};
                for (int b = 0; b < num_he_bands; b++) {
                    
                    double temp_nonion_cell_optd = 0.;
                    double dtauplus = dtaus[b];
                    //double temp_total_cell_optd = 0.;
                    double GammaT = 0.;
                    double GammaHnofactor = 0.;
                    
                    //
                    // dtaus[b]  += species[0].opacity_twotemp(j, b) * (nX[0] + nX[1]) * mX[0] * dx[j];
                    //
                    
                    //cout<<"PCHEM pos 1"<<endl;
                    
                    for(int s=3; s<num_species; s++) {
                        temp_nonion_cell_optd   += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j];
                        dtauplus                += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j];
                    }
                    
                    //cout<<"PCHEM pos 2, temp_optd  = "<<temp_nonion_cell_optd<<endl;
                    //for(int s=0; s<num_species; s++) {
                    //    temp_total_cell_optd += species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j];
                    //}
                    for(int s=3; s<num_species; s++) {
                        species[s].fraction_total_solar_opacity(j,b) = species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / temp_nonion_cell_optd;
                    
                        //cout<<"PCHEM pos 3, s = "<<s<<species[s].fraction_total_solar_opacity(j,b)<<endl;
                    }
                    
                    if(j < num_cells + 1) {
                        GammaHnofactor = 0.25 * solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) / dx[j];
                        GammaT         = 0.25 * solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtauplus)) / dx[j];
                        radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + dtauplus;
                        S_band(j,b)                       = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)); 
                    }
                        
                    else {
                        GammaHnofactor = 0.25 * solar_heating(b) * dtaus[b];
                        GammaT         = 0.25 * solar_heating(b) * dtauplus;
                        radial_optical_depth_twotemp(j,b) = dtauplus;
                        S_band(j,b)                       = solar_heating(b);
                    }
                    
                    //cout<<"PCHEM pos 4, GammaH  = "<<GammaHnofactor<<endl;
                    //
                    // Now recompute heating with higher optical depth, neglecting ionization and ionization potential
                    //
                    
                    dS_band(j,b)  = (GammaT - GammaHnofactor);
                    
                    if(debug >= 2)
                        cout<<"PCHEM pos 6, dGamma  = "<<dS_band(j,b)<<endl;
                    //Distribute GammaResidual to other species according to their optical depth contribution
                    
                    for(int s=3; s<num_species; s++) {
                        //species[s].dS(j) += GammaResidual * species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / temp_total_cell_optd;
                        species[s].dS(j)    += dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                        adddS[s]            += dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                    }
                    
                    if(dS_band(j,b) > 1e5) {
                            for(int s=3; s<num_species; s++) {
                                    cout<<" PHOTOCHEM j/s/dS = "<<j<<"/"<<s<<"/"<<species[s].dS(j)<<endl;
                            }
                    } 
                    
                    //cout<<"PCHEM pos 7, radial_tau  = "<<radial_optical_depth_twotemp(j,b)<<endl;
                    
                    cell_optical_depth_twotemp(j,b)   = dtauplus;
                    total_opacity_twotemp(j,b)        = dtauplus/dx[j];
                    
                    //No factor 0.25 here, as we want to see the unaveraged radiation we put in, consistent with low-energy bands
                    //dS_band(j,b) = 0.25 * solar_heating(b) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) / dx[j];
                    
                    if( (steps == 2 || steps==1000 ) && j == num_cells-10) {
                        
                        double GammaHfinal = 0.;
                        double GammaHfinal2= 0;
                        double Gamma0final = 0.;
                        double one = 0.;
                        double one2 = 0.;
                        double nphot = 1.;
                        double npart = nX[0] + nX[1];
                        char wait;
                        
                        for(int b=0; b < num_he_bands; b++) {
                            
                            nphot = solar_heating(b) / photon_energies[b];
                            
                            //one  = std::exp(-radial_optical_depth_twotemp(j+1,b)) / photon_energies[b]
                            //one2 =
                            
                            GammaHfinal += 0.25 * solar_heating(b)  * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) / dx[j] * dlognu; 
                            
                            //GammaHfinal2 += 0.25 * solar_heating(b)  * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) * dlognu * (-std::expm1(-dtaus[b])) / dx[j] ; 
                            
                            GammaHfinal2 += (photon_energies[b] - 13.6 * ev_to_K * kb) * (1.-std::exp(-dtaus[b]) ) / dx[j] ;
                            
                            //Gamma0[b] = 0.25 * solar_heating(b) / photon_energies[b] * std::exp(-radial_optical_depth_twotemp(j+1,b)) / dx[j];
                            Gamma0final += Gamma0[b] * -std::expm1(-dtaus[b]);
                        }
                        
                        Gamma0final /= (nphot * npart);
                        GammaHfinal /= (nphot * npart);
                        GammaHfinal2 /= (nX[0] + nX[1]);
                        
                        cout<<" In photochem, x[j] = "<<x_i12[j]<<" G0/GH/GH2 = "<<Gamma0final<<" / "<<GammaHfinal<<"/"<<GammaHfinal2<<"/"<<GammaHfinal2*4.<<" nphot/npart = "<<nphot<<"/"<<npart<<" ne/nH , ne/(nH+ne), x_bar = "<<nX[2]/ nX[0]<<" / "<<nX[2]/(nX[0]+nX[1])<<"/"<<x_bar<<" photon energy = "<<photon_energies[1]/(ev_to_K * kb)<<endl;
                        
                        //cin>>wait;
                    }
                    
                    if(debug >= 0 && j >= debug_cell && j <= debug_cell && steps >= debug_steps && b==1 && steps <= debug_steps+10 ) {
                        //if(debug >= 0 && j >= 6 && j <= 6 && steps >= 361 && steps <= 371 && b==1) {
                        
                        double GammaHfinal = 0.;
                        double Gamma0final = 0.;
                        
                        for(int b=1; b < num_he_bands; b++) {
                            GammaHfinal += 0.25 * solar_heating(b)  * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b]) * x_bar) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) / dx[j];
                        }
                        
                        double R = H_radiative_recombination(Te);
                        double C = H_collisional_ionization(Te);
                        double B = H_threebody_recombination(Te);
                        double ne = nX[2];
                        double x_eq = (Gamma0final + ne * C) / (Gamma0final + ne * (C + R + ne * B));
                        double t_rat = dt * (Gamma0final + ne * (C + R + ne * B));
                        double x1 = species[1].u[j].u1 / (species[0].u[j+1].u1+species[1].u[j+1].u1);
                        double x2 = species[1].u[j+1].u1/(species[0].u[j+1].u1+species[1].u[j+1].u1);
                        double x_grad = std::abs( x1/x2 );
                        double T2 = species[2].prim[j+1].temperature;
                        double T0 = species[2].prim[j-1].temperature;
                        
                        //cout<<" steps/G0 = "<<steps<<"/"<<Gamma0final<<" R/C/B = "<<R<<"/"<<C<<"/"<<B<<" T0/Te/T2 ="<<T0<<"/"<<Te<<"/"<<T2<<" x_grad "<<x_grad<<" dtau/xbar/xeq/xmin = "<<dtaus[b]<<"/"<<x_bar<<"/"<<x_eq<<"/"<<xmin<<endl;
                        
                        cout<<" steps/GH = "<<steps<<"/"<<GammaH<<" R/B = "<<R<<"/"<<B<<" T0/Te/T2 ="<<T0<<"/"<<Te<<"/"<<T2<<" x_grad "<<x_grad<<" dS1 dG1 + = "<<species[2].dS(j)<<"/"<<species[2].dG(j)<<"/"<<species[2].dS(j)+species[2].dG(j)<<" dS2 dG2 + = "<<species[2].dS(j+2)<<"/"<<species[2].dG(j+2)<<"/"<<species[2].dS(j+2)+species[2].dG(j+2)<<endl;
                    }
                                
                }
                
            }

            // Finally, lets update the conserved quantities
            for (int s = 0; s < num_species; s++) {
                species[s].eos->update_eint_from_T(&species[s].prim[0], num_cells + 2);
                species[s].eos->update_p_from_eint(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_auxillary(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_conserved(&species[s].prim[0], &species[s].u[0], num_cells + 2);
            }
            //cout<<endl;
        }
    //}
//}

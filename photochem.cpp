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
    double x = 2 * 157807 / T_e;
    // Case B
    return 2.753e-14 * pow(x, 1.500) / pow(1 + pow(x / 2.740, -0.407), 2.242);
}
double H_threebody_recombination(double T_e) {
    using std::pow;

    double x = 2 * 157807 / T_e;
    return (1.005e-14 / (T_e * T_e * T_e)) * pow(x, -1.089) /
           pow(1 + pow(x / 0.354, 0.874), 1.101);
}
double H_collisional_ionization(double T_e) {
    using std::exp;
    using std::pow;

    if (T_e < 220) return 0 ;

    double x = 2 * 157807. / T_e;
    double term = 21.11 * pow(T_e, -1.5) * pow(x, -1.089) /
                  pow(1 + pow(x / 0.354, 0.874), 1.101);
    return term * exp(-0.5 * x);
}

// Cooling rate per electron.
double HOnly_cooling(const std::array<double, 3> nX, double Te) {
    using std::exp;
    using std::pow;
    using std::sqrt;

    const double T_HI = 157807.;
    double x = 2 * T_HI / Te;
    double cooling = 0;

    // Recombination Cooling:
    double term ;
    if (x < 1e5) 
        term = 3.435e-30 * Te * pow(x, 1.970) / pow(1 + pow(x / 2.250, 0.376), 3.720);
    else 
        term = 7.562-27 * std::pow(Te, 0.42872) ;
    cooling += nX[1] * term;
    

    // Collisional ionization cooling:
    term = kb * T_HI * H_collisional_ionization(Te);
    cooling += nX[0] * term;

    // HI Line cooling (Lyman alpha):
    term = 7.5e-19 * exp(-0.75 * T_HI / Te) / (1 + sqrt(Te / 1e5));
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
        C = std::max(C, 1e-20*R) ;
    };

    // Computes <x_e> - xe, the residual in the average electron abundance
    double operator()(double xe) const { return _xe_average(xe) - xe; }

   private:
    // Compute the average electron abundance, given a guess for the average
    // electron abundance.
    double _xe_average(double x) const {
        double ne = ne0 + (x - x0) * nH;

        // Compute equilibrium
        double ion = 0;
        for(int b=0; b < num_he_bands; b++)
            ion += Gamma0[b]/nH * -std::expm1(-tau0[b] * (1 - x));
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
        double ion = 0.;
        
        for(int b=0; b < num_he_bands; b++)
            Gamma0[b]/nH * -std::expm1(-tau0[b] * (1 - x_bar));
        
        double x_eq = (ion + ne * C) / (ion + ne * (C + R + ne * B));
        double t_rat = dt * (ion + ne * (C + R + ne * B));

        double x_new = x_eq + (x0 - x_eq) * std::exp(-t_rat);

        return {nH * (1 - x_new), nH * x_new, nH * (x_new - x0) + ne0};
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
                        const std::array<double, 3>& nX_average,
                        const std::array<double, 3>& TX,
                        const Matrix_t& collisions)
     : GammaH(GammaH_),
       dt(dt_),
       nX(nX_average),
       T(TX),
       coll_mat(Mat3::Identity() - collisions * dt){};

    // Compute the temperature residual.
    double operator()(double Te) const { return _compute_T(Te)(2) - Te; }

    // Compute the cooling rate give the electron temperature.
    std::array<double, 3> net_heating_rate(double Te) const {
        //return heating_rate(Te) - cooling_rate(Te);
        return {0.0, 0.0, GammaH - nX[2]*HOnly_cooling(nX, Te)};
    }
    
    std::array<double, 3> heating_rate(double Te) const {
        return {0.0, 0.0, GammaH};
    }
    
    std::array<double, 3> cooling_rate(double Te) const {
        return {0.0, 0.0, - nX[2]*HOnly_cooling(nX, Te)};
    }
    

   private:
    Vec3 _compute_T(double Te) const {
        // Setup RHS
        Vec3 RHS(T.data());
        RHS(2) += (GammaH/nX[2] - HOnly_cooling(nX, Te)) * dt / (1.5 * kb);

        // Solve for temperature
        return coll_mat.partialPivLu().solve(RHS);
    }

   public:
    double GammaH, dt;

    std::array<double, 3> nX, T;
    Mat3 coll_mat;
};

void c_Sim::do_photochemistry() {
    // cout<<" Doing photochemistry "<<endl;
            
            //double tau0 = 0;
        
            for (int b = 0; b < num_he_bands; b++) {
                radial_optical_depth(num_cells + 1,b)         = 0.;
                radial_optical_depth_twotemp(num_cells + 1,b) = 0.;
            }
            
            for (int j = num_cells + 1; j > 0; j--) {
                
                for (int b = 0; b < num_he_bands; b++) {
                    cell_optical_depth(j,b)           = 0.;
                    total_opacity(j,b)                = 0 ;
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

                // Remove negative densities/temperatures
                for (int s=0; s < 3; s++) {
                    nX[s] = std::max(nX[s],0.0);
                    TX[s] = std::max(TX[s],1.0);
                }

                // The local optical depth to ionizing photons
                // under the "initially all neutral" assumption

                // Assume 20eV photons for now
                //double E_phot = 20 * ev_to_K * kb;
                std::vector<double> dtaus = np_zeros(num_he_bands);
                std::vector<double> Gamma0  = np_zeros(num_he_bands);
                //std::array<double, num_he_bands> Gamma0;
                
                for (int b = 0; b < num_he_bands; b++) {
                    Gamma0[b] = 0.25 * solar_heating(b) / photon_energies[b] * std::exp(-radial_optical_depth_twotemp(j,b)) / dx[j];
                    dtaus[b]  = species[0].opacity_twotemp(j, b) * (nX[0] + nX[1]) * mX[0] * dx[j];
                }
                // Update the ionization fractions
                C2Ray_HOnly_ionization ion(&Gamma0[0], &dtaus[0], dt, nX, TX[2], num_he_bands);
                Brent root_solver;
                double xmin = std::max((nX[1] - nX[2]) / (nX[0] + nX[1]), 0.0);
                
                if(!(ion(xmin)*ion(1)<=0)) {
                    std::cout << j << " " << nX[0] << " " << nX[1] << " " << nX[2] << "\n"
                        << "\t" << xmin << " " << 1 << "," << ion(xmin) << " " << ion(1) << "\n" 
                        << "\t" << TX[0] << " " << TX[1] << " " << TX[2] << "\n" ;
                }

                double x_bar = root_solver.solve(xmin, 1, ion);

                std::array<double, 3> nX_new = ion.nX_new(x_bar);
                std::array<double, 3> nX_bar = ion.nX_average(x_bar);

                // Store the optical depth through this cell:
                for (int b = 0; b < num_he_bands; b++) {
                    dtaus[b] *= (1 - x_bar);
                    if(j<num_cells+1)
                        radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + dtaus[b];
                    else
                        radial_optical_depth_twotemp(j,b) = dtaus[b];
                }
                //tau0 += dtau;

                // Next update the primitive quantities (conserved are done at
                // the end)

                if (nX_new[0] > nX[0]) {
                    // Net recombination - add the energy / momentum to the H atoms
                    double f1 = species[0].cv / species[1].cv;
                    double f2 = species[0].cv / species[2].cv;
                    uX[0] += (1 - nX[0] / nX_new[0]) * (f1*uX[1] + f2*uX[2] - uX[0]);
                    TX[0] += (1 - nX[0] / nX_new[0]) * (   TX[1] +    TX[2] - TX[0]);
                    f1 = mX[1]/mX[0] ; f2 = mX[2]/mX[0] ;
                    vX[0] += (1 - nX[0] / nX_new[0]) * (f1*vX[1] + f2*vX[2] - vX[0]);
                } else {
                    // Net ionization - add the energy / momentum to the p+/e-
                    double f1 = species[1].cv / species[0].cv;
                    uX[1] += (1 - nX[1] / nX_new[1]) * (f1*uX[0] - uX[1]);
                    uX[2] += (1 - nX[2] / nX_new[2]) * ( 0*uX[0] - uX[2]);
                    TX[1] += (1 - nX[1] / nX_new[1]) * (   TX[0] - TX[1]);
                    TX[2] += (1 - nX[2] / nX_new[2]) * (      0  - TX[2]);

                    f1 = mX[0]/(mX[1]+mX[2]);
                    vX[1] += (1 - nX[1] / nX_new[1]) * (f1*vX[0] - vX[1]);
                    vX[2] += (1 - nX[2] / nX_new[2]) * (f1*vX[0] - vX[2]);
                }

                for (int s = 0; s < 3; s++) {
                    species[s].prim[j].number_density = nX_new[s];
                    species[s].prim[j].speed = vX[s] ;
                    species[s].prim[j].density = nX_new[s] * mX[s];
                    species[s].prim[j].internal_energy = uX[s];
                    species[s].prim[j].temperature = TX[s];
                }

                double GammaH = 0.;
                for (int b = 0; b < num_he_bands; b++) {
                    GammaH += 0.25 * solar_heating(b) * (1 - 13.6 / photon_energies[b]) *
                                std::exp(-radial_optical_depth_twotemp(j,b)) * (-std::expm1(-dtaus[b])) / dx[j];
                }

                // Collisional heat exchange:
                fill_alpha_basis_arrays(j);
                compute_collisional_heat_exchange_matrix(j);

                // Solve for radiative cooling implicitly
                C2Ray_HOnly_heating heat(GammaH, dt, nX_bar, TX,
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
                    std::cout << j << " " << Te1 << " " << Te2 << " " << TX[2] << "\n"  
                        << "\t" << heat(Te1) << " " << heat(Te2) << "\n" 
                        << "\t" << nX_bar[0] << " " << nX_bar[1] << " " << nX_bar[2] << "\n" 
                        << "\t" << TX[0] << " " << TX[1] << " " << TX[2] << "\n" 
                        << "\t\t" << heat.coll_mat(0,0) << " " <<  heat.coll_mat(0,1) << " " << heat.coll_mat(0,2) << "\n"
                        << "\t\t" << heat.coll_mat(1,0) << " " <<  heat.coll_mat(1,1) << " " << heat.coll_mat(1,2) << "\n"
                        << "\t\t" << heat.coll_mat(2,0) << " " <<  heat.coll_mat(2,1) << " " << heat.coll_mat(2,2) << "\n";
                }

                double Te = root_solver.solve(Te1, Te2, heat);
                std::array<double, 3> heating = heat.heating_rate(Te);
                std::array<double, 3> cooling = heat.cooling_rate(Te);

                /*
                std::cout << "\t" << TX[2] << " " << Te1 << " " << Te2 << "\n" ; 
                std::cout << "\t" << Te << "\n" ;
                */

                //TODO: This might need to be called for all bands. When H, p+ and e- densities changes, so do their opacities across all bands, including low-energy bands.
                //TODO: This needs to be ~k_rosseland, not k_Planck,solar. However k_rosseland in the UV is needed for this.
                
                for (int b = 0; b < num_he_bands; b++) {   
                    total_opacity(j,b)                = species[0].opacity(j,b)          * (nX[0] + nX[1]) * mX[0]; //TODO: Assign sensible opacities to electrons and protons
                    total_opacity_twotemp(j,b)        = species[0].opacity_twotemp(j, b) * (nX[0] + nX[1]) * mX[0] ;
                    cell_optical_depth(j,b)           = dtaus[b];
                
                    S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j,b)); //No factor 0.25 here, as we want to see the unaveraged radiation we put in, consistent with low-energy bands
                    dS_band(j,b) = 0.25 * solar_heating(b) * (1 - 13.6 / photon_energies[b]) *
                                std::exp(-radial_optical_depth_twotemp(j,b)) * (-std::expm1(-dtaus[b])) / dx[j];
                }
                    
                for (int s = 0; s < 3; s++) { 
                    //Even in cases when radiative transport is not used, we want to document those quanitites in the output
                    species[s].dG(j) = cooling[s];
                    species[s].dS(j) = heating[s];
                    //When the radiative transport is not used, we update the internal energy here in place
                    //if (!use_rad_fluxes) 
                    //    species[s].prim[j].internal_energy += (heating[s] + cooling[s]) * dt / species[s].prim[j].density ;
                }
            }

            // Finally, lets update the conserved quantities
            for (int s = 0; s < num_species; s++) {
                species[s].eos->update_p_from_eint(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_auxillary(&species[s].prim[0], num_cells + 2);
                species[s].eos->compute_conserved(&species[s].prim[0], &species[s].u[0], num_cells + 2);
            }
        }
    //}
//}

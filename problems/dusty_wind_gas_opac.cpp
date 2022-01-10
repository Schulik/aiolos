
#include <stdexcept>

#include "aiolos.h"

#include "problems/condensation.h"

class DustOpacity {
   public:
    DustOpacity(double size_cm, double BETA = 1, double RHO = 3)
     : _k0(3 / (4 * RHO * size_cm)),
       _a(size_cm),
       _T0(1 / (2 * M_PI * 3.475 * size_cm)),
       _BETA(BETA) {}

    double kappa_wle(double wle_micron) const {
        double x = 1e4 * 2 * M_PI * _a / wle_micron;

        if (_BETA != 1) x = std::pow(x, _BETA);

        return _k0 * std::min(1., x);
    }
    double kappa_wle(double wle0, double wle1) const {
        double x0 = wle0 / (1e4 * 2 * M_PI * _a) ;
        double x1 = wle1 / (1e4 * 2 * M_PI * _a) ;

        double fac ;
        if (_BETA == 1) {
            auto f = [](double x) {
                if (x <= 1)
                    return x ;
                else 
                    return 1 + std::log(x) ;
            } ;

            fac = (f(x1) - f(x0)) / (x1 - x0) ;
        } else {
            auto f = [this](double x) {
                if (x <= 1)
                    return x ;
                else 
                    return 1 + (std::pow(x, 1-_BETA)-1)/(1-_BETA)  ;
            } ;

            fac = (f(x1) - f(x0)) / (x1 - x0) ;
        }

        return _k0 * fac ;
    }

    double kappa_Planck(double T) const {
        double k;
        if (_BETA == 1)
            k = _k0 * T / std::sqrt(_T0*_T0 + T*T);
        else
            k = _k0 / std::pow(1 + _T0*_T0/(T*T), 0.5*_BETA);

        return k;
    }

    double kappa_Ross(double T) const { return kappa_Planck(T); }

   private:
    double _k0, _a, _T0, _BETA;
};

class GasOpacity {
   public:
    static const int num_wle_bins = 41;
    static const double wle[];
    static const double kappa_2000[];

    static double kappa_Planck(double T) { return kappa_eval(T, logKpl); }
    static double kappa_Ross(double T) { return kappa_eval(T, logKRoss); }

   private:
    static double kappa_eval(double T, const double* logKappa) {
        double logT = std::log(T);

        double logK;
        if (logT <= logT0)
            logK = logKappa[0];
        else if (logT >= logT1)
            logK = logKappa[num_pl - 1];
        else {
            double f = (num_pl - 1) * ((logT - logT0) / (logT1 - logT0));
            int i = static_cast<int>(f);
            f -= i;
            logK = logKappa[i] * (1 - f) + f * logKappa[i + 1];
        }

        return std::exp(logK);
    }

    static const int num_pl = 300;
    static constexpr double logT0 = 5.703782475;
    static constexpr double logT1 = 8.699514748;

    static const double logKpl[];
    static const double logKRoss[];
};

// Setup the condensation parameters
//    Perez-Becker & Chiang (2013)
Condensible silicate(169, 3.21e10, 6.72e14, 0.1);
GasOpacity kappa_gas;
DustOpacity kappa_dust(1);

static const double RHO_DUST = 3;
static const double THICKNESS = 10;

static double T_core;
static double f_dust = 1e-15;
static double a_grain = -1;

void c_Sim::user_output_function(int output_counter) {
    std::ofstream file;
    std::string filename = workingdir + "user_output.dat";

    if (output_counter == 0) {
        file.open(filename);
        file << "# snap time[s] T_surf[K] L_surf[K]\n";
    } else {
        file.open(filename, std::ios::app);
    }
    file << output_counter << " " << globalTime << " " << T_core << " "
         << 4 * pi * R_core * R_core * F_core << "\n";
    double T_F = std::pow(std::abs(F_core) / sigma_rad, 0.25);
    if (F_core < 0) T_F *= -1;
    cout << "T_core, Net thermal flux[K]: " << T_core << ", " << T_F << "K\n";
}

void c_Species::user_boundary_left(std::vector<AOS>& u) {
    for (int j = 0; j < base->num_ghosts; j++) {
        u[j] = u[2 * base->num_ghosts - (j + 1)];
        u[j].u2 *= -1;
        base->phi[j] = base->phi[base->num_ghosts];
    }
}

void c_Species::user_boundary_right(std::vector<AOS>& u) {
    int N = num_cells + 2 - base->num_ghosts;
    for (int j = 0; j < base->num_ghosts; j++) {
        u[N + j] = u[N - 1];
        base->phi[N + j] = base->phi[N - 1];
    }
}

void c_Species::user_initial_conditions() {
    // Set the input bands:
    assert(base->num_bands_in == kappa_gas.num_wle_bins);
    for (int i = 1; i < kappa_gas.num_wle_bins; i++) {
        base->l_i_in[i] = kappa_gas.wle[i];
        base->l_i12_in[i - 1] =
            std::sqrt(base->l_i_in[i - 1] * base->l_i_in[i]);
    }

    if (is_dust_like && a_grain < 0) {
        a_grain = std::pow(3 / (4 * M_PI) * mass_amu * amu / RHO_DUST, 1 / 3.);
        kappa_dust = DustOpacity(a_grain);
        f_dust = initial_fraction;
    }

    // Set the starting temperature a little below the equilibrium temperature
    T_core = base->T_star;
    double S = 0.25 * pow(T_core, 4.) * std::pow(base->R_star * rsolar, 2.) /
               std::pow(base->planet_semimajor * au, 2.);
    T_core = pow(S / 2, 0.25);

    // Set the planet temperature
    base->T_surface = T_core;
    base->L_core =
        4 * pi * base->R_core * base->R_core * sigma_rad * std::pow(T_core, 4);
    base->use_planetary_temperature = 1;

    if (is_dust_like)
        std::cout << "Dusty wind setup:\n\tT_eq=" << T_core << "K\n"
                  << "\tGrain size=" << a_grain * 1e4 << "micron\n";

    // First get the boundary condition
    AOS_prim p0;
    p0.pres = initial_fraction * (30 / mass_amu) * silicate.P_vap(T_core);
    p0.density = p0.pres * mass_amu / (Rgas * T_core);
    p0.speed = 0;

    // Assume stratification with the gas:
    // Increase Bondi-radius slightly to over-pressure the gas
    double cs2 = Rgas * T_core / 30;
    double Rb = 1.1 * G * base->planet_mass / cs2;
    std::cout << "Bondi Radius:" << Rb / 1.01 << "cm\n";
    std::cout << initial_fraction << "\n";

    for (int j = 0; j <= num_cells + 1; j++) {
        double f = std::min(1., std::exp(Rb / base->x_i12[j] -
                                         Rb / base->x_i12[base->num_ghosts]));
        AOS_prim p = p0;
        p.density *= f;
        p.pres *= f;
        eos->compute_conserved(&p, &u[j], 1);
    }
}

void c_Species::user_opacity() {
    const auto& l_i_in = base->l_i_in;

    if (is_dust_like) {
        for (int j = 0; j < num_cells + 2; j++) {
            for (int b = 0; b < num_bands_in; b++) {
                //double wle = std::sqrt(l_i_in[b] * l_i_in[b + 1]);
                opacity_twotemp(j, b) = kappa_dust.kappa_wle(l_i_in[b], l_i_in[b + 1]);
            }
            for (int b = 0; b < num_bands_out; b++) {
                opacity(j, b) = kappa_dust.kappa_Ross(prim[j].temperature);
                opacity_planck(j, b) =
                    kappa_dust.kappa_Planck(prim[j].temperature);
            }
        }
    } else {
        // Gas
        for (int j = 0; j < num_cells + 2; j++) {
            for (int b = 0; b < num_bands_in; b++) {
                opacity_twotemp(j, b) = kappa_gas.kappa_2000[b];
            }
            for (int b = 0; b < num_bands_out; b++) {
                opacity(j, b) = kappa_gas.kappa_Ross(prim[j].temperature);
                opacity_planck(j, b) =
                    kappa_gas.kappa_Planck(prim[j].temperature);
            }
        }
    }
}
void c_Sim::user_heating_function() {

    if (num_species != 2)
        throw std::runtime_error("Condensation requires num_species=2") ;

    assert(a_grain > 0) ;
    
    SingleGrainCondensation cond(
        silicate, species[0].mass_amu, a_grain, RHO_DUST, 
        species[0].cv, species[1].cv) ;
        
    for (int b = 0; b < num_bands_in; b++) {
        radial_optical_depth_twotemp(num_cells + 1,b) = 0.;
    }

    RadiationProperties rad(num_bands_in, num_bands_out) ;


    for (int j = num_cells + 1; j >= 0; j--) {

        // Get the local conditions:
        std::array<double, 2> rho = 
            { species[0].prim[j].density, species[1].prim[j].density } ;
        std::array<double, 2> T =   
            { species[0].prim[j].temperature, species[1].prim[j].temperature };
        std::array<double, 2> v =   
            { species[0].prim[j].speed, species[1].prim[j].speed } ;

        if (use_rad_fluxes == 1) {

            // Set the stellar radiation properties:
            for (int b=0; b < num_bands_in; b++) {
                rad.kappa_star[b][0] = species[0].opacity_twotemp(j, b) ;
                rad.kappa_star[b][1] = species[1].opacity_twotemp(j, b) ;

                if (j < num_cells + 1)
                    rad.flux_star[b] = 0.25 * solar_heating(b) * 
                        std::exp(-radial_optical_depth_twotemp(j+1,b)) ;
                else
                    rad.flux_star[b] = 0.25 * solar_heating(b) ;
            }

            // Set the thermal radiation properties
            for (int b=0; b < num_bands_out; b++) {
                rad.J_thermal[b] = Jrad_FLD(j, b) ;

                rad.kappa_thermal[b][0] = species[0].opacity_planck(j, b) ;
                rad.kappa_thermal[b][1] = species[1].opacity_planck(j, b) ;

                rad.f_band[b][0] = compute_planck_function_integral3(l_i_out[b], l_i_out[b+1] ,T[0]) ;
                rad.f_band[b][1] = compute_planck_function_integral3(l_i_out[b], l_i_out[b+1] ,T[1]) ;
            }

            // Cell size (for the optical depth)
            rad.dr = dx[j]; 

            // Solve for the new state
            cond.set_state(rho, T, v, rad, dt) ;
        }
        else {
            // No radiation
            cond.set_state(rho, T, v, dt) ;
        }

        double d0, d1 ; 
        std::tie(d0,d1) = cond.bracket_solution() ;
        
        Brent brent(1e-10*rho[1]) ;
        d0 = brent.solve(d0, d1, cond) ;

        std::tie(rho, T, v) = cond.update_T_rho_v(d0) ;

        // Store the new values
        //     Set the new density and velocity
        species[0].prim[j].density = rho[0] ;
        species[1].prim[j].density = rho[1] ;

        species[0].prim[j].speed = v[0] ;
        species[1].prim[j].speed = v[1] ;

        if (T[0] < 0 || T[1] < 0) {
            std::cout << j << " " << dt 
                      << " (" << rho[0] << " " << rho[1] << ")"
                      << " (" << T[0] << " " << T[1] << ")"
                      << " (" << v[0] << " " << v[1] << ")\n" ;
            std::cout << "\t" << species[0].prim[j].temperature << " " << species[1].prim[j].temperature << "\n" ;

            std::array<double,2> heating, cooling ;
            heating = rad.compute_stellar_heating(rho) ;
            cooling = cond.net_heating_rate(rho, T, v) ;

            std::cout << "\t(" << heating[0] << " " << heating[1] << ")"
                      <<  " (" << cooling[0] << " " << cooling[1] << ")\n" ;

            std::cout << "\t(" << (heating[0] + cooling[0])*dt/ (rho[0]*species[0].cv)
                      <<  " " << (heating[1] + cooling[1])*dt / (rho[1]*species[1].cv) << ")\n" ;
        }


        //     Compute the net heating/cooling rate:
        if (use_rad_fluxes) {
            std::array<double,2> heating, cooling ;

            heating = {0, 0 }; // rad.compute_stellar_heating(rho) ;
            cooling = cond.net_heating_rate(rho, T, v) ;

            for (int s=0; s < 2; s++) {
                species[s].dS(j) = heating[s] ;
                species[s].dG(j) = cooling[s] ;
            }

            //     Save the new optical depth etc
            for (int b=0; b < num_bands_in; b++) {

                total_opacity_twotemp(j, b) = rho[0]*rad.kappa_star[b][0] + rho[1]*rad.kappa_star[b][1] ;
                cell_optical_depth_twotemp(j, b) = dx[j]*total_opacity_twotemp(j,b) ;

                if (j < num_cells+1)
                    radial_optical_depth_twotemp(j,b) = cell_optical_depth_twotemp(j, b) +
                        radial_optical_depth_twotemp(j+1, b) ;
                else
                    radial_optical_depth_twotemp(j,b) = 0 ;

                S_band(j,b) = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j,b)) ;
                if (j == num_cells+1)
                    dS_band(j,b) = 0.25 * solar_heating(b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
                else
                    dS_band(j,b) = 0.25 * S_band(j+1,b) * (-std::expm1(-cell_optical_depth_twotemp(j, b)) / dx[j]) ;
            }
            // Same for the thermal bands
            for (int b=0; b < num_bands_out; b++) {
                total_opacity(j, b) = rho[0]*rad.kappa_thermal[b][0] + rho[1]*rad.kappa_thermal[b][1] ;
                cell_optical_depth(j, b) = dx[j]*total_opacity(j,b) ;

                if (j < num_cells+1)
                    radial_optical_depth(j,b) = cell_optical_depth(j, b) + radial_optical_depth(j+1, b) ;
                else
                    radial_optical_depth(j,b) = 0 ;
            }
        }
        else {
            // Without radiation, just set the new temperatues
            for (int s=0; s < 2; s++) 
                species[s].prim[j].temperature = T[s] ;
        }
    }



    // Finally, lets update the conserved quantities
    for (int s = 0; s < num_species; s++) {
        species[s].eos->update_eint_from_T(&species[s].prim[0], num_cells + 2);
        species[s].eos->update_p_from_eint(&species[s].prim[0], num_cells + 2);
        species[s].eos->compute_auxillary(&species[s].prim[0], num_cells + 2);
        species[s].eos->compute_conserved(&species[s].prim[0], &species[s].u[0], num_cells + 2);
    }

    // Recompute radiative heating etc given new densities.
    update_dS() ;
}

void c_Sim::user_loop_function() {

    // Set energy control (efficiency + initial stability)

    if (num_species != 2)
        throw std::runtime_error("Condensation requires num_species=2") ;

    // Compute evaporation and condensation from the surface

    RadiationProperties rad(num_bands_in, 1) ;

    int j = num_ghosts ;
    SurfaceCondensation planet_surf(silicate, species[0].mass_amu, 
                                   species[0].cv, species[1].cv, RHO_DUST*THICKNESS,surf[j-1]/vol[j]);
    if (use_rad_fluxes == 1) {
        // Stellar radiation reaching the planet
        for (int b=0; b < num_bands_in; b++) 
            rad.flux_star[b] = 0.25 * solar_heating(b) * 
                    std::exp(-radial_optical_depth_twotemp(j,b)) ;

        // Thermal radiation reaching the planet (F_down = \pi J in diffusion limit)
        rad.J_thermal[0] = - F_core / pi ;
    } else {
        throw std::runtime_error("Radiation must be turned on") ;
    }

    // Energy due to rain out
    double L_rain = 0 ;
    if (species[1].prim[j].speed < 0) {
        AOS_prim& prim = species[1].prim[j] ;

        double E_tot = prim.density*
            (species[1].cv*prim.temperature + 0.5*prim.speed*prim.speed) ;

        L_rain = - E_tot * prim.speed ;

        prim.density += prim.density * prim.speed * dt * surf[j-1]/vol[j] ;
    } 

    // Get the surface temperature and mass-loss rate
    planet_surf.set_state(T_core,
                          species[0].prim[j].temperature, 
                          species[0].prim[j].pres,
                          rad, L_rain, dt) ;

    double T0, T1 ; 
    std::tie(T0,T1) = planet_surf.bracket_solution() ;
    
    Brent brent(1e-10*T_core) ;

    T_core = brent.solve(T0, T1, planet_surf) ;
    
    L_core = 4*pi*R_core*R_core*sigma_rad*std::pow(T_core, 4) ;

    /*
    double Fstar = 0 ;
    for (int b=0; b < num_bands_in; b++)
        Fstar += rad.flux_star[b] ;

    std::cout 
        << "t, dt, T " << globalTime << " " << dt << " " << T_core << ", " 
        << Fstar << " " << pi*rad.J_thermal[0] << " "
        << L_core/(4*pi*R_core*R_core) << " " << L_rain << "\n" ;
    */

    double drho = dt*planet_surf.mass_flux(T_core)*surf[j-1]/vol[j] ;
    double heat = dt*planet_surf.gas_heating_rate(T_core)*surf[j-1]/vol[j] ;

    std::array<double, 2> frac = {1-f_dust, f_dust} ;

    for (int s=0; s < 2; s++) {
        // Update the mass, momentum and energy of the corresponding cell
        if (frac[s] == 0) continue ;

        AOS_prim& prim = species[s].prim[j] ;

        double E_tot = prim.density*
            (species[s].cv*prim.temperature + 0.5*prim.speed*prim.speed) ;

        E_tot += frac[s]*heat;

        prim.speed /= (1 + frac[s]*drho/prim.density) ;
        prim.density += frac[s]*drho ;
        
        E_tot -= 0.5*prim.density*prim.speed*prim.speed ;

        prim.temperature = E_tot / (prim.density*species[s].cv) ;
    }

    for (int s = 0; s < num_species; s++) {
        species[s].eos->update_eint_from_T(&species[s].prim[j], 1);
        species[s].eos->update_p_from_eint(&species[s].prim[j], 1);
        species[s].eos->compute_auxillary(&species[s].prim[j], 1);
        species[s].eos->compute_conserved(&species[s].prim[j], &species[s].u[j], 1);
    }
} ;

// Opacity tables

const double GasOpacity::wle[] = {
    // num_wle_bins + 1 values (edges)
    0.010000000,  0.242257563,   0.256266246,  0.283554614,  0.291965833,
    0.305052789,  0.318010461,   0.327075819,  0.346768039,  0.365996188,
    0.370967563,  0.382831630,   0.392860424,  0.395075126,  0.419332568,
    0.434205492,  0.474003174,   0.511663427,  0.529811150,  0.537611747,
    0.556054279,  0.598885433,   0.635656740,  0.682314145,  0.721769353,
    0.786154333,  0.915009875,   1.034304556,  1.167838565,  1.262043193,
    1.589175551,  3.036839747,   5.157061894,  7.498942093,  9.295263782,
    10.192970827, 12.119817876,  14.508452552, 20.535250265, 49.192005801,
    61.044217525, 100.000000000,
};
const double GasOpacity::kappa_2000[] = {
    // num_wle_bins values (centers)
     79025.400268919,  753885.289333429,    2524.086512977, 1241175.364814282,
      1377.367925671,      19.015425556,       1.330771286,       0.225716482,
         1.910196628,      54.162601755,    3177.357963578,   64067.646275615,
    531767.077900031,       0.093689448,       1.205861884,      14.624020316,
        76.671136863,     501.438674192,    3919.420932684,      91.578418347,
        27.041404213,      11.825018368,       5.407402349,       2.081863323,
         0.667285222,       0.202218618,       0.633024355,       2.044300897,
         4.608498637,      11.310124638,      30.370486352,       8.874371758,
         1.231658192,     103.849113784,      10.293226342,       0.864521700,
         0.386907335,       1.584981016,       0.154267384,       0.811545320,
        88.453348420,
};

const double GasOpacity::logKpl[] = {
    2.578267188, 2.595515552, 2.612292810, 2.628601556, 2.644445530,
    2.659828708, 2.674754709, 2.689228384, 2.703254462, 2.716838424,
    2.729985621, 2.742701890, 2.754993695, 2.766867279, 2.778329888,
    2.789388502, 2.800050741, 2.810324750, 2.820218705, 2.829740708,
    2.838900346, 2.847706617, 2.856168974, 2.864297734, 2.872102820,
    2.879594994, 2.886785333, 2.893685144, 2.900306275, 2.906660396,
    2.912760028, 2.918618003, 2.924247228, 2.929661083, 2.934873177,
    2.939897347, 2.944747899, 2.949439176, 2.953986287, 2.958403528,
    2.962706332, 2.966910095, 2.971030008, 2.975081553, 2.979080084,
    2.983041541, 2.986981117, 2.990914391, 2.994856594, 2.998822844,
    3.002827907, 3.006886719, 3.011013174, 3.015221084, 3.019524134,
    3.023934901, 3.028465804, 3.033128349, 3.037933315, 3.042891139,
    3.048010910, 3.053301118, 3.058769313, 3.064422087, 3.070265105,
    3.076302733, 3.082538732, 3.088975173, 3.095613756, 3.102454405,
    3.109496424, 3.116737609, 3.124175023, 3.131804317, 3.139620438,
    3.147616963, 3.155786620, 3.164121228, 3.172611628, 3.181247830,
    3.190018872, 3.198913376, 3.207918892, 3.217022562, 3.226210939,
    3.235469930, 3.244785242, 3.254141872, 3.263524799, 3.272918540,
    3.282307591, 3.291676141, 3.301008483, 3.310288760, 3.319501240,
    3.328630242, 3.337660119, 3.346575553, 3.355361344, 3.364002665,
    3.372484853, 3.380793637, 3.388915120, 3.396835650, 3.404542165,
    3.412021944, 3.419262659, 3.426252493, 3.432980119, 3.439434619,
    3.445605612, 3.451483179, 3.457057949, 3.462320994, 3.467263925,
    3.471878852, 3.476158339, 3.480095505, 3.483683944, 3.486917664,
    3.489791256, 3.492299746, 3.494438590, 3.496203694, 3.497591502,
    3.498598781, 3.499222775, 3.499461135, 3.499311958, 3.498773624,
    3.497845000, 3.496525271, 3.494813980, 3.492711015, 3.490216595,
    3.487331298, 3.484055913, 3.480391653, 3.476339906, 3.471902389,
    3.467081058, 3.461878164, 3.456296126, 3.450337637, 3.444005625,
    3.437303225, 3.430233749, 3.422800705, 3.415007869, 3.406859111,
    3.398358548, 3.389510450, 3.380319281, 3.370789684, 3.360926482,
    3.350734725, 3.340219622, 3.329386625, 3.318241415, 3.306789914,
    3.295038325, 3.282993171, 3.270661281, 3.258049907, 3.245166720,
    3.232019925, 3.218618271, 3.204971196, 3.191088889, 3.176982412,
    3.162663862, 3.148146486, 3.133444880, 3.118575311, 3.103555888,
    3.088406679, 3.073150066, 3.057811326, 3.042419242, 3.027006125,
    3.011608166, 2.996266330, 2.981027407, 2.965944442, 2.951077349,
    2.936494126, 2.922271649, 2.908496024, 2.895264202, 2.882686542,
    2.870887297, 2.860002976, 2.850183376, 2.841595977, 2.834428097,
    2.828882351, 2.825173626, 2.823532162, 2.824203878, 2.827444602,
    2.833517635, 2.842690810, 2.855228568, 2.871384112, 2.891394396,
    2.915477690, 2.943824146, 2.976580751, 3.013846807, 3.055679030,
    3.102084378, 3.153019278, 3.208383328, 3.268023933, 3.331765898,
    3.399415177, 3.470721739, 3.545401000, 3.623172975, 3.703724434,
    3.786735337, 3.871942308, 3.959103604, 4.047942961, 4.138199317,
    4.229639433, 4.322034787, 4.415177884, 4.508891623, 4.603016332,
    4.697394827, 4.791873562, 4.886328830, 4.980681781, 5.074847782,
    5.168717873, 5.262202057, 5.355235319, 5.447748751, 5.539689487,
    5.631030426, 5.721732937, 5.811752732, 5.901056017, 5.989615138,
    6.077402881, 6.164398320, 6.250587676, 6.335958296, 6.420468005,
    6.504081288, 6.586826349, 6.668740850, 6.749789678, 6.829907350,
    6.909089623, 6.987400904, 7.064861973, 7.141421908, 7.217047529,
    7.291754584, 7.365554051, 7.438443452, 7.510415242, 7.581469544,
    7.651608275, 7.720827344, 7.789140930, 7.856572793, 7.923117115,
    7.988755613, 8.053494010, 8.117321504, 8.180223170, 8.242239522,
    8.303413206, 8.363729324, 8.423173893, 8.481755377, 8.539447305,
    8.596240615, 8.652167882, 8.707256318, 8.761491763, 8.814852719,
    8.867349062, 8.919007042, 8.969832150, 9.019827162, 9.068999086,
    9.117347516, 9.164872628, 9.211578414, 9.257457552, 9.302513452,
    9.346760947, 9.390200919, 9.432852571, 9.474752399, 9.515872575,

};

const double GasOpacity::logKRoss[] = {
    -1.062616020, -1.045543184, -1.028519475, -1.011545391, -0.994621317,
    -0.977747521, -0.960924158, -0.944151268, -0.927428782, -0.910756515,
    -0.894134177, -0.877561368, -0.861037582, -0.844562210, -0.828134542,
    -0.811753771, -0.795418994, -0.779129214, -0.762883348, -0.746680227,
    -0.730518600, -0.714397136, -0.698314433, -0.682269017, -0.666259348,
    -0.650283825, -0.634340786, -0.618428520, -0.602545263, -0.586689206,
    -0.570858502, -0.555051263, -0.539265571, -0.523499481, -0.507751021,
    -0.492018199, -0.476299009, -0.460591432, -0.444893438, -0.429202997,
    -0.413518073, -0.397836636, -0.382156660, -0.366476128, -0.350793038,
    -0.335105401, -0.319411246, -0.303708624, -0.287995610, -0.272270306,
    -0.256530842, -0.240775378, -0.225002110, -0.209209267, -0.193395118,
    -0.177557970, -0.161696172, -0.145808114, -0.129892232, -0.113947009,
    -0.097970974, -0.081962706, -0.065920833, -0.049844036, -0.033731050,
    -0.017580661, -0.001391713,  0.014836893,  0.031106198,  0.047417184,
     0.063770773,  0.080167822,  0.096609124,  0.113095409,  0.129627336,
     0.146205496,  0.162830406,  0.179502510,  0.196222173,  0.212989682,
     0.229805236,  0.246668950,  0.263580844,  0.280540842,  0.297548767,
     0.314604332,  0.331707137,  0.348856659,  0.366052246,  0.383293108,
     0.400578303,  0.417906736,  0.435277134,  0.452688045,  0.470137816,
     0.487624579,  0.505146238,  0.522700443,  0.540284575,  0.557895722,
     0.575530657,  0.593185808,  0.610857233,  0.628540590,  0.646231103,
     0.663923528,  0.681612117,  0.699290580,  0.716952038,  0.734588986,
     0.752193243,  0.769755905,  0.787267294,  0.804716908,  0.822093366,
     0.839384356,  0.856576575,  0.873655679,  0.890606227,  0.907411625,
     0.924054077,  0.940514533,  0.956772649,  0.972806745,  0.988593772,
     1.004109290,  1.019327451,  1.034220995,  1.048761266,  1.062918229,
     1.076660520,  1.089955503,  1.102769356,  1.115067173,  1.126813095,
     1.137970464,  1.148502005,  1.158370032,  1.167536683,  1.175964179,
     1.183615110,  1.190452736,  1.196441320,  1.201546461,  1.205735447,
     1.208977614,  1.211244694,  1.212511169,  1.212754600,  1.211955939,
     1.210099810,  1.207174759,  1.203173458,  1.198092867,  1.191934338,
     1.184703677,  1.176411138,  1.167071371,  1.156703312,  1.145330020,
     1.132978467,  1.119679291,  1.105466500,  1.090377156,  1.074451032,
     1.057730250,  1.040258915,  1.022082739,  1.003248673,  0.983804551,
     0.963798741,  0.943279819,  0.922296262,  0.900896168,  0.879126999,
     0.857035353,  0.834666758,  0.812065501,  0.789274477,  0.766335062,
     0.743287017,  0.720168407,  0.697015542,  0.673862937,  0.650743288,
     0.627687464,  0.604724510,  0.581881659,  0.559184361,  0.536656308,
     0.514319476,  0.492194169,  0.470299062,  0.448651255,  0.427266323,
     0.406158372,  0.385340091,  0.364822809,  0.344616546,  0.324730074,
     0.305170961,  0.285945631,  0.267059410,  0.248516577,  0.230320414,
     0.212473248,  0.194976500,  0.177830729,  0.161035672,  0.144590284,
     0.128492784,  0.112740686,  0.097330845,  0.082259485,  0.067522242,
     0.053114199,  0.039029915,  0.025263468,  0.011808481, -0.001341836,
    -0.014194662, -0.026757520, -0.039038252, -0.051044984, -0.062786086,
    -0.074270146, -0.085505932, -0.096502356, -0.107268443, -0.117813296,
    -0.128146060, -0.138275886, -0.148211901, -0.157963168, -0.167538659,
    -0.176947214, -0.186197510, -0.195298032, -0.204257036, -0.213082522,
    -0.221782202, -0.230363471, -0.238833384, -0.247198624, -0.255465484,
    -0.263639843, -0.271727147, -0.279732388, -0.287660094, -0.295514315,
    -0.303298611, -0.311016044, -0.318669177, -0.326260070, -0.333790280,
    -0.341260866, -0.348672393, -0.356024948, -0.363318142, -0.370551134,
    -0.377722644, -0.384830972, -0.391874023, -0.398849332, -0.405754087,
    -0.412585162, -0.419339149, -0.426012386, -0.432600997, -0.439100927,
    -0.445507977, -0.451817851, -0.458026189, -0.464128615, -0.470120781,
    -0.475998407, -0.481757333, -0.487393567, -0.492903328, -0.498283102,
    -0.503529693, -0.508640272, -0.513612438, -0.518444268, -0.523134375,
    -0.527681971, -0.532086923, -0.536349818, -0.540472024, -0.544455758,
    -0.548304151, -0.552021318, -0.555612424, -0.559083760, -0.562442811,
    -0.565698331, -0.568860414, -0.571940564, -0.574951772, -0.577908581,
};

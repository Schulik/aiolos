#define EIGEN_RUNTIME_NO_MALLOC

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <tuple>

#include "aiolos.h"
#include "brent.h"


/* class Condensible
 *
 * Represents the thermodynamic properties of a condensible species.
 * 
 * Parameters:
 *   mass   : mass of a single monomer of condensate (in amu)
 *   Lsub   : Specific latent heat of sublimation (erg / g)
 *   P0     : Vapor pressure coefficient 
 *   P_cond : Condensation probability, fraction of collisions that 'stick'
 * 
 *  Note: 
 *    P_vap(T) = P0 * exp( - Lsub*mass/(Rgas*T))
 */
class Condensible {
  public:
    Condensible(double mass_, double Lsub_, double P0_, double P_cond_)
     : mass(mass_), Lsub(Lsub_), P_stick(P_cond_), _P0(P0_)
    { } ;

    double P_vap(double T) const {
        return _P0 * std::exp(- Lsub*mass / (Rgas*T)) ;
    }

    double mass, Lsub, P_stick;
  private:
    double _P0 ;
} ;

/* struct Radiation properties
 *
 * Simple class for holding the local properties of the radiation field
 * and the respective opacities of the two species.
 * 
 *  Notes:
 *     flux_star is the stellar radiation field 
 *     J_thermal is the local thermal radiation field
 *     f_band is the fraction of thermal energy radiated into the band.
 */
struct RadiationProperties {

    RadiationProperties() {} ;
    RadiationProperties(int num_bands_star, int num_bands_thermal)
      : kappa_star(num_bands_star), flux_star(num_bands_star),
        kappa_thermal(num_bands_thermal), f_band(num_bands_thermal),
        J_thermal(num_bands_thermal)
    { } ;

    std::array<double, 2> compute_stellar_heating(std::array<double, 2> rho) const {

        std::array<double,2> S = {0, 0} ;

        // Stellar heating
        int num_stellar = flux_star.size() ;
        for (int b=0; b < num_stellar; b++) {

            double fac = 1 ;
            double tau = dr * 
                (rho[0]*kappa_star[b][0] + rho[1]*kappa_star[b][1]) ;
            if (tau > 0) {
                fac = -std::expm1(-tau) / tau ;
            }
            
            S[0] += flux_star[b] * fac * rho[0]*kappa_star[b][0] ;
            S[1] += flux_star[b] * fac * rho[1]*kappa_star[b][1] ;
        }

        return S ;
    }

    std::array<double, 2> compute_specific_thermal_heating() const {
        std::array<double,2> S = {0, 0} ;
        
        int num_thermal = J_thermal.size() ;
        for (int b=0; b < num_thermal; b++) {
            S[0] += 4*pi * J_thermal[b] * kappa_thermal[b][0] ;
            S[1] += 4*pi * J_thermal[b] * kappa_thermal[b][1] ;
        }

        return S ;
    }

    std::array<double, 2> compute_kappa_thermal() const {
        
        std::array<double,2> kappa = {0, 0} ;
        
        int num_thermal = J_thermal.size() ;
        for (int b=0; b < num_thermal; b++) {
            kappa[0] += f_band[b][0] * kappa_thermal[b][0] ;
            kappa[1] += f_band[b][1] * kappa_thermal[b][1] ;
        }

        return kappa ;
    }

    // Irradiation of a 'flat' surface
    double surface_irradiation() const {

        double S = 0 ;

        int num_stellar = flux_star.size() ;
        for (int b=0; b < num_stellar; b++)
            S += flux_star[b] ;
    
        int num_thermal = J_thermal.size() ;
        for (int b=0; b < num_thermal; b++)
            S += pi*J_thermal[b] ;

        return S ;
    }


    // Properties of stellar irradiation
    std::vector<std::array<double, 2>> kappa_star ;
    std::vector<double> flux_star ;

    // Properties of thermal radiation
    std::vector<std::array<double,2>> kappa_thermal, f_band ;
    std::vector<double> J_thermal;

    // Size of cell
    double dr = 0 ;

} ;


/* class SingleGrainCondensation
 *
 * A model for dust condensation where the dust grains are assumed to have
 * a single pre-specified size. 
 */
class SingleGrainCondensation {
  public:
    typedef std::array<double,2> arr_t ;

    SingleGrainCondensation(Condensible cond, double mu_gas,
                            double grain_size, double grain_density,
                            double C_gas, double C_l)
     : _cond(cond), _mu_gas(mu_gas), _C_v(C_gas), _C_l(C_l),
       _grain_size(grain_size), _grain_density(grain_density),
       _area(3/(grain_density*grain_size)),
       _use_radiation(false) 
    { } ;

    void set_state(arr_t rho, arr_t T, arr_t v, RadiationProperties rad,
                   double dt) {

        _rho0 = rho ;
        _T0 = T ;
        _v0 = v ;
        _dt = dt ;

        _rad = rad ;
        _use_radiation = true ;

        _u0[0] = _rho0[0]*_C_v*_T0[0] ;
        _u0[1] = _rho0[1]*_C_l*_T0[1] ;

        // Cache the radiative terms
        _Stherm = _rad.compute_specific_thermal_heating() ;
        arr_t sKT3 = _rad.compute_kappa_thermal() ;

        sKT3[0] *= sigma_rad*_T0[0]*_T0[0]*_T0[0] ;
        sKT3[1] *= sigma_rad*_T0[1]*_T0[1]*_T0[1] ;

        // Thermal heating + T-linearization term.
        _Stherm[0] += 12*sKT3[0]*_T0[0] ;
        _Stherm[1] += 12*sKT3[1]*_T0[1] ;

        // The effective heat capacity is increased by radiative cooling
        _Ctot = { _C_v + 16*sKT3[0]*dt, _C_l + 16*sKT3[1]*dt } ;
    }

    void set_state(arr_t rho, arr_t T, arr_t v, double dt) {

        _rho0 = rho ;
        _T0 = T ;
        _v0 = v ;
        _dt = dt ;

        _use_radiation = false ;

        _u0[0] = _rho0[0]*_C_v*_T0[0] ;
        _u0[1] = _rho0[1]*_C_l*_T0[1] ;

        _Ctot = {_C_v, _C_l} ;
    }

    arr_t net_heating_rate(arr_t rho, arr_t T, arr_t v) {

        // Compute the condensation rate
        double Pgas = rho[0] * ((Rgas/_mu_gas) * T[0]) ;
        double Pv = _cond.P_vap(T[1]) ;

        double cond = _cond.P_stick * _area * rho[1] * Pgas / std::sqrt(2*pi*(Rgas/_mu_gas)*T[0]) ;
        double evap = _cond.P_stick * _area * rho[1] * Pv   / std::sqrt(2*pi*(Rgas/_mu_gas)*T[1]) ;
        
        // Compute change in kinetic energy
        double dEk_dt = 
            0.5*(rho[0]*v[0]*v[0] - _rho0[0]*_v0[0]*_v0[0] +
                 rho[1]*v[1]*v[1] - _rho0[1]*_v0[1]*_v0[1]) / (_dt+1e-300) ;

        // Term due to change in density ( - drho/dt * C * T0):
        arr_t rho_term = {
            - (evap - cond) * _C_v * _T0[0],
            - (cond - evap) * _C_l * _T0[1] 
        } ;

        return {
            rho_term[0] + _C_v*(evap*T[1] - cond*T[0]) - dEk_dt, 
            rho_term[1] + _C_v*(cond*T[0] - evap*T[1]) + _cond.Lsub * (cond - evap) 
        } ;
    }
    
    double operator()(double rho_l) const {
        arr_t rho, T, v ;
        std::tie(rho, T, v) = update_T_rho_v(rho_l) ;

        double Pgas = rho[0] * ((Rgas/_mu_gas) * T[0]) ;
        double Pv = _cond.P_vap(T[1]) ;

        double cond = _cond.P_stick * _area * rho[1] * Pgas / std::sqrt(2*pi*(Rgas/_mu_gas)*T[0]) ;
        double evap = _cond.P_stick * _area * rho[1] * Pv   / std::sqrt(2*pi*(Rgas/_mu_gas)*T[1]) ;
        
        return  (rho_l - _rho0[1]) - (cond - evap)*_dt ;
    }

       
    std::tuple<arr_t, arr_t, arr_t> update_T_rho_v(double rho_l) const {


        // Update densities (and compute density change) using guess provided
        double drho = rho_l - _rho0[1] ;
        arr_t rho = { _rho0[0] + _rho0[1] - rho_l, rho_l } ;
        

        // Compute condensation rate / collision rate:
        double Pgas = rho[0] * ((Rgas/_mu_gas) * _T0[0]) ;
        double cs = std::sqrt((Rgas/_mu_gas) * _T0[0]) ;

        double cond_rate = _cond.P_stick * _area * rho[1] * 
            Pgas / std::sqrt(2*pi * (Rgas/_mu_gas) * _T0[0]) ;

        double coll_rate = std::sqrt(8/(9*pi)) * _area * rho[1] * rho[0] * (Rgas/_mu_gas) * cs ;

        // Update momentum
        double C = cond_rate*_dt ;
        double E = cond_rate*_dt - drho ;

        // Fix E < 0 (since {C, E} >= 0 by definition)
        if (E < 0) {
           C -= E ; E = 0 ;
        }
 
        double fac = 1 / (E + rho[1]) ;

        arr_t v ;
        v[0] = (_rho0[0]*_v0[0] + _rho0[1]*_v0[1]*E*fac) / (rho[0] + rho[1]*C*fac) ;
        v[1] = (_rho0[1]*_v0[1] +  v[0]*C) * fac  ;

        // Compute change in kinetic energy
        double dEk = 0.5*(rho[0]*v[0]*v[0] - _rho0[0]*_v0[0]*_v0[0] +
                          rho[1]*v[1]*v[1] - _rho0[1]*_v0[1]*_v0[1]) ;


        // Update temperature based on energy conservation
        //    Include radiative heating / cooling using the Commercon+ (2011)
        //    linearization
        arr_t u = { _u0[0] - dEk, _u0[1] + _cond.Lsub*drho } ;

        //    Add the radiation terms
        if (_use_radiation) {
            arr_t heat = _rad.compute_stellar_heating(rho) ;

            // Add the thermal heating
            heat[0] += rho[0]*_Stherm[0] ;
            heat[1] += rho[1]*_Stherm[1] ;

            u[0] += heat[0] * _dt ;
            u[1] += heat[1] * _dt ;
        }

        
        //    Solve for T
        C = (coll_rate + cond_rate*_C_v)*_dt ;
        E = (coll_rate + cond_rate*_C_v)*_dt - drho*_C_v ;

        // Fix E < 0 (since {C, E} >= 0 by definition)
        if (E < 0) {
           C -= E ; E = 0 ;
        }

        fac = 1 / (E + rho[1]*_Ctot[1]) ;

        arr_t T ;
        T[0] = (u[0] + u[1]*E*fac) / (rho[0]*_Ctot[0] + rho[1]*_Ctot[1]*C*fac);
        T[1] = (u[1] + T[0]*C)*fac ;
 
        return std::make_tuple(rho, T, v) ;
    }



    std::tuple<double, double> bracket_solution() const {
        double rho_min, rho_max ;

        double fac = 1 + 0.25 * _C_l*_T0[1]/_cond.Lsub ;

        rho_min = _rho0[1] ;
        if (this->operator()(rho_min) < 0) {
            rho_max = std::min(rho_min*fac, _rho0[0] + _rho0[1]) ;
            while(this->operator()(rho_max) < 0) {
                rho_min = rho_max ;
                rho_max *= fac ;
            }
        }
        else {
            rho_max = rho_min ;
            rho_min /= fac ;
            while(this->operator()(rho_min) > 0) {
                rho_max = rho_min ;
                rho_min /= fac ;
            }
        }
        return std::make_tuple(rho_min, rho_max) ;
    }

  private:
    Condensible _cond ;
    double _mu_gas, _C_v, _C_l, _grain_size, _grain_density, _area ;

    arr_t _rho0, _T0, _v0, _u0, _Stherm, _Ctot ;
    double _dt ;

    RadiationProperties _rad ;
    bool _use_radiation ;
};



/* class SurfaceCondensation
 *
 * A model for dust condensation / evaporation from a planet's surface.
 * 
 * Notes:
 *   The layer_mass = density * thickness, controls the effective thermal
 *   inertia of the planet.
 */
class SurfaceCondensation {
  public:
    typedef std::array<double,2> arr_t ;

    SurfaceCondensation(Condensible cond, double mu_gas,
                        double C_gas, double C_l, double layer_mass) 
     : _cond(cond), _mu_gas(mu_gas), _C_v(C_gas), _C_l(C_l), _m_layer(layer_mass)
    { } ;

    void set_state(double T_surf, double T_gas, double P_gas, 
                   RadiationProperties rad, double dt) {

        _Pgas = P_gas ;
        _Tgas = T_gas ;
        _T0 = T_surf ;
        _dt = dt ;

        _rad = rad ;

        // Compute the internal energy and heat capacity modified by radiation
        double sT3 = sigma_rad*_T0*_T0*_T0 ;
        
        _u0 = _m_layer*_C_l*_T0 + dt*(_rad.surface_irradiation() + 3*sT3*_T0) ;

        _Ctot = _m_layer*_C_l + _dt*4*sT3 ;
    }

    double mass_flux(double T_surf) const {

        double Pv = _cond.P_vap(T_surf) ;

        double cond = _cond.P_stick * _Pgas / std::sqrt(2*pi*Rgas/_mu_gas*_Tgas) ;
        double evap = _cond.P_stick * Pv    / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        return evap - cond ;
    }

    double surface_cooling_rate(double T_surf) const {
        double Pv = _cond.P_vap(T_surf) ;

        double cond = _cond.P_stick * _Pgas / std::sqrt(2*pi*Rgas/_mu_gas*_Tgas) ;
        double evap = _cond.P_stick * Pv    / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double dTdt = (T_surf-_T0)/(_dt + 1e-300) ;

        return (evap-cond)*(_cond.Lsub + (_C_v - _C_l)*T_surf) + _m_layer*_C_l*dTdt ;
    }
    double gas_heating_rate(double T_surf) const {
        double Pv = _cond.P_vap(T_surf) ;

        double cond = _cond.P_stick * _Pgas / std::sqrt(2*pi*Rgas/_mu_gas*_Tgas) ;
        double evap = _cond.P_stick * Pv    / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        return (evap-cond)*_C_v*T_surf ;
    }
    
    double operator()(double T_surf) const {

        double Pv = _cond.P_vap(T_surf) ;

        double cond = _cond.P_stick * _Pgas / std::sqrt(2*pi*Rgas/_mu_gas*_Tgas) ;
        double evap = _cond.P_stick * Pv    / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double Ctot = _Ctot + (_C_v - _C_l)*(evap-cond)*_dt ;

        return T_surf - (_u0 - _cond.Lsub*(evap-cond)*_dt) / Ctot ;
    }


    std::tuple<double, double> bracket_solution() const {
        double T_min, T_max ;

        double fac = 2; 

        T_min = _T0 ;
        if (this->operator()(T_min) < 0) {
            T_max = T_min*fac ;
            while(this->operator()(T_max) < 0) {
                T_min = T_max ;
                T_max *= fac ;
            }
        }
        else {
            T_max = T_min ;
            T_min /= fac ;
            while(this->operator()(T_min) > 0) {
                T_max = T_min ;
                T_min /= fac ;
            }
        }

        return std::make_tuple(T_min, T_max) ;
    }

  private:
    Condensible _cond ;
    double _mu_gas, _C_v, _C_l, _m_layer ;

    double _T0, _Tgas, _Pgas, _u0, _Ctot, _dt ;

    RadiationProperties _rad ;
};



/* class SurfaceCondensation2
 *
 * A model for dust condensation / evaporation from a planet's surface.
 * 
 * Notes:
 * - The layer_mass = density * thickness, controls the effective thermal
 *   inertia of the planet.
 * - Area_Vol is the (inner) surface area to volume ratio of the cell first
 *   active cell
 */
class SurfaceCondensation2 {
  public:
    typedef std::array<double,2> arr_t ;

    SurfaceCondensation2(Condensible cond, double mu_gas,
                        double C_gas, double C_l, double layer_mass, double Area_Vol)
     : _cond(cond), _mu_gas(mu_gas), _C_v(C_gas), _C_l(C_l),
       _m_layer(layer_mass), _A_V(Area_Vol)
    { } ;

    void set_state(double T_surf, double T_gas, double P_gas, 
                   RadiationProperties rad, double L_rain, double dt) {

        _Pgas = P_gas ;
        _Tgas = T_gas ;
        _rho_gas = P_gas*_mu_gas/(Rgas*T_gas) ;

        _T0 = T_surf ;
        _dt = dt ;

        _rad = rad ;
        _Lrain = L_rain; 

        // Compute the internal energy and heat capacity modified by radiation
        double sT3 = sigma_rad*_T0*_T0*_T0 ;
        
        _u0 = _m_layer*_C_l*_T0 + dt*(_rad.surface_irradiation() + 3*sT3*_T0 + _Lrain) ;

        _Ctot = _m_layer*_C_l + _dt*4*sT3 ;
    }

    double mass_flux(double T_surf) const {
        
        double Pv = _cond.P_vap(T_surf) ;
        double evap = _cond.P_stick * Pv / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double T_gas, cond ;
        std::tie(T_gas, cond) = compute_gas_T_and_cond_rate(T_surf, evap) ;

        return evap - cond ;
    }

    double surface_cooling_rate(double T_surf) const {
        
        double Pv = _cond.P_vap(T_surf) ;
        double evap = _cond.P_stick * Pv / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double T_gas, cond ;
        std::tie(T_gas, cond) = compute_gas_T_and_cond_rate(T_surf, evap) ;

        double dTdt = (T_surf-_T0)/(_dt + 1e-300) ;

        return _m_layer*_C_l*dTdt - _Lrain + (evap-cond)*(_cond.Lsub - _C_l*T_surf)
            + evap*_C_v*T_surf - cond*_C_v*T_gas ;
    }
    
    double gas_heating_rate(double T_surf) const {
        
        double Pv = _cond.P_vap(T_surf) ;
        double evap = _cond.P_stick * Pv / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double T_gas, cond ;
        std::tie(T_gas, cond) = compute_gas_T_and_cond_rate(T_surf, evap) ;

        return _C_v*(evap*T_surf - cond*T_gas) ;
    }
    
    double operator()(double T_surf) const {

        double Pv = _cond.P_vap(T_surf) ;
        double evap = _cond.P_stick * Pv / std::sqrt(2*pi*Rgas/_mu_gas*T_surf) ;

        double T_gas, cond ;
        std::tie(T_gas, cond) = compute_gas_T_and_cond_rate(T_surf, evap) ;

        double Ctot = _Ctot + (evap*(_C_v - _C_l) + cond*_C_l)*_dt ;

        return T_surf - (_u0 - _cond.Lsub*(evap-cond)*_dt + cond*_C_v*T_gas*_dt) / Ctot ;
    }


    std::tuple<double, double> bracket_solution() const {
        double T_min, T_max ;

        double fac = 2; 

        T_min = _T0 ;
        if (this->operator()(T_min) < 0) {
            T_max = T_min*fac ;
            while(this->operator()(T_max) < 0) {
                T_min = T_max ;
                T_max *= fac ;
            }
        }
        else {
            T_max = T_min ;
            T_min /= fac ;
            while(this->operator()(T_min) > 0) {
                T_max = T_min ;
                T_min /= fac ;
            }
        }

        return std::make_tuple(T_min, T_max) ;
    }

  private:
    // Solve for the gas temperature / pressure given the surface 
    //  temperature
    std::tuple<double,double> 
    compute_gas_T_and_cond_rate(double T_surf, double evap_rate) const {

        double term1 = _dt*_A_V*evap_rate/_rho_gas ;

        double T = (_Tgas + term1*T_surf) / (1 + term1) ;

        double cond0 = _cond.P_stick * std::sqrt(Rgas*T/(2*pi*_mu_gas)) ;

        double rho = _rho_gas * (1 + term1) / (1 + _dt*_A_V*cond0) ;

        return std::make_tuple(T, cond0*rho) ;
    }

    Condensible _cond ;
    double _mu_gas, _C_v, _C_l, _m_layer, _A_V ;

    double _T0, _Tgas, _Pgas, _rho_gas, _u0, _Ctot, _Lrain, _dt ;

    RadiationProperties _rad ;
};



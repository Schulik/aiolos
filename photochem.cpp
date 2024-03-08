/*
 * photochem.cpp
 * 
 * Implementation of the C2Ray scheme, to solve for the correct average simultaneous ionisation and heating rates, according to Mellema+2006 and Friedrich+2012
 */

#define EIGEN_RUNTIME_NO_MALLOC

#include <array>
#include <cassert>

#include "aiolos.h"
#include "brent.h"

/*
Heating and cooling rates according to Black (1981), the physical state of
primordial intergalactic clouds*/

// Atomic data for Hydrogen:
double H_radiative_recombination(double T_e) {
    //return 2.753e-14;
    
    using std::pow;
    if (T_e < 220) T_e = 220.;
    
    double x = 2 * 157807 / T_e;
    // Case B
    return 2.753e-14 * pow(x, 1.500) / pow(1 + pow(x / 2.740, 0.407), 2.242);
    ////       2.753e-14 * pow(2 * 157807 / x, 1.500) / pow(1 + pow(2 * 157807 / x / 2.740, 0.407), 2.242)
    //return 1.e+0 * 2.7e-13 * pow(T_e/1e4, -0.9) ;//* pow(T_e/1e4, -0.9); //MC2009 value
    
    //JOwens code:
    //  rec_coeff =  8.7d-27*tgas**0.5 * t3**(-0.2)*(1+t6**0.7)**(-1) */
}
double H_threebody_recombination(double T_e) {
    using std::pow;
    
    if (T_e < 220) T_e = 220.;

    double x = 2 * 157807 / T_e;
    return 1.e0 * (1.005e-14 / (T_e * T_e * T_e)) * pow(x, -1.089) /
           pow(1 + pow(x / 0.354, 0.874), 1.101);
}
double H_collisional_ionization(double T_e) {
    using std::exp;
    using std::pow;

    if (T_e < 220) T_e = 220.; //return 1e-50; //T_e = 220.;//return 0 ;

    double x = 2 * 157807. / T_e;
    double term = 21.11 * pow(T_e, -1.5) * pow(x, -1.089) /
                  pow(1 + pow(x / 0.354, 0.874), 1.101);
    return 1e0 * term * exp(-0.5 * x);
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
        term = 7.562e-27 * std::pow(Te, 0.42872) ;
    cooling += 1.* nX[1] * term;
    
    // Collisional ionization cooling:
    term = kb * T_HI * H_collisional_ionization(Te);
    cooling += 1. *nX[0] * term;

    // HI Line cooling (Lyman alpha):
    //term = 7.5e-19 * exp(-0.75 * T_HI / Te) / (1 + sqrt(Te / 1e5));
     term = 7.5e-19 * exp(-0.75 * T_HI / Te); //MC2009
    //term = 7.5e-19 * exp(-0.75 * T_HI / 3000.); //MC2009
    cooling += 1. * nX[0] * term;

    // Free-Free:
    term = 1.426e-27 * 1.3 * sqrt(Te) ;
    cooling += 1. * nX[1] * term  ;

    return 1.*cooling;
} // */
//std::vector<double[4]> lines_O;
std::vector<std::array<double, 4>> lines_O;
std::vector<std::array<double, 4>> lines_Op;
std::vector<std::array<double, 4>> lines_Opp;
std::vector<std::array<double, 4>> lines_O3p;
std::vector<std::array<double, 4>> lines_O4p;
std::vector<std::array<double, 4>> lines_C;
std::vector<std::array<double, 4>> lines_Cp;
std::vector<std::array<double, 4>> lines_Cpp;
std::vector<std::array<double, 4>> lines_C3p;
std::vector<std::array<double, 4>> lines_C4p;
void init_line_cooling_data() {
	lines_O.push_back( {6300*angstroem, 1.15685197E-14, 22830.7, 8.61213387e+05});
	lines_Op.push_back( {834*angstroem, 5.78622E-4, 172421.6, 1.32174E15});
	lines_Op.push_back( {2741*angstroem, 3.81198E-13, 58225.3, 4.48777E7});
	lines_Op.push_back( {3727*angstroem, 4.29901E-16, 38575.0, 5.36461E3});
	lines_Op.push_back( {7320*angstroem, 3.76929E-13, 53063.6, 3.11018E7});
	lines_Opp.push_back( {52*micron, 3.13852E-18, 277.682, 2.5493E3});
	lines_Opp.push_back( {5000*angstroem, 3.386678E-14 ,28728.6, 9.66741E5});
	lines_Opp.push_back( {166*angstroem, 6.59979E-10, 86632.4, 1.475889E10});
	lines_Opp.push_back( {83.5*angstroem, 1.75205E3, 172569.7, 5.405937E21});
	lines_O3p.push_back( {25*micron, 2.58649627E-17, 555.66, 3.19064670E+03});
	lines_O3p.push_back( {1404*angstroem,    1.84651335e-08, 102685.994, 8.75179124E+10});
	lines_O4p.push_back( {1218*angstroem, 1.84651335E-08, 118094.57, 7.84917617e+10});

	lines_C.push_back( {1.*angstroem,0, 1., 1.});  //C line cooling currently unclear, dummy line
	lines_Cp.push_back( {157*micron, 1.78292054E-20, 91.2, 1.38779668E+01});
	lines_Cp.push_back( {2326*angstroem, 1.21471023E-10, 61853.9, 1.20977633e+09});
	lines_Cp.push_back( {1334*angstroem, 2.41304508E-03, 107718.1, 3.74008770E+15});
	lines_Cpp.push_back( {1910*angstroem, 3.84223024E-10, 75460.8, 1.31478953E+9});
	lines_Cpp.push_back( {977*angstroem, 1.79050834E-03, 147263.9, 7.17164800E+14});
        lines_C3p.push_back( {1550*angstroem, 5.52131947E-03, 92934.39, 2.86266976e+15});
	lines_C4p.push_back( {6300*angstroem, 1.15685197E-14, 22830.7, 8.61213387e+0});
}

//std::vector<std::array<double, 4>> lines_O;
//lines_O.push_back( {6300*angstr, 1.15685197E-14, 22830.7, 8.61213387e+05});

/*lines_Op = [
    [834*angstr, 5.78622E-4, 172421.6, 1.32174E15],
[2741*angstr, 3.81198E-13, 58225.3, 4.48777E7],
[3727*angstr, 4.29901E-16, 38575.0, 5.36461E3],
[7320*angstr, 3.76929E-13, 53063.6, 3.11018E7]]

lines_Opp = [
[52*micron, 3.13852E-18, 277.682, 2.5493E3],
[5000*angstr, 3.386678E-14 ,28728.6, 9.66741E5],
[166*angstr, 6.59979E-10, 86632.4, 1.475889E10],
[83.5*angstr, 1.75205E3, 172569.7, 5.405937E21]]

lines_C = [
    [0, 0, 0, 1]]

lines_Cp = [
    [157*micron, 1.78292054E-20, 91.2, 1.38779668E+01],
    [2326*angstr, 1.21471023E-10, 61853.9, 1.20977633e+09],
    [1334*angstr, 2.41304508E-03, 107718.1, 3.74008770E+15]]

lines_Cpp = [
    [1910*angstr, 3.84223024E-10, 75460.8, 1.31478953E+9],
    [977*angstr, 1.79050834E-03, 147263.9, 7.17164800E+14]]

*/

double C_cooling(double Te, double ne) {
    double term = 1e-25;
    for (auto & ln : lines_C) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}

double Cp_cooling(double Te, double ne) {    
    double term = 1e-25;
    for (auto & ln : lines_Cp) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}

double Cpp_cooling(double Te, double ne) {
    double term = 1e-25;
    for (auto & ln : lines_Cpp) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}
double C3p_cooling(double Te, double ne) {
    double term = 1e-25;
    for (auto & ln : lines_C3p) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}
double C4p_cooling(double Te, double ne) {
    double term = 1e-25;
    for (auto & ln : lines_C4p) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}
double O_cooling(double Te, double ne) {

    double term = 1e-25;
    for (auto & ln : lines_O) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
	//cout<<ln[1]<<"/"<<ln[2]<<"/"<<ln[3]<<"/"<<ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne))<<endl;
    }
    //for ln in lines_O:
    //double term = 5.5e-24+1.1e-20*std::exp(-30162/Te)*(1.+std::pow(Te/0.75e4, 0.5));    
    return term;
}

double Op_cooling(double Te, double ne) {
    
    double term = 1e-25;
    for (auto & ln : lines_Op) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }
    return term;
}

double Opp_cooling(double Te, double ne) {
    
    double term = 1e-25;
    for (auto & ln : lines_Opp) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }

    return term;
}
double O3p_cooling(double Te, double ne) {
    
    double term = 1e-25;
    for (auto & ln : lines_O3p) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }

    return term;
}
double O4p_cooling(double Te, double ne) {
    
    double term = 1e-25;
    for (auto & ln : lines_O4p) {
        term += ln[1]*std::exp(-ln[2]/Te) / (ne*(1.+ln[3]/ne));
    }

    return term;
}
double h3plus_cooling(double Te) {
    //cout<<"H3+ cooling active!"<<endl;
    double temp = 0.;
    
    if(Te > 900)
        temp = -8.24045e-21 + 3.54583e-23 * Te - 8.66296e-26 * Te*Te + 9.76608e-29* std::pow(Te,3) - 1.61317e-32 * std::pow(Te,4);
    else
        temp = -6.11904e-21 + 4.96694e-23 * Te - 1.443608e-25 *Te*Te + 1.60926e-28* std::pow(Te,3) - 3.87932e-32 * std::pow(Te,4);
        
    if(temp < 0)
        return 0.;
    else
        return temp;
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
            
        for(int b=0; b < num_he_bands; b++) {
            if( tau0[b] * (1 - x_bar) < 1e-10 )
                ion += Gamma0[b]/(nH ) * tau0[b] * (1. - 0.5 * tau0[b] * (1 - x_bar));
            else
                ion += Gamma0[b]/(nH * (1-x_bar) ) * -std::expm1(-tau0[b] * (1 - x_bar));
        }
                
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
    //using Matrix_t =
     //   Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

    C2Ray_HOnly_heating(double GammaH_, double dt_,
                        const std::array<double, 3>& nX,
                        const std::array<double, 3>& TX,
                        double ion_rate, double recomb_rate,
                        const Mat3& collisions)
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

/**
 * Simple by-hand solver using a linearised implicit scheme. Solves for steady-state value by following the time evolution.
 * 
 * @return ion fraction
 */
double get_implicit_photochem() { //double dt, double n1_old, double n2_old, double F, double dx, double kappa, double alpha

    
    int cnt = 0;
    long double dt = 1e-2;
    long double tmax = 1e4;
    long double globalT = 0.;
    long double n1_old = 1.;
    long double n2_old = 0.;
    
    long double F = 1e16;
    long double dx = 1e7;
    long double kappa = 2.2e-18;
    long double alpha = 2.7e-13;
     /*
    long double F = 1e16;
    long double dx = 1e-5;
    long double kappa = 1e-5;
    long double alpha = 1.; */
    
    cout << "Size of double : " << sizeof(double) << " bytes" << endl;
    cout << "Size of long double : " << sizeof(long double) << " bytes" << endl;
     
    char a;
    cout<<"Starting implicit photochem time evolution. Press any button to continue."<<endl;
    cin>>a;
    
    while(globalT < tmax) {
    
        long double dtau = dx*kappa*n1_old;
        
        //long double b1 = n1_old - dt * ( n2_old*n2_old * alpha + F/dx*(1. - std::exp(-dtau)*(1+dtau)) );
        //long double b2 = n2_old + dt * ( n2_old*n2_old * alpha + F/dx*(1. - std::exp(-dtau)*(1+dtau)) );
        
        long double b1 = n1_old - dt * ( n2_old*n2_old * alpha + F/dx*(- std::expm1(-dtau)-dtau*std::exp(-dtau) ) );
        long double b2 = n2_old + dt * ( n2_old*n2_old * alpha + F/dx*(- std::expm1(-dtau)-dtau*std::exp(-dtau)) );
        
        long double dA00 = 1. + dt*F*kappa*std::exp(-dtau);
        long double dA11 = 1. + dt*2.*n2_old*alpha;
        
        long double dA01 =    - dt*F*kappa*std::exp(-dtau);
        long double dA10 =    - dt*2*n2_old*alpha;
        
        long double n2_new = (b2 * dA00 - b1 * dA01) / (dA00*dA11-dA01*dA10);
        long double n1_new = (b1 - n2_new * dA10 )/dA00;
        
        n1_old=n1_new;
        n2_old=n2_new;
        
        cnt++;
        globalT += dt;
    }
    
    cout<<"Done with photochem time evolution. Final nvalues = "<<n1_old<<"/"<<n2_old<<"."<<endl;
    cout<<"1- Final nvalues = "<<1.-n1_old<<"/"<<1.-n2_old<<". Press any button to continue."<<endl;
    cin>>a;
    
    return n2_old/(n1_old+n2_old);
}


/**
 * Simple by-hand solver using a linearised implicit scheme. Only one timestep.
 * 
 * @return ion fraction.
 */
double get_implicit_photochem2(double dt, double n1_old, double n2_old, double F, double dx, double kappa, double alpha) {


        long double dtau = dx*kappa*n1_old;
        
        long double b1 = n1_old - dt * ( n2_old*n2_old * alpha + F/dx*(- std::expm1(-dtau)-dtau*std::exp(-dtau) ) );
        long double b2 = n2_old + dt * ( n2_old*n2_old * alpha + F/dx*(- std::expm1(-dtau)-dtau*std::exp(-dtau)) );
        
        long double dA00 = 1. + dt*F*kappa*std::exp(-dtau);
        long double dA11 = 1. + dt*2.*n2_old*alpha;
        
        long double dA01 =    - dt*F*kappa*std::exp(-dtau);
        long double dA10 =    - dt*2*n2_old*alpha;
        
        long double n2_new = (b2 * dA00 - b1 * dA01) / (dA00*dA11-dA01*dA10);
        long double n1_new = (b1 - n2_new * dA10 )/dA00;
        
        n1_old=n1_new;
        n2_old=n2_new;
    
    return n2_old/(n1_old+n2_old);
}

/**
 * Simple by-hand solver using a linearised implicit scheme. Only one timestep, using ne+np=const.
 * 
 * @return ion fraction.
 */
double get_implicit_photochem3(double dt, double n1_old, double n2_old, double F, double dx, double kappa, double alpha) {
    
    long double ntot = (n1_old + n2_old);
    long double x    = n2_old/ntot;
    long double dtau = dx*kappa*ntot*(1-x);
    
    long double nom   = 1 + dt*F*kappa*std::exp(-dtau) +   dt*ntot*alpha * x;
    long double denom = 1 + dt*F*kappa*std::exp(-dtau) + 2*dt*ntot*alpha * x;
    long double source= dt*F/(ntot*dx) * (-std::expm1(-dtau));
    
    long double xnew = (x*nom + source) / denom;
    
    //long double n2_new = xnew * ntot;
    //long double n1_new = (1-xnew) * ntot;
    
    return xnew; //n2_old/(n1_old+n2_old);
}

/**
 * The actual C2Ray scheme.
 */
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
                    dtaus[b]  = species[0].opacity_twotemp(j, b) * (nX[0] + nX[1]) * mX[0] * dx[j]; 
                    
                    if(debug >= 0 && j >= 1053 && steps <= 2) {
                        //cout<<" PH j= "<<j<<" b == "<<b<<" Gamma0[b] = "<<Gamma0[b]<<" dtau = "<<dtaus[b]<<" opa/nX0+nX1/dx[j] = "<<species[0].opacity_twotemp(j, b)<<"/"<<(nX[0] + nX[1])<<"/"<<dx[j]<<endl;
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
                double x_bar = 1e-10;
                try {
                        //x_bar = get_implicit_photochem3(dt, nX[0], nX[1], Gamma0[1], dx[j], species[0].opacity_twotemp(j, 1) * mX[0], 2.7e-13);
                        //x_bar = get_implicit_photochem(); //dt, n1_old, n2_old, F, dx, kappa, alpha
                        x_bar = root_solver_ion.solve(xmin, 1, ion);
                }
                catch(...) { }

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
                double dn_R = 1.*(ion.R + ion.B*ne)*ne*ne*dt ;
                double dn_I = 1.*(ion.C*ne + ion.photoionization_rate(x_bar))*nH*dt ;

                std::array<double, 3> mom = { 
                    nX[0]*mX[0]*vX[0], nX[1]*mX[1]*vX[1], nX[2]*mX[2]*vX[2] } ;

                double fe = 1/(1 + mX[2]/mX[1]), fp = 1/(1 + mX[1]/mX[2]) ;
                double f = nX_new[0] + dn_I - dn_I*dn_R*(fp/(nX_new[1]+dn_R) + fe/(nX_new[2] + dn_R)) ;
                mom[0] = (mom[0] + dn_R*(mom[1]/(nX_new[1]+dn_R) + mom[2]/(nX_new[2]+dn_R))) / f ;
                mom[1] = (mom[1] + dn_I*mom[0]*fp)/(nX_new[1]+dn_R) ;
                mom[2] = (mom[2] + dn_I*mom[0]*fe)/(nX_new[2]+dn_R) ;

                // Compute the change in kinetic energy:
                double dEk = 0.0*0.5*(mom[0]*mom[0]*nX_new[0] / mX[0] - mX[0]*nX[0]*vX[0]*vX[0] +
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
                    double heatingrate = 0.25 * solar_heating(b)  * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) / dx[j] * dlognu; 
                    GammaH += heatingrate; 
                    dS_band(j,b) = heatingrate;
                    dS_band_special(j,b) = heatingrate;
                }
                
                if(use_rad_fluxes) {
                        
                    Te = TX[2];
                    
                    heating = {0.0, 0.0, GammaH - dEk/(dt+1e-300)}; //    heat.heating_rate(Te);
                    cooling = {0., 0., - nX[2]*HOnly_cooling(nX, Te) }; //
                    
                    //
                    // Enforce no cooling spikes in the boundary cell
                    //
                    if(std::abs(cooling[2]) > heating[2])
                        cooling[2] = -heating[2];
                    
                } else {
                    
                    // Collisional heat exchange:
                    //if(!use_rad_fluxes) {
                    fill_alpha_basis_arrays(j);
                    compute_collisional_heat_exchange_matrix(j);
                    //}
                    
                    // Solve for radiative cooling implicitly
                    //Comment back in for 3 species used
                    /*
                    C2Ray_HOnly_heating heat(GammaH - dEk/(dt+1e-300), dt, nX_new, TX,
                           (ion.R + ion.B*ne)*ne*ne, (ion.C*ne + ion.photoionization_rate(x_bar))*nH,
                           friction_coefficients);
                    */
                    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> dummy;
                    C2Ray_HOnly_heating heat(GammaH - dEk/(dt+1e-300), dt, nX_new, TX,
                           (ion.R + ion.B*ne)*ne*ne, (ion.C*ne + ion.photoionization_rate(x_bar))*nH,
                           dummy );
                    
                    
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
                    //dS_band(j,0)     += heating[s];
                }
                
                if(j==num_cells-10e99) {
                        cout<<"Pos 1.1 dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                
                
                //char a;
                //cin>>a;
                
                double tau = total_opacity(j,0)*(x_i12[j+1]-x_i12[j])*1e2;
                double red = (1.+tau*tau);
                
                if( C_idx!=-1 && e_idx!=-1 ) { species[C_idx].dG(j)     -=  C_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( Cp_idx!=-1 && e_idx!=-1 ) { species[Cp_idx].dG(j)   -=  Cp_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( Cpp_idx!=-1 && e_idx!=-1 ) { species[Cpp_idx].dG(j) -=  Cpp_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( C3p_idx!=-1 && e_idx!=-1 ) { species[C3p_idx].dG(j) -=  C3p_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( C4p_idx!=-1 && e_idx!=-1 ) { species[C4p_idx].dG(j) -=  C4p_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( O_idx!=-1 && e_idx!=-1 ) { species[O_idx].dG(j)     -=  O_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( Op_idx!=-1 && e_idx!=-1 ) { species[Op_idx].dG(j)   -=  Op_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( Opp_idx!=-1 && e_idx!=-1 ) { species[Opp_idx].dG(j) -=  Opp_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( O3p_idx!=-1 && e_idx!=-1 ) { species[O3p_idx].dG(j) -=  O3p_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                if( O4p_idx!=-1 && e_idx!=-1 ) { species[O4p_idx].dG(j) -=  O4p_cooling(species[e_idx].prim[j].temperature, ne)/red; }
                
                //Correction in heating function for non-ionising radiation. 
                double adddS[3] = {0.,0.,0.};
                for (int b = 0; b < num_he_bands; b++) {
                    
                    double temp_nonion_cell_optd = 0.;
                    double dtauplus = dtaus[b];
                    //double temp_total_cell_optd = 0.;
                    double GammaTotal = 0.;
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
                        GammaHnofactor = 0.25 * solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) / dx[j]; //dtaus[b] contains only ionising opacities
                        GammaTotal         = 0.25 * solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtauplus)) / dx[j]; //dtauplus contains all opacities
                        radial_optical_depth_twotemp(j,b) = radial_optical_depth_twotemp(j+1,b) + dtauplus;
                        cell_optical_depth_highenergy(j,b)= radial_optical_depth_twotemp(j+1,b) + dtauplus;
                        S_band(j,b)                       = solar_heating(b) * std::exp(-radial_optical_depth_twotemp(j+1,b)); 
                    }
                        
                    else {
                        GammaHnofactor = 0.25 * solar_heating(b) * dtaus[b];
                        GammaTotal         = 0.25 * solar_heating(b) * dtauplus;
                        radial_optical_depth_twotemp(j,b) = dtauplus;
                        cell_optical_depth_highenergy(j,b)= dtauplus;
                        S_band(j,b)                       = solar_heating(b);
                    }
                    
                    //cout<<"PCHEM pos 4, GammaH  = "<<GammaHnofactor<<endl;
                    //
                    // Now recompute heating with higher optical depth, neglecting ionization and ionization potential
                    //
                    
                    dS_band(j,b)  += 0.*(GammaTotal - 1.*GammaHnofactor);
                    
                    if(debug >= 2)
                        cout<<"PCHEM pos 6, dGamma  = "<<dS_band(j,b)<<endl;
                    //Distribute GammaResidual to other species according to their optical depth contribution
                    
                    for(int s=3; s<num_species; s++) {
                        //species[s].dS(j) += GammaResidual * species[s].opacity_twotemp(j,b) * species[s].u[j].u1 * dx[j] / temp_total_cell_optd;
                        species[s].dS(j)    += 1.* dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                        adddS[s]            += 1.* dS_band(j,b) * species[s].fraction_total_solar_opacity(j,b);
                    
                        //cout<<" in highe-correction  dS_band(j,b) = "<< dS_band(j,b)<<" to species["<<s<<"]."<<endl;
                    }
                    //char d;
                    //cin>>d;
                    
                    if(dS_band(j,b) > 1e15) {
                            for(int s=3; s<num_species; s++) {
                                    cout<<" PHOTOCHEM j/s/dS = "<<j<<"/"<<s<<"/"<<species[s].dS(j)<<endl;
                            }
                    } 
                    
                    //cout<<"PCHEM pos 7, radial_tau  = "<<radial_optical_depth_twotemp(j,b)<<endl;
                    
                    cell_optical_depth_twotemp(j,b)   = dtauplus;
                    total_opacity_twotemp(j,b)        = dtauplus/dx[j];
                    
                    //No factor 0.25 here, as we want to see the unaveraged radiation we put in, consistent with low-energy bands
                    //dS_band(j,b) = 0.25 * solar_heating(b) * (1 - 13.6 * ev_to_K * kb / photon_energies[b]) * std::exp(-radial_optical_depth_twotemp(j+1,b)) * (-std::expm1(-dtaus[b])) / dx[j];
                                
                }
                
                if(j==num_cells-10e99) {
                        cout<<"Pos 1.2 dS_UV = "<<dS_band(num_cells-10,0)<<endl;
                }
                
                //Previous loop was to take into account the difference between non-ionising and ionising highenergy absorption, but it will report 0 heating in the diagnostic file, as dS_band is zero for cases of no difference. The actual highenergy heating has been written directly into dS. So we now just reset dS_band so that its documentation is consistent with dS. 
                for (int b = 0; b<num_he_bands;b++) {
                        //if(BAND_IS_HIGHENERGY[b]) {
                        //dS_band(j,b) = 0.;
                        for (int s = 0; s < 3; s++) {
                            //dS_band(j,b)     += heating[s];
                            
                        }
                }
                //cout<<"adding"<<dS_band(j,0)<<endl;
                
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

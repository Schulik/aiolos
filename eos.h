/**
 * eos.h
 * 
 * Declarations on the equation of state. Assumes a global equation of state. 
 * Spatially changing EOS, or adiabatic EOS with spatially varying adiabatic index would require modifications.
 */
#ifndef _AIOLOS_EOS_H_
#define _AIOLOS_EOS_H_

struct AOS ;
struct AOS_prim ;

class EOS_Base {

  public:
    virtual void compute_primitive(const AOS*, AOS_prim*, int) const = 0;
    virtual void compute_conserved(const AOS_prim*, AOS*, int) const = 0 ;
    virtual void compute_auxillary(AOS_prim*, int) const = 0;
    virtual void update_p_from_eint(AOS_prim*, int) const =0;
    virtual void update_eint_from_T(AOS_prim*, int) const =0;
    virtual void get_p_over_rho_analytic(const double*, double*) const = 0;
    virtual void get_p_from_rhoT(const double*,const double*, double*) const = 0;
    
    // Base classes must have virtual destructors
    virtual ~EOS_Base(){} ; 
} ;

class IdealGas_EOS final: public EOS_Base {

  public:
    IdealGas_EOS(double dof, double cv, double mass) 
      : _gamma_m1(2./dof), _cv(cv), _mass(mass)
    { } ;

    void compute_primitive(const AOS*, AOS_prim*, int) const ;
    void compute_conserved(const AOS_prim*, AOS*, int) const ;
    void compute_auxillary(AOS_prim*, int) const ;
    void update_p_from_eint(AOS_prim*, int) const ;
    void update_eint_from_T(AOS_prim*, int) const ;
    void get_p_over_rho_analytic(const double*, double*) const ;
    void get_p_from_rhoT(const double*,const double*, double*) const;
    
  private:
    double _gamma_m1, _cv, _mass ;
} ;

class Polytropic_EOS final: public EOS_Base {

  public:
    Polytropic_EOS(double gamma, double cv, double mass, double K0=1) 
      : _gamma(gamma), _cv(cv), _mass(mass), _K0(K0) 
    { } ;

    void compute_primitive(const AOS*, AOS_prim*, int) const ;
    void compute_conserved(const AOS_prim*, AOS*, int) const ;
    void compute_auxillary(AOS_prim*, int) const ;
    void update_p_from_eint(AOS_prim*, int) const ;
    void update_eint_from_T(AOS_prim*, int) const ;
    void get_p_over_rho_analytic(const double*, double*) const ;
    void get_p_from_rhoT(const double*,const double*, double*) const;
    
  private:
    double _gamma, _cv, _mass, _K0 ;
} ;


#endif//_AIOLOS_EOS_H_

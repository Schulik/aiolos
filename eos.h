
#ifndef _AIOLOS_EOS_H_
#define _AIOLOS_EOS_H_

struct AOS ;
struct AOS_prim ;

class EOS_Base {

  public:
    virtual void compute_primitive(const AOS*, AOS_prim*, int) const = 0;
    virtual void compute_conserved(const AOS_prim*, AOS*, int) const = 0 ;
    virtual void compute_auxillary(AOS_prim*, int) const = 0;
    virtual void get_p_over_rho_analytic(const double*, double*) const = 0;
    virtual void get_p_from_rhoT(const double*,const double*, double*) const = 0;
    
    // Base classes must have virtual destructors
    virtual ~EOS_Base(){} ; 
} ;

class Adiabatic_EOS final: public EOS_Base {

  public:
    Adiabatic_EOS(double gamma, double cv) 
      : _gamma(gamma), _cv(cv) 
    { } ;

    void compute_primitive(const AOS*, AOS_prim*, int) const ;
    void compute_conserved(const AOS_prim*, AOS*, int) const ;
    void compute_auxillary(AOS_prim*, int) const ;
    void get_p_over_rho_analytic(const double*, double*) const ;
    void get_p_from_rhoT(const double*,const double*, double*) const;
    
  private:
    double _gamma, _cv ;
} ;

class Polytropic_EOS final: public EOS_Base {

  public:
    Polytropic_EOS(double gamma, double cv, double K0=1) 
      : _gamma(gamma), _cv(cv), _K0(K0)
    { } ;

    void compute_primitive(const AOS*, AOS_prim*, int) const ;
    void compute_conserved(const AOS_prim*, AOS*, int) const ;
    void compute_auxillary(AOS_prim*, int) const ;
    void get_p_over_rho_analytic(const double*, double*) const ;
    void get_p_from_rhoT(const double*,const double*, double*) const;
    
  private:
    double _gamma, _cv, _K0 ;
} ;


#endif//_AIOLOS_EOS_H_

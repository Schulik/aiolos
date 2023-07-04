#!/usr/bin/python3
import sys, os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params

class Riemann(object):
    '''Exact Solution of a 1D Riemann problem for an ideal gas with a fixed
    adiabatic index.

    args:
       L : Left State, dict. contains 'rho', 'u' and 'P'
       R : Right State, as L.
       
       gamma : default 5/3. Adiabatic index
    '''
    def __init__(self, L, R, gamma=5/3.):
        self._L = L
        self._R = R

        self._L['a'] = (gamma * L['P'] / L['rho'])**0.5
        self._R['a'] = (gamma * R['P'] / R['rho'])**0.5
        
        self._g  = gamma
        self._gp1 = gamma + 1.
        self._gm1 = gamma - 1.
        self._mu  = (gamma-1.) / (gamma+1.)
        self._alp = 0.5 * (gamma+1.) / gamma
        self._alm = 0.5 * (gamma-1.) / gamma

        # Compute the star values
        self._pstar = self._ps()
        self._vstar = self._vs(self._pstar)

        # Construct the solutions between the waves
        if (self._pstar > L['P']):
            # Left shock
            x = self._pstar / L['P']
            self._Sl = L['u'] - L['a'] * (self._alp*x + self._alm)**0.5
            self._rhosl = L['rho'] * (x + self._mu)/ (1 + x*self._mu)
            self._SlT = None
        else:
            # Left rarefaction
            x = self._pstar / L['P']
            self._rhosl = L['rho'] * x**(1./self._g)
            self._asl = L['a'] * x**self._alm
            self._Sl  = L['u'] - L['a'] 
            self._SlT = self._vstar - self._asl

        # Construct the solutions between the waves
        if (self._pstar > R['P']):
            # Right shock
            x = self._pstar / R['P']
            self._Sr = R['u'] + R['a'] * (self._alp*x + self._alm)**0.5
            self._rhosr = R['rho'] * (x + self._mu)/ (1 + x*self._mu)
            self._SrT = None
        else:
            # Right rarefaction
            x = self._pstar / R['P']
            self._rhosr = R['rho'] * x**(1./self._g)
            self._asr = R['a'] * x**self._alm
            self._Sr  = R['u'] + R['a'] 
            self._SrT = self._vstar + self._asr

        self.density  = np.vectorize(self._rho)
        self.velocity = np.vectorize(self._vel)
        self.pressure = np.vectorize(self._pres)

    def internal_energy(self,x, t):
        '''Internal energy per unit mass'''
        return self.pressure(x,t) / (self._gm1 * self.density(x,t))

    # Base functional form of the density in rarefactions
    def _rare_L(self, vt):
        L = self._L
        return (2./ self._gp1 + self._mu*(L['u'] - vt) / L['a'])**(2./self._gm1)

    def _rare_R(self, vt):
        R = self._R
        return (2./ self._gp1 - self._mu*(R['u'] - vt) / R['a'])**(2./self._gm1)

    def _rho(self, x, t):
        '''Density'''
        vt = x / t
        if vt < self._vstar:
            L = self._L
            if vt < self._Sl:
                return L['rho']
            elif self._pstar > L['P'] or vt > self._SlT:
                return self._rhosl
            else:
                return L['rho']*self._rare_L(vt)
        else:
            R = self._R 
            if vt > self._Sr:
                return R['rho']
            elif self._pstar > R['P'] or vt < self._SrT:
                return self._rhosr
            else:
                return R['rho']*self._rare_R(vt)

    def _vel(self, x, t):
        '''Velocity'''
        vt = x / t
        if vt < self._vstar:
            L = self._L
            if vt < self._Sl:
                return L['u']
            elif self._pstar > L['P'] or vt > self._SlT:
                return self._vstar
            else:
                return (2./ self._gp1)*(+L['a'] + 0.5*self._gm1*L['u'] + vt)
        else:
            R = self._R 
            if vt > self._Sr:
                return R['u']
            elif self._pstar > R['P'] or vt < self._SrT:
                return self._vstar
            else:
                return (2./ self._gp1)*(-R['a'] + 0.5*self._gm1*R['u'] + vt)

    def _pres(self, x, t):
        '''Pressure'''
        vt = x / t
        if vt < self._vstar:
            L = self._L
            if vt < self._Sl:
                return L['P']
            elif self._pstar > L['P'] or vt > self._SlT:
                return self._pstar
            else:
                return L['P'] * self._rare_L(vt)**self._g
        else:
            R = self._R 
            if vt > self._Sr:
                return R['P']
            elif self._pstar > R['P'] or vt < self._SrT:
                return self._pstar
            else:
                return R['P'] * self._rare_R(vt)**self._g


    def _fp(self, p, P):
        '''pressure functionm for determing p star'''
        if p > P['P']: # Shock
            A = 2. / (self._gp1 * P['rho'])
            B = self._mu * P['P']
            return (p - P['P']) * (A / (p+B))**0.5 
        else: # Rarefaction
            x = p / P['P']
            return 2.*P['a'] * (x**self._alm - 1.) / self._gm1 
        
    def _ps(self, tol=1e-7):
        '''Compute central pressure'''
        du = self._R['u'] - self._L['u']
        
        # p* is given by the root here:
        def f(p):
            return self._fp(p,  self._L) + self._fp(p, self._R) + du
        
        # Find the bounds
        pmin = min(self._L['P'], self._R['P'])
        pmax = max(self._L['P'], self._R['P'])

        fmin = f(pmin)
        fmax = f(pmax)
        if fmin == 0:
            return pmin
        elif fmax == 0:
            return pmax 
        
        
        if fmin > 0 and fmax > 0:
            pmin, pmax = 0, pmin
        elif fmin < 0 and fmax > 0:
            pass
        else:
            pmin = pmax
            pmax *= 2 
            while f(pmax) < 0:
                pmin = pmax     
                pmax *= 2

        # Use bisection to find the root
        while (pmax - pmin) > tol:
            pm = 0.5*(pmin + pmax)
            fp = f(pm)
            if fp == 0:
                return pm
            elif fp > 0:
                pmax = pm
            else:
                pmin = pm

        return 0.5*(pmax + pmin)
                
    def _vs(self, ps):
        '''Central Celocity'''
        return 0.5*(self._L['u'] + self._R['u'] + 
                    self._fp(ps, self._R) - self._fp(ps, self._L))


def setup_riemann_solver(param_file):
    params = load_aiolos_params(param_file)
    
    UL = list(map(float, [params['PARI_INIT_DATA_U1L'],
                          params['PARI_INIT_DATA_U2L'],
                          params['PARI_INIT_DATA_U3L']]))
    UR = list(map(float, [params['PARI_INIT_DATA_U1R'],
                          params['PARI_INIT_DATA_U2R'],
                          params['PARI_INIT_DATA_U3R']]))

    GAMMA= float(params['PARI_GAMMA'])

    L = { 'rho' : UL[0], 'u': UL[1], 'P' : UL[2] }
    R = { 'rho' : UR[0], 'u': UR[1], 'P' : UR[2] }
                                              
    return Riemann(L, R, GAMMA)                                             

def plot_riemann_solution(par_num):
    param_file = 'shock_tube{}.par'.format(par_num)
    snap_file = 'output_shock_tube{}_H2_t-1.dat'.format(par_num)

    #Get the Riemann problem
    exact = setup_riemann_solver(param_file)

    data = load_aiolos_snap(snap_file)

    params = load_aiolos_params(param_file)
    x0 = float(params['PARI_INIT_SHOCK_MID'])
    t  = float(params['PARI_TIME_TMAX'])

    x = data['x']
    xa = x - x0
   
    import matplotlib.pyplot as plt

    f, axes = plt.subplots(2,2)
    f.suptitle('Riemann Problem {}'.format(par_num))

    axes[0][0].plot(x, data['density'], '+', c='k')
    axes[0][0].plot(x, exact.density(xa, t), c='k')
    axes[0][0].set_ylabel('density')   
    axes[0][0].xaxis.set_ticklabels([])


    axes[0][1].plot(x, data['velocity'], '+', c='k')
    axes[0][1].plot(x, exact.velocity(xa, t), c='k')
    axes[0][1].set_ylabel('velocity')
    axes[0][1].xaxis.set_ticklabels([])

    axes[1][0].plot(x, data['pressure'], '+', c='k')
    axes[1][0].plot(x, exact.pressure(xa, t), c='k')
    axes[1][0].set_ylabel('pressure')
    axes[1][0].set_xlabel('x')

    u= (data['energy'] - 0.5*data['momentum']*data['velocity'])/data['density']
    axes[1][1].plot(x, u, '+', c='k')
    axes[1][1].plot(x, exact.internal_energy(xa, t), c='k')
    axes[1][1].set_ylabel('internal energy')
    axes[1][1].set_xlabel('x')

    f.subplots_adjust(top=0.9, bottom=0.11,
                      left=0.095, right=0.965,
                      hspace=0.05, wspace=0.235)
    return f
    
def check_riemann_problem(par_num, L1_target=0,
                          make_plots=False, raise_error=False):
    param_file = 'shock_tube{}.par'.format(par_num)
    snap_file = 'output_shock_tube{}_H2_t-1.dat'.format(par_num)

    #Get the Riemann problem
    exact = setup_riemann_solver(param_file)

    data = load_aiolos_snap(snap_file)

    params = load_aiolos_params(param_file)
    x0 = float(params['PARI_INIT_SHOCK_MID'])
    t  = float(params['PARI_TIME_TMAX'])

    L1 = np.abs(data['density'] - exact.density(data['x']-x0, t)).mean()


    if raise_error:
        np.testing.assert_array_less(L1, L1_target)
    
    
    if L1 <= L1_target:
        print('Riemann test {} L1 check passed'.format(par_num))
    else:
        print("Riemann test {} L1 check failed. ".format(par_num) + 
              "L1={}, expected={}".format(L1,L1_target))

    if make_plots:
        os.makedirs('plots', exist_ok=True)
        fig_name = os.path.join('plots', 'shock_tube{}.png'.format(par_num))
        plot_riemann_solution(par_num).savefig(fig_name)
        

def test_shock_tube(make_plots=False):
    L1_errors = [np.nan, 0.016, 0.033, 0.25,
                 1.06, 0.078, 0.0040, 0.0035]
    for i in range(1, 8):
        check_riemann_problem(i, L1_errors[i], make_plots, not make_plots)
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser("Check shock tube test results")
    parser.add_argument("-p", "--make_plots", default=None,
                        action="store_true",
                        help="Make plots of the results")
    args = parser.parse_args()

    try:
        test_shock_tube(args.make_plots)
    except AssertionError:
        pass

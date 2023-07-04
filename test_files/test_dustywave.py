import os
import numpy as np

from load_aiolos import load_aiolos_snap, load_aiolos_params

class _DustyWaveSolution(object):
    def __init__(self, t, k,
                 rho_g, rho_d, P0, drho_g, drho_d, v_g, v_d, dP, v_f):
        self._t = t
        self._k = k
        self._rho_g = rho_g
        self._rho_d = rho_d
        self._drho_g = drho_g
        self._drho_d = drho_d
        self._v_g = v_g
        self._v_d = v_d
        self._P0 = P0
        self._dP = dP
        self._v_f = v_f


    def v_gas(self, x):
        return (self._v_f + self._v_g * np.exp(1j*self._k*x)).real
    def v_dust(self, x):
        return (self._v_f + self._v_d * np.exp(1j*self._k*x)).real

    def rho_gas(self, x):
        return (self._rho_g + self._drho_g * np.exp(1j*self._k*x)).real
    def rho_dust(self, x):
        return (self._rho_d + self._drho_d * np.exp(1j*self._k*x)).real

    def P(self, x):
        return (self._P0 + self._dP * np.exp(1j*self._k*x)).real

    @property
    def time(self):
        return self._t

class DustyWaveSolver(object):

    def __init__(self, rho_g=1., rho_d=1., cs=1., K=1., delta=1e-6, vf=0,
                 GAMMA=1.4, wavelength=1., feedback=True):
        
        self.GAMMA = GAMMA

        self.vf = vf
        self.rho_g = rho_g
        self.rho_d = rho_d
        self.cs = cs
        self.K = K
        self.delta = delta
        self.wavelength=wavelength
        self.feedback=feedback

    def _solve_system(self, times):
        '''Solve the dusty-wave problem up to specified times'''
        k = 2*np.pi / self.wavelength
        cs2 = self.cs**2
        gammaP0 = cs2 * self.rho_g 

        ts_inv = self.K * self.rho_g / self.rho_d
        if self.feedback:
            Kg = (self.rho_d/self.rho_g) * ts_inv
            Kd = ts_inv
        else:
            Kg = 0
            Kd = ts_inv

        # Construct the jacobian of the linearized system
        jac = np.zeros([5,5], dtype='c8')
        jac[0,:] = [-1j*k*self.vf, -1j*k*self.rho_g, 0, 0, 0]
        jac[1,:] = [0, -Kg - 1j*k*self.vf, 0,  Kg, -1j*k/self.rho_g]
        jac[2,:] = [0, 0, -1j*k*self.vf, -1j*k*self.rho_d, 0]
        jac[3,:] = [0,  Kd, 0, -Kd - 1j*k*self.vf,      0         ]
        jac[4,:] = [0, - 1j*k*gammaP0, 0, 0, -1j*k*self.vf]

        # Do the right going part of the wave (we can get the left going wave
        # for free by taking the real part).
        e = -1j * self.delta

        IC = np.array([self.rho_g*e, self.cs*e, 
                       self.rho_d*e, self.cs*e,
                       gammaP0*e],
                      dtype='c8')

        # Diagonalize the equations
        l, U = np.linalg.eig(jac)
        u0 = np.linalg.solve(U, IC)

        sol = []
        for ti in times:
            # Use the analytical solution and project back to 
            # the primitive variables
            sol.append(np.dot(U, u0 * np.exp(l*ti)))

        return np.array(sol)

    def __call__(self, time):
        '''Solve the dustwave problem at a given time or list of times'''
        try:
            iter(time)
        except TypeError:
            time = [time,]
        time = sorted(time)
        sol = self._solve_system(list(time))

        def _create_solution(t, drho_g, v_g, drho_d, v_d, dP):
            return _DustyWaveSolution(t, 2*np.pi/self.wavelength,
                                      self.rho_g, self.rho_d, 
                                      self.rho_g*self.cs**2/self.GAMMA,
                                      drho_g, drho_d, v_g, v_d, dP, self.vf)

        sol = [_create_solution(t, *s) for t, s in zip(time,sol)]
        if len(sol) == 1:
            return sol[0]
        else:
            return sol

def solve_dustywave(param_file):
    params = load_aiolos_params(param_file)
    dw = DustyWaveSolver(K = float(params["PARI_ALPHA_COLL"]))

    time = float(params['PARI_TIME_TMAX'])
    
    return dw(time)

def plot_dustywave(name):
    param_file = name + '.par'

    gas  = load_aiolos_snap('output_' + name + '_H2_t-1.dat')
    dust = load_aiolos_snap('output_' + name + '_Dust_t-1.dat')
    
    sol = solve_dustywave(param_file)

    x = gas['x']
    
    import matplotlib.pyplot as plt

    plt.figure()

    c = plt.plot(x, gas['velocity'], '+', label='gas')
    plt.plot(x, sol.v_gas(x), ls='-', c=c[0].get_color())
    c = plt.plot(x, dust['velocity'], 'x', label='dust')
    plt.plot(x, sol.v_dust(x), ls='--', c=c[0].get_color())
    
    plt.xlabel('x')
    plt.ylabel('velocity')
    plt.legend()
    plt.title(name)
    plt.savefig(os.path.join('plots',name+'.png'))

def check_dustywave_problem(name, L1s=None, make_plots=False):

    param_file = name + '.par'

    gas  = load_aiolos_snap('output_' + name + '_H2_t-1.dat')
    dust = load_aiolos_snap('output_' + name + '_Dust_t-1.dat')
    
    sol = solve_dustywave(param_file)

    x = gas['x']
    
    L1 = [
        np.abs(gas['velocity'] - sol.v_gas(x)).mean(),
        np.abs(dust['velocity'] - sol.v_dust(x)).mean()
        ]

    if L1s is not None:
        passed = True
        for l1, target in zip(L1, L1s):
            passed &= l1 <= target

        if passed:
            print('Test {} L1 check passed'.format(name))
        else:
            print('Test {} L1 checked failed:'.format(name))
            print('\tL1={}, target={}'.format(L1, L1s))
        
        np.testing.assert_array_less(L1, L1s)

    else:
         print('Test {} L1 values:'.format(name))
         print('\tL1={}'.format(L1))

    if make_plots:
         plot_dustywave(name)
    
def test_dustywave(make_plots=False):

    check_dustywave_problem("dustywave_nonstiff", [1.25e-8, 1.34e-8],
                            make_plots=make_plots)
    check_dustywave_problem("dustywave_stiff", [4.78e-8, 4.78e-8],
                            make_plots=make_plots)

if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser("Check dusty wave test results")
    parser.add_argument("-p", "--make_plots", default=None,
                        action="store_true",
                        help="Make plots of the results")
    args = parser.parse_args()

    try:
        test_dustywave (args.make_plots)
    except AssertionError:
        pass

import numpy as np

def load_aiolos_snap(filename):
    dtype = [('x',          'f8'),
             ('density',    'f8'),
             ('momentum',   'f8'),
             ('energy',     'f8'),
             ('pressure',   'f8'),
             ('velocity',    'f8'),
             ('temperature', 'f8'),
             ('soundspeed', 'f8'),]
    cols=[0,1,2,3,10,11,12,15]
    
    data = np.genfromtxt(filename, dtype=dtype, usecols=cols)
    return data

def load_aiolos_params(filename):
    params = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if len(line) == 0 or line[0] == '#':
                continue

            key, val = line.split()[:2]
            
            params[key.strip()] = val.strip()
    return params

def load_aiolos_diag(filename):
    # Work out the number of bands:
    with open(filename) as f:
        line = f.readline()
        cols = len(line.strip().split())
        bands = (cols - 6) // 10

    # Construct the dtype:
    def add_cols(name, dt='f8'):
        return [ (name + str(i), dt) for i in range(bands)]
    def add_kappa_cols():
        cols = []
        for i in range(bands):
            cols += [('kappa_R' + str(i), 'f8')]
            cols += [('kappa_P' + str(i), 'f8')]
            cols += [('kappa_P(Tsun)' + str(i), 'f8')]
        return cols
        
    dtype = [('x',    'f8'),
             ('Jtot', 'f8'),
             ('Stot', 'f8')]
    dtype += add_cols('J')
    dtype += add_cols('S')
    dtype += add_cols('dS')
    dtype += add_cols('rho_kappaR')
    dtype += add_cols('tau_cell')
    dtype += add_cols('tau_radial')
    dtype += add_kappa_cols()
    dtype += [('T_rad', 'f8'),
              ('T_gas', 'f8'),
              ('pressure', 'f8')]
    dtype += add_cols('Flux')

    # Load the data
    return np.genfromtxt(filename, dtype=dtype)
    
class Species:
    def __init__(self, name, mass, dof, is_dust):
        self.name = name
        self.mass = mass
        self.dof = dof
        self.gamma = (dof + 2.)/dof
        self.dust_like = is_dust 

def load_aiolos_species(filename):
    species = {}
    with open(filename, 'r') as f:
        for l in f:
            if l.startswith('@'):
                data = l[1:].split()
                species[int(data[0].strip())] = Species(data[1].strip(),
                                                    float(data[2].strip()),
                                                    float(data[3].strip()),
                                                    int(data[7].strip()))

    return species
            
    

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    for snap in sys.argv[1:]:
        data = load_aiolos_snap(snap)
        plt.plot(data['x'], data['density'],
                 marker='+', ls='', label=snap)
        

    plt.xlabel('x')
    plt.ylabel('Mach number')
    plt.legend()

    plt.show()

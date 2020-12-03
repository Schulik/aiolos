import numpy as np

def load_aiolos_snap(filename):
    dtype = [('x',          'f8'),
             ('density',    'f8'),
             ('momentum',   'f8'),
             ('energy',     'f8'),
             ('pressure',   'f8'),
             ('velocity',    'f8'),
             ('temperature', 'f8')]
    cols=[0,1,2,3,10,11,12]
    
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
        plt.plot(data['x'], data['velocity'], marker='+', ls='', label=snap)

    plt.xlabel('x')
    plt.ylabel('density')
    plt.legend()

    plt.show()

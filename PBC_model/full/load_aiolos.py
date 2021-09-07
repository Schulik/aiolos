import numpy as np

def load_aiolos_snap(filename):
    dtype = [('x',          'f8'),
             ('density',    'f8'),
             ('momentum',   'f8'),
             ('energy',     'f8'),
             ('pressure',   'f8'),
             ('velocity',    'f8'),
             ('temperature', 'f8'),
             ('soundspeed',  'f8'),
             ('potential',   'f8'),
             ('cooling',     'f8'),
             ('heating',     'f8')]
    cols=[0,1,2,3,10,11,12,15,19,21,22]
    
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

def load_aiolos_diag(filename, bands=None, species=None):
    if bands is None or species is None:
        # Work out the number of bands:
        with open(filename) as f:
            line = f.readline()
            cols = len(line.strip().split())

        def c(b, s): return 5 + 7*b + 3*b*s + s
        res = []
        b = 1
        while True:
            s = 1
            if c(b,s) > cols: break

            while c(b,s) < cols:
                s += 1
            if c(b,s) == cols:
                res.append((b,s))

            b += 1
        if len(res) == 0:
            raise ValueError("Error when finding bands and species")
        if len(res) > 1:
            for b, s in res:
                if b == bands or s == species:
                    bands = b
                    species = s
                    break
            else:
                raise ValueError("Error determining the number of bands and"
                                 "species, help out by specifying them")
        else:
            bands, species = res[0]
            
    # Construct the dtype:
    def add_cols(name, dt='f8'):
        return [ (name + str(i), dt) for i in range(bands)]
    def add_kappa_cols():
        cols = []
        for i in range(bands):
            for j in range(species):
                cols += [('kappa_P(Tsun)' + str(i) + '_' + str(j), 'f8')]
        for i in range(bands):
            for j in range(species):
                cols += [('kappa_R' + str(i) + '_' + str(j), 'f8')]
            for j in range(species):
                cols += [('kappa_P' + str(i) + '_' + str(j), 'f8')]
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
    dtype += [('T_rad', 'f8')]
    dtype += [ ('T_gas' + str(i), 'f8') for i in range(species)]
    dtype += [('pressure' , 'f8')]
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

    f, ax = plt.subplots(4,1, sharex=True)

    snaps = [s for s in sys.argv[1:] if s[0] == 'o'] + [s for s in sys.argv[1:] if s[0] != 'o']
    
    for snap in snaps:
        if snap.startswith('output'): 
            data = load_aiolos_snap(snap)
            if 'Gas' in snap:
                m='+'
            else:
                m ='x'

            ax[0].loglog(data['x'], data['density'],
                        marker=m, ls='', label=snap)


            ax[1].semilogx(data['x'], data['temperature'], ls='',
                        marker=m, label=snap)
            #ax[1].semilogx(data['x'], data['pressure']/data['density']**1.29, ls='',
            #            marker=m, label=snap)

            ax[2].semilogx(data['x'], data['velocity'], ls='',
                marker=m, label=snap)

            l,= ax[3].loglog(data['x'], data['heating'], ls='-',
                label=snap)
            ax[3].semilogx(data['x'], np.abs(data['cooling']), ls='--',
                c=l.get_color())
            #ax[2].semilogx(data['x'], data['soundspeed'], ls='--',
            #    marker=m)
        elif snap.startswith('diagnostic'):
            data = load_aiolos_diag(snap)
            ax[1].semilogx(data['x'], data['T_rad'], ls='',
                        marker='^', label=snap)
            #ax[3].semilogx(data['x'], np.abs(data['dS0']), marker='^', ls='')

    ax[-1].set_xlabel('r')
    ax[0].set_ylabel('density')
    ax[1].set_ylabel('temperature')
    ax[2].set_ylabel('velocity')
    ax[3].set_ylabel('heating / cooling')


    ax[1].legend(ncol=3)
    ax[2].legend(ncol=2)

    plt.show()

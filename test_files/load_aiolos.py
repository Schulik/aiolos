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

#!/home/devon/anaconda2/bin/python

import numpy as np
import sys
from scipy.ndimage.filters import gaussian_filter
sys.path.insert(0, '.')
import hydroPlayground as hp

# global stuff
N = 128
L = 1.0
dim = 2

# used to initialize the grid cells
def init_grid(grid):

    # set the flat background field
    grid['rho'] = 0.001 
    grid['mom'] = 0.0 
    grid['etot'] = 1.0 

    #grid['rho'][N/2-10+2:N/2+10+2,N/2+20:N/2+40] = 10

    # add some blobbiness
    for sample in np.random.randint(0, high=N, size=(80, dim)):
        grid['rho'][tuple(sample)] = 10.0 
    gaussian_filter(grid['rho'], 0.01*N, output=grid['rho'], mode='reflect')
    gaussian_filter(grid['etot'], 0.01*N, output=grid['etot'], mode='reflect')
    
    # large internal energy at the center cell for a blast wave
    #grid['etot'][dim*(N/2,)] *= 1000


# set up our parametercluding our callback functions
params = {
   'name': 'scratch',
   'dim': dim,
   'eos_gamma': 1.4,
   'nx': dim*(N,),
   'bctype': 'wall',
   'dx': L/N,
   'init_grid': init_grid,
}

# set up the test
hydroprob = hp.HydroProblem(params_dict=params)
#hydroprob = hp.HydroProblem()

# run in user-specified time chunks. Nice for plotting and dumping output!
for i in xrange(200):
    hydroprob.run(tstop=0.1, max_steps=1, output_every=50)
    hydroprob.plot()

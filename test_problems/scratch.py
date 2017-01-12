#!/home/devon/anaconda2/bin/python

import numpy as np
import sys
from scipy.ndimage.filters import gaussian_filter
sys.path.insert(0, '.')
import hydroPlayground as hp

# global stuff
N = 48 #128 
L = 1.0
dim = 2

# used to initialize the grid cells
def init_bw(self):
    grid = self.getGrid()
    grid['rho'] = 0.1 
    grid['mom'] = 0.0 
    grid['etot'] = 1.0 
    grid['com'] = 0.0

    #nblobs = 200 
    #for sample in np.random.sample((nblobs, self.dim)):
        #hg = tuple(int(s/self.dx) for s in sample)
        #grid['rho'][hg] = 1
    #gaussian_filter(grid['rho'], 0.05/self.dx, output=grid['rho'], mode='reflect')
    #gaussian_filter(grid['etot'], 0.05/self.dx, output=grid['etot'], mode='reflect')
    hg = tuple(n/2 for n in self.nx[:self.dim])
    grid['etot'][hg] *= 100/(self.dx**self.dim) 


# set up our parametercluding our callback functions
params = {
   'name': 'scratch',
   'dim': dim,
   'eos_gamma': 1.4,
   'nx': dim*(N,),
   'bctype': 'wall',
   'dx': L/N,
   'init_grid': init_bw,
}

# set up the test
sod_test = hp.HydroProblem(params_dict=params)
#sod_test = hp.HydroProblem()

# run in user-specified time chunks. Nice for plotting and dumping output!
for i in xrange(20):
    sod_test.run(tstop=0.1, max_steps=5, output_every=50)
    sod_test.plot()

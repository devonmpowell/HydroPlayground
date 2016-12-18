import numpy as np
import hydroPlayground as hp

# plotting function
def make_plots(hp, fsz=6):
    if hp.dim == 1:
        print 'yay'

# set up a callback function
# used to initialize the grid cells
def init_bw(self):
    grid = self.getGrid()
    grid['rho'] = 0.1
    grid['mom'] = 0.0
    grid['etot'] = 0.1

    # this function gives two grid cells a  internal energy
    hg = tuple(n/2 for n in self.nx[:self.dim])
    grid['etot'][hg] = 10 
    hg = tuple(3*n/4 for n in self.nx[:self.dim])
    grid['etot'][hg] = 20 

# set up our parameters
N = 256 
L = 1.0
dim = 2
params = {
   'name': 'blast_waves',
   'dim': dim,
   'eos_gamma': 1.4,
   'cfl_fac': 0.4,
   'nx': dim*(N,),
   'bctype': 'periodic',
   'dx': L/N,
   'init_grid': init_bw,
}

# set up the test
sod_test = hp.HydroProblem(params_dict=params)

# run in user-specified time chunks. Nice for plotting and dumping output!
for i in xrange(10):
    sod_test.run(tstop=0.2, max_steps=1000, output_every=50)
    sod_test.plots()

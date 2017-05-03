#!/home/devon/anaconda2/bin/python

import numpy as np
import matplotlib.pyplot as plt 
import sys
from scipy.ndimage.filters import gaussian_filter
sys.path.insert(0, '.')
import hydroPlayground as hp

# global stuff
N = 64
L = 6.6 # kpc
dim = 3

# used to initialize the grid cells
def init_grid(grid):

    # set the flat background field
    grid['rho'] = 2.94e1 # in 1.0e60 H/kpc
    grid['x'] = 1.0 # ionization fraction 

    # a random set of unit vectors on the sphere
    centers = np.random.normal(size=(200,dim))
    len = np.sum(centers**2, axis=1)**0.5;
    centers /= len[:,None]

    centers *= 1.5

    dx = L/N
    for c in centers:
        #print c
        grid['x'][int(c[0]/dx)+N/2+1,int(c[1]/dx)+N/2+1,int(c[2]/dx)+N/2+1] = 0.0 

    #grid['x'][N/3:2*N/3,2*N/3:5*N/6,N/3:2*N/3] = 0.0
    # add some blobbiness
    #for sample in np.random.randint(0, high=N, size=(100, dim)):
        #grid['x'][sample[0],sample[1],grid['x'].shape[2]/2] = 0.0 
    #gaussian_filter(grid['x'], 0.01*N, output=grid['x'], mode='reflect')
    #grid['x'][grid['x']<0.5] = 0.0
    #grid['x'][grid['x']>0.5] = 1.0

def plot_stromgren(time, radii, fsz=8):


    # From Abel+Wise MORAY
    Rs = 5.4 # Kpc
    Trec = 122.4 # Myr

    fig, ax = plt.subplots(1, 1, figsize=(fsz, fsz))
    ax.set_ylabel(r'$r_s$')
    ax.set_xlabel(r'$t$')
    ax.set_title('$\log{[E_\mathrm{rad}]}$')
    ax.set_xlim(0, 500.0)
    ax.set_ylim(0, 1.1*Rs)

    filt = (time > 0.0)
    t = time[filt]
    ax.plot(t, radii[filt],'k--', label='Numerical')
    t = np.linspace(0, 500.0, 128) 
    ax.plot(t, Rs*(1.-np.exp(-t/Trec))**(1./3.),'r-', label='Analytic')
    plt.show()
    
# default plot-maker
def plots(hypr, fsz=8):
    plgrid = hypr.getGrid()
    etot = plgrid['etot']
    if hypr.dim == 3:
        hypr.fig, ax = plt.subplots(1, 3, figsize=(3*fsz, fsz))
    hypr.fig.subplots_adjust(wspace=0.5)
    if hypr.dim == 3:

        imargs = {'interpolation': 'nearest', 'origin': 'lower', 'extent': [0,1,0,1], 'cmap':
                plt.cm.RdBu_r}

        # check isotropy
        ipts = np.mgrid[:hypr.nx[0],:hypr.nx[1],:hypr.nx[2]].T
        r2 = hypr.dx*hypr.dx*((ipts[:,:,:,0])**2+(ipts[:,:,:,1])**2+(ipts[:,:,:,2])**2)
        iofrac = hypr.pygrid['x'][hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts]
        etot = hypr.pyradgrid['E'][hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts]
        #rhorad = ax[0].imshow(np.log10((etot*r2)[:,:,hypr.nx[2]/2-1]*(4*np.pi*1000)/hypr.dx), **imargs)
        rhorad = ax[0].imshow(np.log10(etot[:,:,hypr.nx[2]/2]), **imargs)
        cax = hypr.fig.add_axes([ax[0].get_position().x1+0.002, 0.21, 0.01, 0.58])
        hypr.fig.colorbar(rhorad, cax=cax)
       
        ax[1].imshow(iofrac[:,:,hypr.nx[2]/2], **imargs)


        #rhorad = ax[0].imshow(np.sum(hypr.pyradgrid['E'][hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts,:],axis=0), **imargs)
        #rhorad = ax[1].imshow(np.sum(hypr.pyradgrid['E'][hypr.nghosts:-hypr.nghosts,:,hypr.nghosts:-hypr.nghosts],axis=1), **imargs)
        #rhorad = ax[2].imshow(np.sum(hypr.pyradgrid['E'][:,hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts],axis=2), **imargs)


        #ax[1].imshow(vel, **imargs)
        #ax[1].tripcolor(triang, hypr.npverts['vel'][:,0], shading='gouraud', cmap=plt.cm.RdBu_r)
        #ax[1].triplot(triang, lw=0.2, c='k', markersize=0, marker=None)
        ax[0].set_ylabel(r'$y$')
        ax[0].set_xlabel(r'$x$')
        ax[0].set_xlim(0, 1)
        ax[0].set_ylim(0, 1)
        ax[0].set_title('$\log{[E_\mathrm{rad}]}$')

        ax[1].set_ylabel(r'$y$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(0, 1)
        ax[1].set_title('Ionized fraction')


        print 'mean ionization fraction =', np.mean(iofrac)
        rstrom = 2*L*(3.0*np.mean(iofrac)/(4*np.pi))**(1./3.)
        print 'rstrom =', rstrom

        ax[2].plot((r2.flatten())**0.5,
                iofrac.flatten(),'r.', label='$E_\mathrm{rad}$', markersize=2)
        #ax[2].semilogy(1.0*(r2.flatten())**0.5,
                #etot.flatten()*r2.flatten()*(4*np.pi*1000)/hypr.dx,'k.', label='$E_\mathrm{rad}$', markersize=2)

        #ax[2].semilogy([0, 0.5], [1, 1],'k--')
        #ax[2].plot([rstrom, rstrom], [-0.1, 1.1],'b--')
        ax[2].set_ylim(-0.1, 1.1)
        ax[2].set_xlim(0.0, L)
        ax[2].set_xlabel(r'$r$')
        ax[2].set_ylabel(r'$x$')
        ax[2].set_title('Ionized fraction $x(r)$')

    plt.show()
    hypr.fig.suptitle('%s, step = %d, t = %.03f' % (hypr.name, hypr.step, hypr.time))

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

rstr = np.zeros(1024)
tstr = np.zeros(1024)

# run in user-specified time chunks. Nice for plotting and dumping output!
for i in xrange(1024):
    hydroprob.run(tstop=10.0, max_steps=6, output_every=50)

    # make a list of the Stromgren radii and the simulation time
    # plot it up
    #if i%50 == 0:
        #plot_stromgren(tstr, rstr)
        ## also pl
    plots(hydroprob)


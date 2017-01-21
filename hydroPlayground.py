import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
import ctypes
from ctypes import *
import warnings
import os 
import time
import inspect 
import itertools
from scipy.stats import norm

# nice plot styling
matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.rcParams['font.size'] = 20.0
matplotlib.rcParams['text.usetex'] = False#True
matplotlib.rcParams['figure.figsize'] = (15,10)
matplotlib.rcParams['lines.linewidth'] = 1.0
matplotlib.rcParams['grid.linewidth'] = 1.5

# load the underlying C shared library
_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
_hp = ctypes.CDLL('%s/hydroplay.so'%_path)

# TODO: add EOS C functions or something 
# default (Sod shock tube in 1D)
# TODO: Figure out where the rho/etot scaling went wrong?? 
def _default_init(self):
    grid = self.getGrid()
    hg = grid.shape[0]/2
    grid['rho'][:hg] = 1.0
    grid['mom'][:hg] = 0.0
    grid['etot'][:hg] = 1.0/(1.0*(self.eos_gamma-1)) 
    grid['rho'][hg:] = 0.125
    grid['mom'][hg:] = 0.0
    grid['etot'][hg:] = 0.01/(0.125*(self.eos_gamma-1)) 
    grid['com'] = 0.0

def _default_output(self):
    print "  - Step %d, t = %f, dt = %.2e" % (self.step, self.time, self.dt)

# default plot-maker
def _default_plots(self, fsz=8):
    plgrid = self.getGrid()
    rho = plgrid['rho']
    mom = plgrid['mom']
    com = plgrid['com']
    if self.dim == 3: # TODO: something more elegant here...
        vel = mom[:,:,:,0]/rho
    if self.dim == 2: # TODO: something more elegant here...
        vel = mom[:,:,0]/rho
        xmean = com/rho[:,:,None]
        xmean[:,:,0] += self.dx*(0.5 + np.arange(self.nx[0]))
        xmean[:,:,1] += self.dx*(0.5 + np.arange(self.nx[1]))[:,None]
    elif self.dim == 1:
        vel = mom[:,0]/rho
        xmean = com[:,0]/rho
        xmean += self.dx*(0.5 + np.arange(self.nx[0]))
    etot = plgrid['etot']
    p = (self.eos_gamma-1.0)*(etot-0.5*rho*vel**2)
    if self.dim == 1:
        self.fig, ax = plt.subplots(4, 1, figsize=(fsz, fsz), sharex=True)
    if self.dim == 2:
        self.fig, ax = plt.subplots(1, 3, figsize=(3*fsz, fsz))
    if self.dim == 3:
        self.fig, ax = plt.subplots(1, 3, figsize=(3*fsz, fsz))
    if self.dim == 1:
        #if hasattr(self, 'fig'):
            #ax = self.fig.get_axes()
        #else:
        x = (0.5+np.arange(self.nx[0]))*self.dx




        ax[0].plot(x, rho) #, drawstyle='steps-mid')
        ax[0].set_ylabel(r'$\rho$')
        ax[0].set_ylim(0.0, 1.2*np.max(rho))

        ax[1].plot(x, vel)
        ax[1].set_ylabel('$v$')
        ax[1].set_ylim(min(-1.0, 1.2*np.min(vel)), max(1.0, 1.2*np.max(vel)))

        ax[2].plot(x, p)
        ax[2].set_ylabel('$p$')
        ax[2].set_ylim(0.0, 1.2*np.max(p))

        dcm = (xmean[1:]-xmean[:-1])/self.dx
        ax[3].plot(x, xmean/self.dx)
        ax[3].set_ylabel('$\Delta_\mathrm{cm}$')

        ax[0].plot(xmean, rho)
        ax[1].plot(xmean, vel)
        ax[2].plot(xmean, p)


    if self.dim == 2:
        imargs = {'interpolation': 'nearest', 'origin': 'lower', 'extent': [0,1,0,1], 'cmap':
                plt.cm.RdBu_r}

        tetinds = self.nptets['verts'][:,:3]
        triang = mtri.Triangulation(self.npverts['pos'][:,0], self.npverts['pos'][:,1], tetinds)

        #ax[0].imshow(np.log10(radsymm), **imargs)
        #print np.min(1.0-radsymm[radsymm > 0]), np.max(1.0-radsymm)
        ax[0].imshow(rho, **imargs)
        #ax[0].tripcolor(triang, self.npverts['rho'], shading='gouraud', cmap=plt.cm.RdBu_r)
        #ax[0].triplot(triang, lw=0.2, c='k', markersize=0, marker=None)
        ax[0].set_ylabel(r'$y$')
        ax[0].set_xlabel(r'$x$')
        ax[0].set_xlim(0, 1)
        ax[0].set_ylim(0, 1)
        ax[0].set_title('Matter density')

        radsymm = 0.25*(self.pyradgrid['E'] + np.roll(self.pyradgrid['E'][::-1,:], 1, axis=0) + np.roll(self.pyradgrid['E'][:,::-1], 1, axis=1) +
            np.roll(np.roll(self.pyradgrid['E'][::-1,::-1], 1, axis=0), 1, axis=1))
        #ax[1].imshow(np.log10(0.0000001+radsymm[self.nghosts:-self.nghosts,self.nghosts:-self.nghosts]), **imargs)
        ax[1].imshow(np.log10(self.pyradgrid['E'][self.nghosts:-self.nghosts,self.nghosts:-self.nghosts]), **imargs)
        #ax[1].imshow(vel, **imargs)
        #ax[1].tripcolor(triang, self.npverts['vel'][:,0], shading='gouraud', cmap=plt.cm.RdBu_r)
        #ax[1].triplot(triang, lw=0.2, c='k', markersize=0, marker=None)
        ax[1].set_ylabel(r'$y$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_xlim(0, 1)
        ax[1].set_ylim(0, 1)
        #ax[1].set_title('velocity (x)')
        ax[1].set_title('radiation flux')

        ax[2].imshow(p, **imargs)
        #ax[2].tripcolor(triang, self.npverts['etot'], shading='gouraud', cmap=plt.cm.RdBu_r)
        #ax[2].tripcolor(triang, np.sum()self.npverts['rho'][self.nptets['verts']], shading='gouraud', cmap=plt.cm.RdBu_r)
        #ax[2].triplot(triang, lw=0.2, c='k', markersize=0, marker=None)
        ax[2].set_ylabel(r'$y$')
        ax[2].set_xlabel(r'$x$')
        ax[2].set_title('pressure')
        ax[2].set_xlim(0, 1)
        ax[2].set_ylim(0, 1)
        #ax[2].set_title('etot')
        #ax[2].set_title('mass error')
        #print "rhomin =", np.min(self.npverts['rho'])
        #print "rhomax =", np.max(self.npverts['rho'])

        #ax[3].scatter(xmean[:,:,0], xmean[:,:,1], c=rho, s = 4, lw = 0)
        #ax[3].set_aspect('equal')
        #ax[3].set_xlim(0,1)
        #ax[3].set_ylim(0,1)
        #ax[3].set_ylabel(r'$y$')
        #ax[3].set_xlabel(r'$x$')
        #ax[3].set_title('cell center of mass')


        #dcm = np.sqrt(np.sum((xmean[1:,1:,:]-xmean[:-1,1:,:])**2 +
            #(xmean[1:,1:,:]-xmean[1:,:-1,:])**2, axis=-1))/self.dx
        #ax[3].imshow(dcm.T, **imargs)
        #ax[3].set_ylabel(r'$y$')
        #ax[3].set_xlabel(r'$x$')
        #ax[3].set_title(r'$\Delta_\mathrm{cm}$')

        #print "Plotting radiation field!"
        from matplotlib.patches import Wedge
        from matplotlib.collections import PatchCollection
        #print self.pyrad

        allrays = [Wedge(self.dx*(0.5+ray['orcell'][:2]), ray['rmax'], ray['angle_id']*360./128, \
                (ray['angle_id']+1)*360./128, width=(ray['rmax']-ray['rmin']), alpha = ray['N']) for ray in self.pyrad[:self.nrays]]
        #ax[1].add_collection(PatchCollection(allrays, lw = 1, facecolor='white', alpha = 0.09))

        # verticies by subtracting random offsets from those center-points
        #numpoly, numverts = 100, 4
        #centers = np.random.random((numpoly,2))
        #offsets = 0.1 * np.random.random((numverts,numpoly,2))
        #verts = centers + offsets
        #verts = np.swapaxes(verts, 0, 1)
        
        ## Color scalar...
        ## in as "facecolors=colorval" when you create the PolyCollection
        #z = np.random.random(numpoly) * 500
        #coll = PolyCollection(verts, array=z, cmap=plt.cm.jet, edgecolors='none')
        #ax[3].add_collection(coll)
        #ax[3].set_title('radiation')

    if self.dim == 3:
        imargs = {'interpolation': 'nearest', 'origin': 'lower'}

        #if hasattr(self, 'fig'):
            #ax = self.fig.get_axes()
        #else:
        plt.show()

        ax[0].imshow(rho[:,:,32], **imargs)
        ax[0].set_ylabel(r'$y$')
        ax[0].set_xlabel(r'$x$')
        ax[0].set_title('rho')

        ax[1].imshow(vel[:,:,32], **imargs)
        ax[1].set_ylabel(r'$y$')
        ax[1].set_xlabel(r'$x$')
        ax[1].set_title('vel')

        ax[2].imshow(p[:,:,32], **imargs)
        ax[2].set_ylabel(r'$y$')
        ax[2].set_xlabel(r'$x$')
        ax[2].set_title('pressure')

    #print 'delta cm [ %f , %f]' % (np.min(dcm), np.max(dcm))
    plt.show()
    ax[-1].set_xlabel(r'$x$')
    self.fig.suptitle('%s, step = %d, t = %.03f' % (self.name, self.step, self.time))

# parameters for a 1D problem with 256 grid cells
# this defines every parameter that HydroPlayground recognizes
_default_params = {
        'name': 'default_hp_problem',
        'dim': 1,
        'eos_gamma': 7./5,
        'cfl_fac': 0.4,
        'nx': (256,),
        'bctype': 'wall',
        'dx': 1.0/256,
        'init_grid': _default_init, 
        'on_output': _default_output, 
}

class HydroProblem(Structure):

    # c struct fields 
    _fields_ = [("name", c_char_p), ("time", c_double), ("dt", c_double), ('step', c_int),
        ("eos_gamma", c_double), ("eos_p", c_void_p), ("cfl_fac", c_double),
        ("dx", c_double), ("bctype", c_int), ("grid", c_void_p), 
        ("dim", c_int), ("nx", 3*c_int), ("strides", 3*c_int),
        ("on_output", CFUNCTYPE(c_void_p, c_void_p)),
        ("rays", c_void_p), ("nrays", c_int), ("base_ray_lvl", c_int), ("max_ray_lvl", c_int), 
        ("rad_grid", c_void_p),
        ("mesh_verts", c_void_p), ("mesh_tets", c_void_p), ("nverts", c_int), ("ntets", c_int)] 

    def __init__(self, params_dict=_default_params, verbose=True):

        self.verbose = verbose
        if self.verbose:
            print "\n------------------------------"
            print "                              "
            print " Welcome to HydroPlayground..."
            print "   have fun!"
            print "                              "
            print "------------------------------"

        # load all of the supplied parameter files
        if self.verbose:
            print "\nLoading parameters..."
        cfields = [f[0] for f in self._fields_]
        bctypes = {'wall': 0, 'periodic': 1, 'free': 2}
        for dp, val in _default_params.iteritems():
            if dp in params_dict:
               val = params_dict[dp]
            if dp == 'bctype':
                self.bcstr = dp
                setattr(self, dp, bctypes[val])
            elif dp == 'on_output':
                output_fn = val
                setattr(self, dp, CFUNCTYPE(c_void_p, c_void_p)(lambda x: output_fn(self)))
            else:
                setattr(self, dp, val)
            if self.verbose:
                print ' -', dp, '=', val 

        # allocate and initialize the grid
        # setup the problem by calling the grid setup callback
        if self.verbose:
            print "\nSetting up the grid..."
        self.nghosts = 2 
        self.dtype = np.dtype([('rho', np.float64), ('mom', np.float64, (3,)), ('com', np.float64, (3,)), ('etot', np.float64)]);
        self.pygrid = np.ones(tuple(axsz+2*self.nghosts for axsz in self.nx[:self.dim]), dtype=self.dtype)
        tstr = tuple(np.cumprod(np.append(1,self.pygrid.shape))[:self.dim])
        self.strides = tstr

        # allocate and initialize the grid
        # setup the problem by calling the grid setup callback
        if self.verbose:
            print "\nSetting up the radiation field buffer..."
        self.rdtype = np.dtype([('angle_id', np.int64), ('rmin', np.float64), ('rmax', np.float64), \
            ('orcell', np.int32, (4,)), ('N', np.float64),]) 
        self.pyrad = np.zeros(16000, dtype=self.rdtype)
        self.nrays = 0

        # Testing the radiation energy density
        self.pyradgrid = np.ones(tuple(axsz+2*self.nghosts for axsz in self.nx[:self.dim]), dtype=np.dtype([('E', np.float64),]))

        if self.verbose:
            print " - underlying numpy.ndarray shape =", self.pygrid.shape
            print " - radiation field buffer shape =", self.pyrad.shape
            print " - strides =", tstr 
            print " - nghosts =", self.nghosts
            print " - calling %s... " % self.init_grid.__name__,
        self.init_grid(self.pygrid)
        self.grid = self.pygrid.ctypes.data_as(c_void_p) 
        self.rays = self.pyrad.ctypes.data_as(c_void_p) 
        self.rad_grid = self.pyradgrid.ctypes.data_as(c_void_p) 
        if self.verbose:
            print " Success!"

        # kick it off!
        self.step = 0
        self.wallstart = time.time()

        # EXPERIMENTAL: initialize the Lagrangian tet mesh
        # TODO: 2D only for now...
        # create a Delaunay triangulation
        # with some boundary guards
        ldim = 48 
        vpos = (1.0/ldim*np.mgrid[:ldim+1,:ldim+1]).T
        vpos[1:-1,1:-1] += 0.6/ldim*(0.5-np.random.sample((ldim-1, ldim-1, 2)))
        vpos = vpos.reshape((ldim+1)**2, 2)
        self.nverts = vpos.shape[0]
        triang = mtri.Triangulation(vpos[:,0], vpos[:,1])
        self.ntets = triang.triangles.shape[0]

        # initialize the vertices
        vtype = np.dtype([('pos', np.float64, (3,)), ('vel', np.float64, (3,)), ('rho', np.float64),  ('etot', np.float64)]);
        self.npverts = np.zeros((self.nverts), dtype=vtype)
        self.npverts['rho'] = 1.0
        self.npverts['pos'][:,:2] = vpos
        self.npverts['vel'] = 0.0
        self.npverts['etot'] = 1.0
        self.mesh_verts = self.npverts.ctypes.data_as(c_void_p) 

        len = np.sqrt(np.sum((vpos-0.5)**2, axis=1))
        norm = (vpos-0.5)
        #norm[:,0] /= len
        #norm[:,1] /= len
        self.npverts['vel'][:,0] = 10.0*norm[:,0]*np.exp(-((vpos[:,0]-0.5)**2+(vpos[:,1]-0.5)**2)/(2*0.1**2)) 
        #self.npverts['vel'][:,1] = 100.0*norm[:,1]*np.exp(-((vpos[:,0]-0.5)**2+(vpos[:,1]-0.5)**2)/(2*0.1**2)) 

        # initialize the tets
        ttype = np.dtype([('rho', np.float64), ('mom', np.float64, (3,)), ('com', np.float64, (3,)),
            ('etot', np.float64), ('verts', np.int32, (4,))])
        self.nptets = np.zeros(self.ntets, dtype=ttype)
        self.nptets['verts'][:,:3] = triang.triangles
        self.mesh_tets = self.nptets.ctypes.data_as(c_void_p) 

        # do C things like find grid strides, etc
        _hp.init(byref(self))


    def __del__(self):
        if self.verbose:
            print "\n------------------------------"
            print "                              "
            print " Goodbye!                     "
            print "                              "
            print "------------------------------"

    def run(self, tstop=0.2, max_steps=1000, output_every=10, output_callback=None):
        if self.verbose:
            print "\nRunning..."
            print " - tstop =", tstop
            print " - current steps = %d, max = %d" % (self.step, max_steps) 
            print " - callback every %d steps" % output_every 
        tstart = time.time()
        _hp.evolve(c_double(tstop), max_steps, output_every, byref(self))
        tstop = time.time()
        if self.verbose:
            print " - wall time = %.03f s (total elapsed = %.03f s)" % (tstop-tstart, tstop-self.wallstart)

    def plot(self):
        # TODO: clean this up...
        _default_plots(self)
    
    def getGrid(self):
        if self.dim == 3:
            return self.pygrid[self.nghosts:-self.nghosts,self.nghosts:-self.nghosts,self.nghosts:-self.nghosts]
        if self.dim == 2:
            return self.pygrid[self.nghosts:-self.nghosts,self.nghosts:-self.nghosts]
        if self.dim == 1:
            return self.pygrid[self.nghosts:-self.nghosts]

    def getRays(self):
        return self.pyrad[:self.nrays]

    def write(self):
        # TODO: figure out how to do output callbacks!!!
        # TODO: do the output dumps all in Python??
    	_psi.write(byref(self))



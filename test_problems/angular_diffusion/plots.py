#!/home/devon/anaconda2/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import sys
sys.path.insert(0, '.')
import hydroPlayground as hp

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'computer modern roman'

# global stuff
N = 64
L = 6.6 # kpc
dim = 3

def plot_stromgren(time, radii, fsz=8):

    # From Abel+Wise MORAY
    Rs = 5.4 # Kpc
    Trec = 122.4 # Myr

    fig, ax = plt.subplots(2, 1, figsize=(8, 5), gridspec_kw={'height_ratios': [2.7,1.3]},sharex=True)
    fig.subplots_adjust(bottom = 0.15, hspace=0.05)
    ax[0].set_ylabel(r'$r_s$ (kpc)')
    ax[1].set_xlabel(r'$t$ (Myr)')
    #ax[0].set_title('$R_s = 5.4$ kpc, $t_\mathrm{rec} = 122$ Myr')
    ax[0].set_xlim(0, 500.0)
    ax[0].set_ylim(0, 1.2*Rs)

    filt = (time > 0.0)
    t = time[filt]
    #t = np.linspace(0, 500.0, 128) 
    analytic = Rs*(1.-np.exp(-t/Trec))**(1./3.)
    ax[0].plot(t, analytic,'k-', 
            label='$r_s(t)=R_s[1-\exp(-t/t_\mathrm{rec})]^{1/3}$')
    ax[0].plot(t, radii[filt],'k--', label='This work')
    ax[0].plot(t, Rs*np.ones_like(t), 'k.', markersize=0.1) 
    ax[0].text(0.2*Trec, Rs*1.01, '$R_s$')
    ax[0].legend(loc='lower right', fontsize=18)

    ax[1].set_ylabel(r'\% Err.')
    ax[1].set_ylim(-9,9)
    ax[1].plot(t, np.zeros_like(t), 'k.', markersize=0.1) 
    ax[1].plot(t, (radii[filt]/analytic-1.0)*100,'k-') 

    plt.show()
    fig.savefig('test_problems/stromgren/stromgren_time.png')
    
# default plot-maker
def plots(hypr, fsz=8):
    plgrid = hypr.getGrid()
    etot = plgrid['etot']
    if hypr.dim == 3:
        hypr.fig, ax = plt.subplots(1, 3, figsize=(3*fsz, fsz))
    hypr.fig.subplots_adjust(wspace=0.5)
    if hypr.dim == 3:

        imargs = {'interpolation': 'nearest', 'origin': 'lower', 'extent': [0,1,0,1], 'cmap':
                plt.cm.RdBu_r, 'vmin': -8}

        # check isotropy
        ipts = np.mgrid[:hypr.nx[0],:hypr.nx[1],:hypr.nx[2]].T
        r2 = hypr.dx*hypr.dx*((ipts[:,:,:,0])**2+(ipts[:,:,:,1])**2+(ipts[:,:,:,2])**2)
        iofrac = hypr.pygrid['x'][hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts]
        etot = hypr.pyradgrid['E'][hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts,hypr.nghosts:-hypr.nghosts]
        #rhorad = ax[0].imshow(np.log10((etot*r2)[:,:,hypr.nx[2]/2-1]*(4*np.pi*1000)/hypr.dx), **imargs)
        rhorad = ax[0].imshow(np.log10(etot[:,:,0]), **imargs)
        cax = hypr.fig.add_axes([ax[0].get_position().x1+0.002, 0.21, 0.01, 0.58])
        hypr.fig.colorbar(rhorad, cax=cax)
       
        ax[1].imshow(iofrac[:,:,0], **imargs)


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

i = 200
rstr = np.fromfile('test_problems/stromgren/rstr_%04d.np'%i)
tstr = np.fromfile('test_problems/stromgren/time_%04d.np'%i)
plot_stromgren(tstr, rstr)

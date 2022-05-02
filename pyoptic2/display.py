import matplotlib.pyplot as plt
import numpy as np


def plot2d(rays, sli='xz'):
    

    def plot_ray(ray):
        ra = np.concatenate([np.atleast_3d(ri.p0) for ri in ray], 2)
        #ra = np.array([ri.p0 for ri in ray])
        if not ray[-1].p1 == None:
            ra = np.concatenate([ra, np.atleast_3d(ray[-1].p1)], 2)
        
        if sli == 'xz':
            plt.plot(ra[:, 2,:].squeeze().T, ra[:, 0,:].squeeze().T, c=ray[0].color)
        elif sli == 'zx':
            plt.plot(ra[:, 0, :].squeeze().T, ra[:, 2, :].squeeze().T, c=ray[0].color)
        elif sli == 'xy':
            plt.plot(ra[:, 1, :].squeeze().T, ra[:, 0, :].squeeze().T, c=ray[0].color)
        else:
            # sli is a tuple of (right, up) vectors 
            assert len(sli)==2

            right = np.array(sli[0])
            up = np.array(sli[1])

            x = (ra[:,:,:]*right[None,:,None]).sum(1).squeeze().T
            y = (ra[:,:,:]*up[None,:,None]).sum(1).squeeze().T

            plt.plot(x, y, c=ray[0].color)
    
    
    for ray in rays:
        plot_ray(ray)


def plot_system2d(sys, sli='xz'):
    if sli == 'xz':
        def plot_surf(s):
            sf = s.surface(proj='y')
            plt.plot(sf[2], sf[0], 'k')
    
    elif sli == 'zx':
        def plot_surf(s):
            sf = s.surface(proj='y')
            plt.plot(sf[0], sf[2], 'k')
    
    elif sli == 'xy':
        def plot_surf(s):
            sf = s.surface(proj='z')
            #print('sf.shape:', sf.shape)
            plt.plot(sf[1], sf[0], 'k')
            
    elif sli == 'yx':
        def plot_surf(s):
            sf = s.surface(proj='z')
            plt.plot(sf[0], sf[1], 'k')

    else:
        # slice is a tuple (right, up)
        assert len(sli) == 2
        def plot_surf(s):
            right = np.array(sli[0])
            up = np.array(sli[1])

            sf = s.surface(proj=(right, up))

            x = (sf*right[:,None,None]).sum(0)
            y = (sf*up[:,None,None]).sum(0)
            plt.plot(x, y, 'k')

    
    for s in sys:
        plot_surf(s)


def dotplt(rays):
    NA = max([np.sin(np.arccos(np.dot(rays[-1].d[0], d_i))) for d_i in rays[-1].d[1:]])
    #print(rays[-1].d[0])
    #print([np.sin(np.arccos(np.dot(rays[-1].d[0], d_i))) for d_i in rays[-1].d[1:]])
    #print(rays[-1])
    #print(NA)
    xh = np.cross(rays[-1].d[0], rays[-1].d[1])
    yh = np.cross(xh, rays[-1].d[0])
    
    xh = xh/np.linalg.norm(xh)
    yh = yh/np.linalg.norm(yh)
    
    r_dl = rays[-1].wavelength * 1e-6 / (2 * NA)
    pts = rays[-1].p0 #np.array([r_i[-1].p0 for r_i in rays])
    
    x1 = (pts*xh[None,:]).sum(1)
    y1 = (pts * yh[None, :]).sum(1)
    plt.plot(x1 -x1.mean(), y1 - y1.mean(), '.')
    
    t = np.linspace(0, 2 * np.pi)
    plt.plot(r_dl * np.sin(t), r_dl * np.cos(t))
    plt.axis('equal')
    #xlim(-.06, .06)
    plt.ylim(-r_dl*3, r_dl*3)
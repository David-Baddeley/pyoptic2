import numpy as np

class RayBundle :
    ''' Class describing a light ray from an initial point to a final point '''

    def __init__(self, point, dir, material, wavelength=635., color=(1.0, 0, 0), cumulativePath=0, intensities=None, stokes=None) :
        '''Construct a ray from a start point and direction dir'''
        #print "Ray:__init__>"
        
        # initial points of the ray
        self.p0 = np.array(point)
        self.d  = np.array(dir)
        self.d  = self.d/np.linalg.norm(self.d, axis=-1, keepdims=True)
        
        self.material = material
        self.wavelength = wavelength

        # values after the ray has been propagated (empty to start with)
        self.p1 = None
        if (intensities is None):
            if (self.d.ndim == 1):
                self.intensities = np.array([1,])
            else:
                self.intensities = np.ones(self.d.shape[0], 'f')
        else:
            self.intensities = intensities

        self.stokes = stokes
        self.color = color
        
        self.prev_pathlength = cumulativePath

    def propagate(self, lam) :
        #print 'lam = ', lam
        if lam is None:
            return None
        if isinstance(lam, np.ndarray):
            return lam[:,None]*self.d + self.p0
        else:
            return lam*self.d + self.p0
        
    def prop(self, lam):
        """ Like propagate, but returns a Placement object"""
        
        from . import placement
        return placement.Placement(np.squeeze(self.propagate(lam)), np.squeeze(self.d))

    def __str__(self) :
        s  = ''
        s += 'Ray                      : \n'
        s += 'Ray p0                   : '+str(self.p0)+'\n'
        s += 'Ray p1                   : '+str(self.p1)+'\n'
        s += 'Ray d                    : '+str(self.d)
        return s
        
    @property
    def length(self):
        return np.linalg.norm(self.p1-self.p0, axis=-1)
    
    @property
    def pathLength(self):
        return self.length*self.material.n(self.wavelength)
    
    @property
    def cumulativePath(self):
        return self.prev_pathlength + self.pathLength
    
    @property
    def focus(self):
        return self._focus()

    def _focus(self, apodisation=None):
        # cumulative path based ... 
        d_ = self.d[0]
        
        #print('d_:', d_)
        #print (self.p0.shape, self.d.shape)
        
        p_ = self.p0[1:, :] - self.p0[0, :][None, :]
        
        #print p_.shape
        # axial component of position
        p_axial = (p_*d_[None,:]).sum(axis=-1)

        pl = self.prev_pathlength[1:] - self.prev_pathlength[0] #- p_axial


        #radial component of distance from centre ray
        p_radial = p_ - p_axial[:,None]*d_[None,:]

        p2 = np.linalg.norm(p_radial, axis=-1)
        pd_ = p_radial/p2[:,None]
        
        #how fast are we propagating in the radial direction?
        d_radial = (pd_*self.d[1:,:]).sum(-1)

        theta = np.arcsin(d_radial)

        if apodisation:
            weights=apodisation(theta)
        else:
            weights = np.ones_like(theta)
        
        #print(p2, d_radial)
        #l = np.linalg.lstsq(d_radial.reshape(-1, 1), p2)

        #print(p_axial)
        #print((p2/d_radial))

        #print(pl, self.d)

        return (weights*(p_axial + pl/(p_axial + 1-1.0/(self.d[1:,:]*d_[None, :]).sum(-1)))).sum()/weights.sum()

        #return np.mean(p2/d_radial)
        
        #l = np.linalg.lstsq(d_radial.reshape(-1, 1), p2)
        
        #return l[0]
        
        #r =  (pd_*self.d[1:, :]).sum(axis=-1)
    
    def __focus(self, apodisation=None):
        d_ = self.d[0]
        
        #print self.p0.shape, self.d.shape
        
        p_ = self.p0[1:, :] - self.p0[0, :][None, :]
        
        #print p_.shape
        # axial component of position
        p_axial = (p_*d_[None,:]).sum(axis=-1)
        #radial component of distance from centre ray
        p_radial = p_ - p_axial[:,None]*d_[None,:]

        p2 = np.linalg.norm(p_radial, axis=-1)
        pd_ = p_radial/p2[:,None]
        
        #how fast are we propagating in the radial direction?
        d_radial = (pd_*self.d[1:,:]).sum(-1)

        theta = np.arcsin(d_radial)

        if apodisation:
            weights=apodisation(theta)
        else:
            weights = np.ones_like(theta)
        
        #print(p2, d_radial)
        #l = np.linalg.lstsq(d_radial.reshape(-1, 1), p2)

        #print(p_axial)
        #print((p2/d_radial))

        return (weights*(p2/d_radial - p_axial)).sum()/weights.sum()

        #return np.mean(p2/d_radial)
        
        #l = np.linalg.lstsq(d_radial.reshape(-1, 1), p2)
        
        #return l[0]
        
        #r =  (pd_*self.d[1:, :]).sum(axis=-1)
        
        

def RayExample() :
    return RayBundle([0, 0, 0], [0, 0, 1])

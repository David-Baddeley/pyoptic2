import numpy as np

class RayBundle :
    ''' Class describing a bundle of light rays from initial points to final points '''

    def __init__(self, point, dir, material, wavelength=635., color=(1.0, 0, 0), cumulativePath=0, intensities=None, stokes=None) :
        '''Construct a raybundle from a start points and direction dir
        
        Parameters
        ----------
        point: array-like, shape (N, 3)
            Starting point(s) of the ray(s)
        dir: array-like, shape (N, 3)
            Direction(s) of the ray(s)
        material: Material
            Material through which the ray(s) propagate
        wavelength: float, optional
            Wavelength of the ray(s) in nm (default is 635)
        color: tuple, optional
            RGB color of the ray(s) (default is red) (display only, does not affect ray propagation)
        cumulativePath: float, optional
            Initial cumulative path length (default is 0)
        intensities: array-like, optional
            Intensities of the ray(s) (default is 1)
        stokes: array-like, optional
            Stokes parameters of the ray(s) (default is None)
        '''
        #print "Ray:__init__>"
        
        # initial points of the ray
        self.p0 = np.array(point)
        if self.p0.ndim == 1:
            self.p0 = self.p0[None,:]
        self.d  = np.array(dir)
        if self.d.ndim == 1:
            self.d = self.d[None,:]
            
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
    
    def copy(self):
        """Create a deep copy of the RayBundle
        
        Returns
        -------
        RayBundle
            A new RayBundle with copied arrays and attributes
        """
        rb = RayBundle(
            point=self.p0.copy(),
            dir=self.d.copy(),
            material=self.material,  # Materials are typically immutable, so no need to copy
            wavelength=self.wavelength,
            color=self.color,
            cumulativePath=self.prev_pathlength.copy() if isinstance(self.prev_pathlength, np.ndarray) else self.prev_pathlength,
            intensities=self.intensities.copy() if self.intensities is not None else None,
            stokes=self.stokes.copy() if self.stokes is not None else None
        )
        
        # Copy the propagated endpoint if it exists
        if self.p1 is not None:
            rb.p1 = self.p1.copy()
        
        return rb
        
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
        """ Find the best compromise focus for the ray bundle, by finding the average distance 
        each ray would need to propagate to reach a plane where the path length difference to 
        the principle ray is zero.
        
        Parameters
        ----------

        apodisation: function (ndarray angles) -> ndarray weights
            Function of angle to principle ray, to weight the rays in the focus calculation.
            This is important for high NA systems where both the transmission losses (apodisation) 
            at high angles can be significant, the abberrations / focal error of the extreme rays 
            can be large, and the number of rays at these angles can also be large (if using constant
            angular spacing, which is the default). If we didn't account for this the high angle
            rays would dominate the calculation. 
        """
        
        # principle ray direction
        d_ = self.d[0]
        
        # relative postion of start of each ray to principle ray
        p_ = self.p0[1:, :] - self.p0[0, :][None, :]
    
        # axial component of position
        p_axial = (p_*d_[None,:]).sum(axis=-1)

        # accumulated path at the start of the ray, relative to the principle ray
        pl = self.prev_pathlength[1:] - self.prev_pathlength[0] #- p_axial

        n_current = self.material.n(self.wavelength)

        if apodisation:
            # Find the radial component to estimate the angle of the ray to the priciple ray (for apodisation if required)
            # NB - radial component is not used for the focus calculation in this version.
            # radial component of distance from centre ray
            p_radial = p_ - p_axial[:,None]*d_[None,:]

            p2 = np.linalg.norm(p_radial, axis=-1)
            pd_ = p_radial/p2[:,None]
            
            #how fast are we propagating in the radial direction?
            d_radial = (pd_*self.d[1:,:]).sum(-1)

            theta = np.arcsin(d_radial)
        
            weights=apodisation(theta)
        else:
            weights = np.ones_like(theta)

       
        # Find the (weighted) average distance each ray would need to propagate to that the difference
        # in path length to the principle ray is zero. 
        return (weights*(p_axial + pl/(p_axial + 1-1.0/(self.d[1:,:]*d_[None, :]).sum(-1)))).sum()/n_current*weights.sum()


    
    def __focus(self, apodisation=None):
        """ Find the best compromise focus for the ray bundle, on a geometric basis by finding the 
        average distance each ray would need to propagate to to have a radial distance of zero to the
        principle ray. 
        
        If we have a full pathlength history, it is likely better to use the pathlength based focus
        calculation above (especially for high NA systems).
        """
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

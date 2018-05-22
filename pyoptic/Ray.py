import numpy as np

class RayBundle :
    ''' Class describing a light ray from an initial point to a final point '''

    def __init__(self, point, dir, material, wavelength=635., color=(1.0, 0, 0), cumulativePath=0) :
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
        self.pha = 0 
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

def RayExample() :
    return RayBundle([0, 0, 0], [0, 0, 1])

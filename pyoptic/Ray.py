import pylab as pl

class Ray :
    ''' Class describing a light ray from an initial point to a final point '''

    def __init__(self, point, dir, material, wavelength=635., color=(1.0, 0, 0)) :
        '''Construct a ray from a start point and direction dir'''
        #print "Ray:__init__>"
        
        # initial points of the ray
        self.p0 = pl.array(point)
        self.d  = pl.array(dir)
        self.d  = self.d/pl.linalg.norm(self.d)
        
        self.material = material
        self.wavelength = wavelength

        # values after the ray has been propagated (empty to start with)
        self.p1 = None
        self.pha = 0 
        self.color = color

    def propagate(self, lam) :
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
        return pl.norm(self.p1-self.p0)
    
    @property
    def pathLength(self):
        return self.length*self.material.n(self.wavelength)

def RayExample() :
    return Ray([0,0,0],[0,0,1])

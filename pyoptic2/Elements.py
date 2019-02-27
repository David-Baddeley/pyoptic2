import numpy as np

import material
from rays       import RayBundle

SHAPE_RECT, SHAPE_CIRC = range(2)


class Element(object) :
    ''' Orientable volume, base class to describe extent, material and orientation of a volume containing '''
            
    def __init__(self, placement, material, name='', shape=SHAPE_CIRC, dimension = [12.5, 12.5, 12.5], **kwargs) :
        self.name      = name
        self.shape     = shape
        self.dimension = np.array(dimension)
        self.placement = placement
        self.material  = material

    def  __str__(self) :
        s  = 'Volume\n'
        s += 'Volume.shape             : '+str(self.shape)+'\n'
        s += 'Volume.dimension         : '+str(self.dimension)+'\n'
        s += 'Volume.placement         : \n'+str(self.placement)
        s += 'Volume.material          : \n'+str(self.material)
        return s
    
    @property
    def orientation(self):
        return self.placement.orientation
    
    @property
    def location(self):
        return self.placement.location
    
    def _surface_coordinates(self, proj=None):
        if False:#proj in  ['x', 'y', 'z']:
            xx = np.linspace(-self.dimension[0], self.dimension[0])[:,None]
            yy = 0 * xx
        #elif proj == 'y':
        #    yy = np.linspace(-self.dimension[0], self.dimension[0])
        #    xx = 0 * yy
        elif self.shape == SHAPE_RECT:
            if not proj is None:
                xd, yd, _ = self.dimension
                xx = np.array([-xd, xd, xd, -xd, -xd])[:,None]
                yy = np.array([-yd, -yd, yd, yd, -yd])[:,None]
            else:
                x = np.arange(-self.dimension[0], self.dimension[0] + 1e-8, self.dimension[0] / 5)
                y = np.arange(-self.dimension[1], self.dimension[1] + 1e-8, self.dimension[1] / 5)
                xx, yy = np.meshgrid(x, y)
        else:
            if not proj is None:
                r = self.dimension[0] * np.sqrt(np.arange(0, 1.01, .1))
            else:
                r = self.dimension[0] * np.sqrt(np.arange(0, 1.01, .1))
            t = 2 * np.pi * np.arange(0, 1.01, .1)
        
            rr, tt = np.meshgrid(r, t)
            pol = rr * np.exp(1j * tt)
            xx = pol.real
            yy = pol.imag
            
        return xx, yy
    
    def _orientate_surface(self, xx, yy, zz, proj=None):
        #if np.argmax(np.abs(self.placement.orientation)) == 2:
            #we are mostly orientated along z
            #if proj == 'x':
            #    N = np.array([1,0,0])
            #else:
            #    N = np.array([1., 0, 0])
            #print 'oriented along z'
        #else:
        #    N = np.array([0, 0, 1.])
        
        if proj == 'z':
            N = np.array([0,0,1])
        elif proj == 'y':
            N = np.array([0,1,0])
        elif proj == 'x':
            N = np.array([1,0,0])
        else:
            N = np.array([0,0,1])
            
            if np.abs(np.dot(N, self.placement.orientation)) >.9:
                #object is aligned with z axis, use x axis instead
                N = np.array([1,0,0])
            
            
        ax = np.cross(N, self.placement.orientation)
        ax = ax / np.linalg.norm(ax)
        
        N2 = np.cross(self.placement.orientation, ax)
        N2 = N2/np.linalg.norm(N2)
        
        coords = self.placement.orientation[:,None,None]*zz[None,:, :] + ax[:,None,None]*xx[None,:,:] + N2[:,None, None]*yy[None,:,:]
        return coords + self.placement.location[:,None, None]
    
       

############################################################################
# General optical surface
############################################################################            
class OpticalSurface(Element) :
    def __init__(self) :
        pass

    def propagate(self,inray) :
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        # compute out going ray
        
        if self.material.type == material.Material.REFLECT :
            outray = reflect(inray,sn)
        elif self.material.type == material.Material.REFRACT :
            outray = snell(inray,sn,inray.material,self.material)

        return outray
        
    def intersection(self,inray) :
        pass
    def surfaceNormal(self,inray) :
        pass

############################################################################
# Plane surface
############################################################################            
class PlaneSurface(OpticalSurface) :
    def __init__(self,*args, **kwargs) :
        Element.__init__(self, *args, **kwargs)

    def surface(self, proj=None) :
        xx, yy = self._surface_coordinates(proj)
        zz = np.zeros_like(xx)
        
        return self._orientate_surface(xx, yy, zz, proj)
        

    def intersection(self,inray) :
        lam = abs(((self.placement.location-inray.p0)*(self.placement.orientation)).sum(-1)/(inray.d* self.placement.orientation).sum(-1))
        inray.p1 = inray.propagate(lam)
                
    def surfaceNormal(self,inray) :
        return self.placement.orientation
    
    
    def __str__(self) :
        s  = 'PlaneSurface             : '+self.name+'\n'
        s += 'PlaneSurface.Volume      :\n' + Element.__str__(self)
        return s

############################################################################
# SphericalSurface
############################################################################            
class SphericalSurface(OpticalSurface) :
    def __init__(self, *args, **kwargs) :
        Element.__init__(self, *args, **kwargs)
        self.radcurv   = kwargs['curvature_radius']

    def surface(self, proj=None) :
        xx, yy = self._surface_coordinates(proj)
        zz = np.sign(self.radcurv)*(-np.sqrt(self.radcurv**2-xx**2-yy**2) +np.abs(self.radcurv))
        
        return self._orientate_surface(xx, yy, zz, proj)

    def intersection(self,ray) :
        dv = ray.p0 - self.placement.orientation*self.radcurv - self.placement.location

        a = 1
        b = 2*(ray.d*dv).sum(-1)
        c = (dv*dv).sum(-1)-self.radcurv**2
        
        qs  = b**2-4*a*c
   
        qsc = np.maximum(qs, 0)
        lamp = (-b+np.sqrt(qsc))/(2*a)
        lamn = (-b-np.sqrt(qsc))/(2*a)
        
        if self.radcurv > 0 :
            lam = np.minimum(lamp,lamn)
        elif self.radcurv < 0 :
            lam = np.maximum(lamp,lamn)
            
        lam = lam*(qs >= 0)
            
        # assign intersection
        ray.p1 = ray.propagate(lam)
    
    def surfaceNormal(self, p1) :
        cv = self.placement.location+self.placement.orientation*self.radcurv

        sn = np.sign(self.radcurv)*(cv-p1)
        sn = sn/np.linalg.norm(sn, axis=-1, keepdims=True)
        return sn

    def __str__(self) :
        s  = 'SphericalSurface         : '+self.name+'\n'
        s += 'SphericalSurface.volume  : \n' + Element.__str__(self) + '\n'
        s += 'SphericalSurface.radcurv : '+str(self.radcurv)+'\n'
        return s

############################################################################
# ParabolicSurface
############################################################################
class ParabolicSurface(Element) :
    def __init__(self,name,shape,dimension,placement,material,a,b) :
        pass
    
    def surfaceNorma(self, p1) :
        pass

############################################################################
# Cylindrical surface
############################################################################                
class CylindricalSurface(OpticalSurface) :
    def __init__(self,*args, **kwargs) :
        Element.__init__(self, *args, **kwargs)
        self.radcurv   = kwargs['curvature_radius']
        self.axiscurve = np.array(kwargs['curvature_axis'])

    def surface(self, proj=None) :
        xx, yy = self._surface_coordinates(proj)
        #cv = self.placement.location+self.placement.orientation*self.radcurv
        cv = self.placement.orientation*self.radcurv        
        #zz = -np.sign(self.radcurv)*np.sqrt(self.radcurv**2-(xx-cv[0])**2-(yy-cv[1])**2)+cv[2]
        zz = -(-np.sign(self.radcurv)*np.sqrt(self.radcurv**2-xx**2-yy**2) +self.radcurv)
        
        return self._orientate_surface(xx, yy, zz, proj)


    def intersection(self,ray) :
        #cv = self.placement.location+self.placement.orientation*self.radcurv
        dv = ray.p0 - self.placement.orientation*self.radcurv - self.placement.location
        #print dv.shape, self.axiscurve.shape
        dv_ = dv - (dv*self.axiscurve).sum(-1)[:,None]*self.axiscurve
        d_ = ray.d - (ray.d*self.axiscurve).sum(-1)[:,None]*self.axiscurve
        a = 1
        b = 2*(d_*dv_).sum(-1)
        c = (dv_*dv_).sum(-1)-self.radcurv**2

        qs = b ** 2 - 4 * a * c
        # if qs == 0 :
        #     lam = -b/(2*a)
        # elif qs < 0 :
        #     lam = None
        # else :

        qsc = np.maximum(qs, 0)
        lamp = (-b + np.sqrt(qsc)) / (2 * a)
        lamn = (-b - np.sqrt(qsc)) / (2 * a)
        #pd   = np.linalg.norm(ray.propagate(lamp)-ray.p0)
        #nd   = np.linalg.norm(ray.propagate(lamn)-ray.p0)
        #            lam = min(lamp,lamn)

        if self.radcurv > 0:
            lam = np.minimum(lamp, lamn)
        elif self.radcurv < 0:
            lam = np.maximum(lamp, lamn)

        lam = lam * (qs >= 0)
            
            # assign intersection
        ray.p1 = ray.propagate(lam)
    
    def surfaceNormal(self, p1) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
#        sn = p1-cv
#        sn = -sn/np.linalg.norm(sn)
        sn = np.sign(self.radcurv)*(cv-p1)
        sn = sn - (sn*self.axiscurve).sum(-1)[:,None]*self.axiscurve
        sn = sn/np.linalg.norm(sn, axis=-1, keepdims=True)
        return sn


############################################################################
# ThinLens
############################################################################            
class ThinLens(PlaneSurface) :
    def __init__(self,*args, **kwargs) :
        PlaneSurface.__init__(self,*args, **kwargs)
        self.focalLength  = kwargs['focal_length']

    def propagate(self,inray) :
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)
        
        #virtual 'in focus' object'
        o1 = inray.p1 - inray.d*self.focalLength/(inray.d*sn).sum(-1)[:,None]
        
        #output ray will be parallel to the ray from point through origin
        r = self.placement.location - o1
        
        d2 = r/np.linalg.norm(r)
    
        
        return RayBundle(inray.p1, d2, inray.material, inray.wavelength, inray.color,
                           cumulativePath=inray.prev_pathlength, intensities=inray.intensities)




############################################################################
# ThinLens (High NA)
############################################################################
class ThinLensH(ThinLens) :
    def __init__(self,*args, **kwargs) :
        ThinLens.__init__(self,*args, **kwargs)
        self.focalLength  = kwargs['focal_length']
        self.focalPoint = self.focalLength*self.orientation + self.location

    def propagate(self,inray) :
        # compute intersections
        self.intersection(inray)

        #virtual 'in focus' object'
        l = np.dot(self.focalPoint - inray.p0, self.placement.orientation)/np.dot(inray.d, self.placement.orientation)
        o1 = inray.p0 + l*inray.d
    
        
        if np.dot(inray.d, self.placement.orientation) <=0:   
            #output ray will be parallel to the ray from point through origin
            r = self.placement.location - o1
            
            d2 = r/np.linalg.norm(r)
        else:
            #find intercept with lens plane
            l = np.dot(self.placement.location - inray.p0, self.placement.orientation)/np.dot(inray.d, self.placement.orientation)
            o2 = inray.p0 + l*inray.d

            of = o1 - (o2-self.placement.location)

            #ray will run from focal plane curve through of
            r = of - inray.p1
            d2 = r/np.linalg.norm(r)
     
        
        return RayBundle(inray.p1, d2, inray.material, inray.wavelength, inray.color,
                           cumulativePath=inray.prev_pathlength, intensities=inray.intensities)

    def intersection(self,inray) :
        oc = inray.p0 - self.focalPoint
        doc = np.dot(inray.d, oc)
        sqterm = np.sqrt(doc*doc - np.dot(oc, oc) + self.focalLength*self.focalLength)
        lam1 = -doc + sqterm
        p11 = inray.propagate(lam1)
        lam2 = -doc - sqterm
        
        p12 = inray.propagate(lam2)

        if np.dot(p11 - self.focalPoint, self.placement.orientation) > 0:
            inray.p1 = p12
        else:
            inray.p1 = p11
    

############################################################################
# Circular aperture 
############################################################################
class Aperture(PlaneSurface) :
    def __init__(self,*args, **kwargs) :
        PlaneSurface.__init__(self,*args, **kwargs)
        self.radius = kwargs['radius']
        
   
    def propagate(self,inray) :
        # compute intersections
        self.intersection(inray)

        of = (inray.p1-self.placement.location)
        #print self.radius, np.linalg.norm(of)

        mask = np.linalg.norm(of, axis=-1) < self.radius

        return RayBundle(inray.p1, inray.d, inray.material, inray.wavelength, inray.color,
                           cumulativePath=inray.prev_pathlength, intensities=inray.intensities*mask)



##############################################
# Mirror predefined object
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

class Mirror(PlaneSurface):
    def __init__(self, *args, **kwargs):
        
        kwargs['material'] = material.MIRROR
        
        PlaneSurface.__init__(self, *args, **kwargs)


    @classmethod
    def to_rotate(cls, placement, axis=[0, 0, 1], angle=np.pi / 2., **kwargs):
        nm = np.dot(rotation_matrix(axis, angle / 2.), placement.orientation)
        return cls(placement.offset(orientation=nm), **kwargs)

    @classmethod
    def to_rotate_deg(cls, placement, axis=[0, 0, 1], angle=90., **kwargs):
        return cls.to_rotate(placement, axis=axis, angle=angle * np.pi / 180., **kwargs)

        
        



############################################################################
# Snells' law 
############################################################################            
def snell(ray,sn,material1,material2) :
    n1 = material1.n(ray.wavelength)
    n2 = material2.n(ray.wavelength)

    nr = n1/n2    
    ct1 = -(sn*ray.d).sum(-1)
    ct2 = np.sqrt(1-(nr**2)*(1-ct1**2))

    ct2 = ct2*np.sign(ct1)
    
    d2 = nr*ray.d+np.atleast_1d(nr*ct1-ct2)[:,None]*sn
    r = RayBundle(ray.p1, d2, material2, ray.wavelength, ray.color,
                           cumulativePath=ray.prev_pathlength, intensities=ray.intensities)

    return r
    
############################################################################
# Reflection law
############################################################################            
def reflect(ray,sn) :
    d2 = ray.d-2*np.atleast_1d((ray.d*sn).sum(-1))[:,None]*sn
    
    #print d2
    r = RayBundle(ray.p1, d2, ray.material, ray.wavelength, ray.color,
                           cumulativePath=ray.prev_pathlength, intensities=ray.intensities)
    
    #print r.d
    return r


import numpy as np

from material  import *
from placement import *
from rays       import *

class Volume() :
    ''' Orientable volume, base class to describe extent, material and orientation of a volume containing ''' 

    # shape of element
    rect = 0
    circ = 1 
            
    def __init__(self,name, shape,dimension,placement,material) :
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
        elif self.shape == self.rect:
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
        
        #ang = np.arccos(np.dot(N, self.placement.orientation))
    
        #ct = np.cos(ang)
        #st = np.sin(ang)
    
        #print ax, ang
        
        coords = self.placement.orientation[:,None,None]*zz[None,:, :] + ax[:,None,None]*xx[None,:,:] + N2[:,None, None]*yy[None,:,:]
        return coords + self.placement.location[:,None, None]
    
        # xn = (ct + ax[0] ** 2 * (1 - ct)) * xx + (ax[0] * ax[1] * (1 - ct) - ax[2] * st) * yy + (ax[0] * ax[2] * (
        # 1 - ct) + ax[1] * st) * zz
        # yn = (ax[0] * ax[1] * (1 - ct) + ax[2] * st) * xx + (ct + ax[1] ** 2 * (1 - ct)) * yy + (ax[1] * ax[2] * (
        # 1 - ct) + ax[0] * st) * zz
        # zn = (ax[0] * ax[2] * (1 - ct) - ax[1] * st) * xx + (ax[2] * ax[1] * (1 - ct) + ax[0] * st) * yy + (ct + ax[
        #     2] ** 2 * (1 - ct)) * zz
        #
        # return [xn + self.placement.location[0], yn + self.placement.location[1], zn + self.placement.location[2]]
        

############################################################################
# General optical surface
############################################################################            
class OpticalSurface(Volume) :
    def __init__(self) :
        pass

    def propagate(self,previous,inray) :
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        # compute out going ray
        
        if self.material.type == Material.mirror :
            outray = reflect(inray,sn)
        elif self.material.type == Material.refract :
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
    def __init__(self,name,shape,dimension,placement,material) :
        Volume.__init__(self,name,shape,dimension,placement,material)

    def surface(self, proj=None) :
        xx, yy = self._surface_coordinates(proj)
            
        zz = np.zeros_like(xx)
        
        return self._orientate_surface(xx, yy, zz, proj)
        
        # N = np.array([0,0,1.])
        # ax = np.cross(N, self.placement.orientation)
        # ax = ax/np.linalg.norm(ax)
        # ang = np.arccos(np.dot(N, self.placement.orientation))
        #
        # ct = np.cos(ang)
        # st = np.sin(ang)
        #
        # #print ax, ang
        #
        # xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        # yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        # zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        #
        # return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]
        

    def intersection(self,inray) :
#        lam = np.linalg.dot(self.placement.location,inray.p0)/np.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(((self.placement.location-inray.p0)*(self.placement.orientation)).sum(-1)/(inray.d* self.placement.orientation).sum(-1))
        #print 'PlaneSurface.intersection',lam
        inray.p1 = inray.propagate(lam)
                
    def surfaceNormal(self,inray) :
        return self.placement.orientation
    
    
    def __str__(self) :
        s  = 'PlaneSurface             : '+self.name+'\n'
        s += 'PlaneSurface.Volume      :\n'+Volume.__str__(self)
        return s

############################################################################
# SphericalSurface
############################################################################            
class SphericalSurface(OpticalSurface) :
    def __init__(self,name,shape,dimension,placement,material,radcurv) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radcurv   = radcurv

    def surface(self, proj=None) :
        xx, yy = self._surface_coordinates(proj)
        #cv = self.placement.location+self.placement.orientation*self.radcurv
        cv = self.placement.orientation*self.radcurv        
        #zz = -np.sign(self.radcurv)*np.sqrt(self.radcurv**2-(xx-cv[0])**2-(yy-cv[1])**2)+cv[2]
        zz = np.sign(self.radcurv)*(-np.sqrt(self.radcurv**2-xx**2-yy**2) +np.abs(self.radcurv))
        #print zz.shape, zz, xx, yy
        
        #return [xx,yy,zz]
        
        return self._orientate_surface(xx, yy, zz, proj)

    def intersection(self,ray) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
        dv = ray.p0 - self.placement.orientation*self.radcurv - self.placement.location
        #print dv.shape, ray.d.shape
        a = 1
        b = 2*(ray.d*dv).sum(-1)
        c = (dv*dv).sum(-1)-self.radcurv**2
        
        #print c
        
        qs  = b**2-4*a*c
        # if qs == 0 :
        #     lam = -b/(2*a)
        # elif qs < 0 :
        #     lam = None
        # else :
            
        qsc = np.maximum(qs, 0)
        lamp = (-b+np.sqrt(qsc))/(2*a)
        lamn = (-b-np.sqrt(qsc))/(2*a)
        #pd   = np.linalg.norm(ray.propagate(lamp)-ray.p0)
        #nd   = np.linalg.norm(ray.propagate(lamn)-ray.p0)
#            lam = min(lamp,lamn)
        
        if self.radcurv > 0 :
            lam = np.minimum(lamp,lamn)
        elif self.radcurv < 0 :
            lam = np.maximum(lamp,lamn)
            
        lam = lam*(qs >= 0)
            
            # assign intersection
        ray.p1 = ray.propagate(lam)
    
    def surfaceNormal(self, p1) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
#        sn = p1-cv
#        sn = -sn/np.linalg.norm(sn)
        sn = np.sign(self.radcurv)*(cv-p1)
        sn = sn/np.linalg.norm(sn, axis=-1, keepdims=True)
        return sn

    def __str__(self) :
        s  = 'SphericalSurface         : '+self.name+'\n'
        s += 'SphericalSurface.volume  : \n'+Volume.__str__(self)+'\n'
        s += 'SphericalSurface.radcurv : '+str(self.radcurv)+'\n'
        return s

############################################################################
# ParabolicSurface
############################################################################
class ParabolicSurface(Volume) :
    def __init__(self,name,shape,dimension,placement,material,a,b) :
        pass
    
    def surfaceNorma(self, p1) :
        pass

############################################################################
# Cylindrical surface
############################################################################                
class CylindricalSurface(OpticalSurface) :
    def __init__(self,name,shape,dimension,placement,material,radcurv, axiscurve) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radcurv   = radcurv
        self.axiscurve = np.array(axiscurve)

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
# EvenAsphericalSurface
############################################################################            
class EvenAsphericalSurface(Volume) :
    def __init__(self,volume) :
        print 'AsphericalSurface.__init__>'
        self.volume = volume

############################################################################
# OddAsphericalSurface
############################################################################            
class OddAsphericalSurface(Volume) :
    def __init__(self,volume) :
        print 'AsphericalSurface.__init__>'
        self.volume = volume

############################################################################
# ThinLens
############################################################################            
class ThinLens(Volume) :
    def __init__(self,name,shape,dimension,placement,material,focalLength) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.focalLength   = focalLength
        
    def surface(self, proj=None):
        xx, yy = self._surface_coordinates(proj)
        #zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = np.zeros_like(xx)
        
        return self._orientate_surface(xx, yy, zz, proj)

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        #if self.material.type == Material.mirror :
        #    outray = reflect(inray,sn)
        #elif self.material.type == Material.refract :
        #    outray = snell(inray,sn,previous.material,self.material)  
            
        #r = inray.p1 - self.placement.location
        
        print inray.p1.shape, inray.d.shape, sn.shape
        
        #virtual 'in focus' object'
        o1 = inray.p1 - inray.d*self.focalLength/(inray.d*sn).sum(-1)[:,None]
        
        #output ray will be parallel to the ray from point through origin
        r = self.placement.location - o1
        
        d2 = r/np.linalg.norm(r)
        
         
#        sn = self.surfaceNormal(inray)
#        mr = np.linalg.norm(r)
#        rhat = r/mr
#        d1x = np.dot(inray.d, sn)
#        d1y = np.dot(inray.d, rhat)
#        f = self.focalLength
#        c = np.sqrt(d1x**2 + d1y**2)
#        mrft = mr/f - d1y/d1x
#        d2y = -np.sign(mrft)*c/np.sqrt(1./mrft**2 + 1)
#        d2x = -np.sign(mrft)*np.sqrt(c**2 - d2y**2)
#        #print d2x, d1x, 
#        
#        d2 = inray.d + (d2y-d1y)*rhat + (d2x-d1x)*sn
#        
#        if (mr < 1e-6):
#            #catch rays passing through center
#            d2 = inray.d
#        
        #print inray.p1, o1, inray.d, d2
#        print d1x, d2x, d1y, d2y
        
        outray = RayBundle(inray.p1, d2, inray.material, inray.wavelength, inray.color, cumulativePath=inray.prev_pathlength)

        # compute out going ray


        return outray

    def intersection(self,inray) :
#        lam = np.linalg.dot(self.placement.location,inray.p0)/np.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(((self.placement.location-inray.p0)*(self.placement.orientation)).sum(-1)/(self.placement.orientation*inray.d).sum(-1))
        #print 'PlaneSurface.intersection',lam
        inray.p1 = inray.propagate(lam)
    
    def surfaceNormal(self,inray) :
        return self.placement.orientation

############################################################################
# ThinLens (High NA)
############################################################################
class ThinLensH(Volume) :
    def __init__(self,name,shape,dimension,placement,material,focalLength) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.focalLength   = focalLength
        self.focalPoint = focalLength*placement.orientation + placement.location
        
    def surface(self, proj=None):
        xx, yy = self._surface_coordinates(proj)
        #zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = np.zeros_like(xx)
    
        return self._orientate_surface(xx, yy, zz, proj)

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersections
        self.intersection(inray)
        

        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        #if self.material.type == Material.mirror :
        #    outray = reflect(inray,sn)
        #elif self.material.type == Material.refract :
        #    outray = snell(inray,sn,previous.material,self.material)  
            
        #r = inray.p1 - self.placement.location
        #virtual 'in focus' object'
        l = np.dot(self.focalPoint - inray.p0, self.placement.orientation)/np.dot(inray.d, self.placement.orientation)
        o1 = inray.p0 + l*inray.d
    
        #o1 = inray.p1 - inray.d*self.focalLength/np.dot(inray.d, sn)
        
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
    
         
#        sn = self.surfaceNormal(inray)
#        mr = np.linalg.norm(r)
#        rhat = r/mr
#        d1x = np.dot(inray.d, sn)
#        d1y = np.dot(inray.d, rhat)
#        f = self.focalLength
#        c = np.sqrt(d1x**2 + d1y**2)
#        mrft = mr/f - d1y/d1x
#        d2y = -np.sign(mrft)*c/np.sqrt(1./mrft**2 + 1)
#        d2x = -np.sign(mrft)*np.sqrt(c**2 - d2y**2)
#        #print d2x, d1x, 
#        
#        d2 = inray.d + (d2y-d1y)*rhat + (d2x-d1x)*sn
#        
#        if (mr < 1e-6):
#            #catch rays passing through center
#            d2 = inray.d
#        
        #print inray.p1, o1, inray.d, d2
#        print d1x, d2x, d1y, d2y
        
        outray = RayBundle(inray.p1, d2, inray.material, inray.wavelength, inray.color, cumulativePath=inray.prev_pathlength)

        # compute out going ray


        return outray

    def intersection(self,inray) :
#        lam = np.linalg.dot(self.placement.location,inray.p0)/np.linalg.dot(self.placement.orientation,self.placement.location)
        #lam = abs(np.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/np.linalg.dot(self.placement.orientation,inray.d))
        oc = inray.p0 - self.focalPoint
        doc = np.dot(inray.d, oc)
        sqterm = np.sqrt(doc*doc - np.dot(oc, oc) + self.focalLength*self.focalLength)
        lam1 = -doc + sqterm
        p11 = inray.propagate(lam1)
        lam2 = -doc - sqterm
        #print 'PlaneSurface.intersection',lam
        p12 = inray.propagate(lam2)

        if np.dot(p11 - self.focalPoint, self.placement.orientation) > 0:
            inray.p1 = p12
        else:
            inray.p1 = p11
    
    def surfaceNormal(self,inray) :
        return self.placement.orientation

############################################################################
# Circular aperture 
############################################################################
class Aperture(Volume) :
    def __init__(self,name,shape,dimension,placement,material,radius) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radius = radius
        
    def surface(self, proj=None) :
        xx, yy, zz = self._surface_coordinates(proj)
        #zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = np.zeros_like(xx)
    
        return self._orientate_surface(xx, yy, zz, proj)

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersections
        self.intersection(inray)

        of = (inray.p1-self.placement.location)
        #print self.radius, np.linalg.norm(of)

        if np.linalg.norm(of) < self.radius:
            outray = RayBundle(inray.p1, inray.d, inray.material, inray.wavelength, inray.color, cumulativePath=inray.cumulativePath)

            # compute out going ray
            return outray
        else:
            return None

    def intersection(self,inray) :
#        lam = np.linalg.dot(self.placement.location,inray.p0)/np.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(np.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/np.linalg.dot(self.placement.orientation,inray.d))
        
        inray.p1 = inray.propagate(lam)
    
    def surfaceNormal(self,inray) :
        return self.placement.orientation

############################################################################
# Circular aperture 
############################################################################
class AtAperture(Volume) :
    def __init__(self,name,shape,dimension,placement,material,radius) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radius = radius
        
    def surface(self, proj=None) :
        xx, yy, zz = self._surface_coordinates(proj)
        #zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = np.zeros_like(xx)
    
        return self._orientate_surface(xx, yy, zz, proj)

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersections
        self.intersection(inray)

        of = (inray.p1-self.placement.location)
        #print self.radius, np.linalg.norm(of)

        if np.linalg.norm(of) < self.radius:
            outray = RayBundle(inray.p1, inray.d, inray.material, inray.wavelength, inray.color, cumulativePath=inray.cumulativePath)

            # compute out going ray
            return outray
        else:
            outray = RayBundle(inray.p1, inray.d, inray.material, inray.wavelength, np.array(inray.color) * .2, cumulativePath=inray.cumulativePath)

            # compute out going ray
            return outray
            return None

    def intersection(self,inray) :
#        lam = np.linalg.dot(self.placement.location,inray.p0)/np.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(np.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/np.linalg.dot(self.placement.orientation,inray.d))
        
        inray.p1 = inray.propagate(lam)
    
    def surfaceNormal(self,inray) :
        return self.placement.orientation

############################################################################
# Snells' law 
############################################################################            
def snell(ray,sn,material1,material2) :
    n1 = material1.n(ray.wavelength)
    n2 = material2.n(ray.wavelength)

#    nr = n1/n2
#    dp  = np.dot(ray.d,sn)
#    gam = ((nr*dp)**2-(n1/n2)+1)**0.5 - nr*dp
#    # new direction
#    d2  = ray.d*nr+gam*sn
#    d2 = d2/np.linalg.norm(d2)
 #   r = Ray(ray.p1,d2)
    nr = n1/n2    
    ct1 = -(sn*ray.d).sum(-1)
    ct2 = np.sqrt(1-(nr**2)*(1-ct1**2))
    #if ct1 < 0 :
    #    ct2 = -ct2
    ct2 = ct2*np.sign(ct1)
    
    d2 = nr*ray.d+np.atleast_1d(nr*ct1-ct2)[:,None]*sn
    r = RayBundle(ray.p1, d2, material2, ray.wavelength, ray.color, cumulativePath=ray.cumulativePath)
#    print 'snell> in\n',ray
#    print 'snell> normal ',sn
#    print 'snell> material 1',material1
#    print 'snell> material 2',material2
#    print 'snell> cos theta',ct1
#    print 'snell> cos theta',ct2    
#    print 'snell> out',r
    return r
    
############################################################################
# Reflection law
############################################################################            
def reflect(ray,sn) :
    d2 = ray.d-2*np.atleast_1d((ray.d*sn).sum(-1))[:,None]*sn
    
    #print d2
    r = RayBundle(ray.p1, d2, ray.material, ray.wavelength, ray.color, cumulativePath=ray.cumulativePath)
    
    #print r.d
    return r


import pylab as pl
import numpy as np

#from numpy.linalg import dot

from Material  import *
from Placement import *
from Ray       import *

class Volume() :
    ''' Orientable volume, base class to describe extent, material and orientation of a volume containing ''' 

    # shape of element
    rect = 0
    circ = 1 
            
    def __init__(self,name, shape,dimension,placement,material) :
        self.name      = name
        self.shape     = shape
        self.dimension = pl.array(dimension)
        self.placement = placement
        self.material  = material

    def  __str__(self) :
        s  = 'Volume\n'
        s += 'Volume.shape             : '+str(self.shape)+'\n'
        s += 'Volume.dimension         : '+str(self.dimension)+'\n'
        s += 'Volume.placement         : \n'+str(self.placement)
        s += 'Volume.material          : \n'+str(self.material)
        return s

############################################################################
# General optical surface
############################################################################            
class OpticalSurface(Volume) :
    def __init__(self) :
        pass

    def proagate(self,previous,inray) :
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        # compute out going ray
        outray = snell(inray,sn,previous.material,self.material)
        
    def intersection(self,inray) :
        pass
    def surfaceNormal(self,inray) :
        pass

############################################################################
# Plane surface
############################################################################            
class PlaneSurface(Volume) :
    def __init__(self,name,shape,dimension,placement,material) :
        Volume.__init__(self,name,shape,dimension,placement,material)

    def surface(self) :
#        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
#        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
#        xx,yy = pl.meshgrid(x,y)
#        zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
#        
#        return [xx + self.placement.location[0],yy+ self.placement.location[1],zz]  
        
        if self.shape == self.rect:
            x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
            y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
            xx,yy = pl.meshgrid(x,y)
        else:
            r = self.dimension[0]*pl.sqrt(pl.arange(0, 1.01,1))
            t =  2*pl.pi*pl.arange(0, 1.01,.1)
            
            rr, tt = pl.meshgrid(r,t)
            pol = rr*pl.exp(1j*tt)
            xx = pol.real
            yy = pol.imag
            
        zz = pl.zeros_like(xx)
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ax = ax/pl.norm(ax)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]
        
    def propagate(self,previous,inray) :        
        
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        
        if self.material.type == Material.mirror :
            outray = reflect(inray,sn)
        elif self.material.type == Material.refract :
            outray = snell(inray,sn,previous.material,self.material)            

        return outray

    def intersection(self,inray) :
#        lam = pl.linalg.dot(self.placement.location,inray.p0)/pl.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d))
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
class SphericalSurface(Volume) :
    def __init__(self,name,shape,dimension,placement,material,radcurv) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radcurv   = radcurv

    def surface(self) :
        if self.shape == self.rect:
            x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
            y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
            xx,yy = pl.meshgrid(x,y)
        else:
            r = self.dimension[0]*pl.sqrt(pl.arange(0, 1.01,.2))
            t =  2*pl.pi*pl.arange(0, 1.01,.1)
            
            rr, tt = pl.meshgrid(r,t)
            pol = rr*pl.exp(1j*tt)
            xx = pol.real
            yy = pol.imag
        #cv = self.placement.location+self.placement.orientation*self.radcurv
        cv = self.placement.orientation*self.radcurv        
        #zz = -pl.sign(self.radcurv)*pl.sqrt(self.radcurv**2-(xx-cv[0])**2-(yy-cv[1])**2)+cv[2]
        zz = pl.sign(self.radcurv)*(-pl.sqrt(self.radcurv**2-xx**2-yy**2) +np.abs(self.radcurv))
        #print zz.shape, zz, xx, yy
        
        #return [xx,yy,zz]
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz

        #print xn, yn, zn
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        if self.material.type == Material.mirror :
            outray = reflect(inray,sn)
        elif self.material.type == Material.refract :
            outray = snell(inray,sn,previous.material,self.material)            

        # compute out going ray


        return outray

    def intersection(self,ray) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
        dv = ray.p0 - self.placement.orientation*self.radcurv - self.placement.location        
        a = 1
        b = 2*pl.linalg.dot(ray.d,dv)
        c = pl.linalg.dot(dv,dv)-self.radcurv**2
        
        qs  = b**2-4*a*c
        if qs == 0 :
            lam = -b/(2*a)
        elif qs < 0 :
            lam = None
        else :
            lamp = (-b+pl.sqrt(b**2-4*a*c))/(2*a)
            lamn = (-b-pl.sqrt(b**2-4*a*c))/(2*a)
            pd   = pl.linalg.norm(ray.propagate(lamp)-ray.p0)
            nd   = pl.linalg.norm(ray.propagate(lamn)-ray.p0)
#            lam = min(lamp,lamn)
            
            if self.radcurv > 0 :
                lam = min(lamp,lamn)
            elif self.radcurv < 0 :
                lam = max(lamp,lamn)
            
            # assign intersection
        ray.p1 = ray.propagate(lam)
    
    def surfaceNormal(self, p1) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
#        sn = p1-cv
#        sn = -sn/pl.linalg.norm(sn)
        sn = pl.sign(self.radcurv)*(cv-p1)
        sn = sn/pl.linalg.norm(sn)
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
class CylindricalSurface(Volume) :
    def __init__(self,name,shape,dimension,placement,material,radcurv, axiscurve) :
        Volume.__init__(self,name,shape,dimension,placement,material)
        self.radcurv   = radcurv
        self.axiscurve = axiscurve

    def surface(self) :
        if self.shape == self.rect:
            x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
            y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
            xx,yy = pl.meshgrid(x,y)
        else:
            r = self.dimension[0]*pl.sqrt(pl.arange(0, 1.01,.2))
            t =  2*pl.pi*pl.arange(0, 1.01,.1)
            
            rr, tt = pl.meshgrid(r,t)
            pol = rr*pl.exp(1j*tt)
            xx = pol.real
            yy = pol.imag
        #cv = self.placement.location+self.placement.orientation*self.radcurv
        cv = self.placement.orientation*self.radcurv        
        #zz = -pl.sign(self.radcurv)*pl.sqrt(self.radcurv**2-(xx-cv[0])**2-(yy-cv[1])**2)+cv[2]
        zz = -(-pl.sign(self.radcurv)*pl.sqrt(self.radcurv**2-xx**2-yy**2) +self.radcurv)
        
        #return [xx,yy,zz]
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersection
        self.intersection(inray)                    
        # compute normal
        sn = self.surfaceNormal(inray.p1)            
        if self.material.type == Material.mirror :
            outray = reflect(inray,sn)
        elif self.material.type == Material.refract :
            outray = snell(inray,sn,previous.material,self.material)            

        # compute out going ray


        return outray

    def intersection(self,ray) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
        dv = ray.p0 - self.placement.orientation*self.radcurv - self.placement.location        
        a = 1
        b = 2*pl.linalg.dot(ray.d,dv)
        c = pl.linalg.dot(dv,dv)-self.radcurv**2
        
        qs  = b**2-4*a*c
        if qs == 0 :
            lam = -b/(2*a)
        elif qs < 0 :
            lam = None
        else :
            lamp = (-b+pl.sqrt(b**2-4*a*c))/(2*a)
            lamn = (-b-pl.sqrt(b**2-4*a*c))/(2*a)
            pd   = pl.linalg.norm(ray.propagate(lamp)-ray.p0)
            nd   = pl.linalg.norm(ray.propagate(lamn)-ray.p0)
#            lam = min(lamp,lamn)
            
            if self.radcurv > 0 :
                lam = min(lamp,lamn)
            elif self.radcurv < 0 :
                lam = max(lamp,lamn)
            
            # assign intersection
        ray.p1 = ray.propagate(lam)
    
    def surfaceNormal(self, p1) :
        cv = self.placement.location+self.placement.orientation*self.radcurv
#        sn = p1-cv
#        sn = -sn/pl.linalg.norm(sn)
        sn = pl.sign(self.radcurv)*(cv-p1)
        sn = sn/pl.linalg.norm(sn)
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
        
    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        #zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = pl.ones_like(xx)
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

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
        
        #virtual 'in focus' object'
        o1 = inray.p1 - inray.d*self.focalLength/pl.dot(inray.d, sn)
        
        #output ray will be parallel to the ray from point through origin
        r = self.placement.location - o1
        
        d2 = r/pl.norm(r)
        
         
#        sn = self.surfaceNormal(inray)
#        mr = pl.norm(r)
#        rhat = r/mr
#        d1x = pl.dot(inray.d, sn)
#        d1y = pl.dot(inray.d, rhat)
#        f = self.focalLength
#        c = pl.sqrt(d1x**2 + d1y**2)
#        mrft = mr/f - d1y/d1x
#        d2y = -pl.sign(mrft)*c/pl.sqrt(1./mrft**2 + 1)
#        d2x = -pl.sign(mrft)*pl.sqrt(c**2 - d2y**2)
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
        
        outray = Ray(inray.p1,d2, inray.material, inray.wavelength, inray.color)

        # compute out going ray


        return outray

    def intersection(self,inray) :
#        lam = pl.linalg.dot(self.placement.location,inray.p0)/pl.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d))
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
        
    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        #zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = pl.ones_like(xx)
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

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
        l = pl.dot(self.focalPoint - inray.p0, self.placement.orientation)/pl.dot(inray.d, self.placement.orientation)
        o1 = inray.p0 + l*inray.d
    
        #o1 = inray.p1 - inray.d*self.focalLength/pl.dot(inray.d, sn)
        
        if pl.dot(inray.d, self.placement.orientation) <=0:   
            #output ray will be parallel to the ray from point through origin
            r = self.placement.location - o1
            
            d2 = r/pl.norm(r)
        else:
            #find intercept with lens plane
            l = pl.dot(self.placement.location - inray.p0, self.placement.orientation)/pl.dot(inray.d, self.placement.orientation)
            o2 = inray.p0 + l*inray.d

            of = o1 - (o2-self.placement.location)

            #ray will run from focal plane curve through of
            r = of - inray.p1
            d2 = r/pl.norm(r)
    
         
#        sn = self.surfaceNormal(inray)
#        mr = pl.norm(r)
#        rhat = r/mr
#        d1x = pl.dot(inray.d, sn)
#        d1y = pl.dot(inray.d, rhat)
#        f = self.focalLength
#        c = pl.sqrt(d1x**2 + d1y**2)
#        mrft = mr/f - d1y/d1x
#        d2y = -pl.sign(mrft)*c/pl.sqrt(1./mrft**2 + 1)
#        d2x = -pl.sign(mrft)*pl.sqrt(c**2 - d2y**2)
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
        
        outray = Ray(inray.p1,d2, inray.material, inray.wavelength, inray.color)

        # compute out going ray


        return outray

    def intersection(self,inray) :
#        lam = pl.linalg.dot(self.placement.location,inray.p0)/pl.linalg.dot(self.placement.orientation,self.placement.location)
        #lam = abs(pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d))
        oc = inray.p0 - self.focalPoint
        doc = pl.dot(inray.d, oc)
        sqterm = np.sqrt(doc*doc - pl.dot(oc, oc) + self.focalLength*self.focalLength)
        lam1 = -doc + sqterm
        p11 = inray.propagate(lam1)
        lam2 = -doc - sqterm
        #print 'PlaneSurface.intersection',lam
        p12 = inray.propagate(lam2)

        if pl.dot(p11 - self.focalPoint, self.placement.orientation) > 0:
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
        
    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        #zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = pl.ones_like(xx)
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersections
        self.intersection(inray)

        of = (inray.p1-self.placement.location)
        #print self.radius, pl.norm(of)

        if pl.norm(of) < self.radius:        
            outray = Ray(inray.p1,inray.d, inray.material, inray.wavelength, inray.color)

            # compute out going ray
            return outray
        else:
            return None

    def intersection(self,inray) :
#        lam = pl.linalg.dot(self.placement.location,inray.p0)/pl.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d))
        
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
        
    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        #zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        zz = pl.ones_like(xx)
        
        N = pl.array([0,0,1.])         
        ax = pl.cross(N, self.placement.orientation)
        ang = pl.arccos(pl.dot(N, self.placement.orientation))
        
        ct = pl.cos(ang)
        st = pl.sin(ang)
        
        #print ax, ang
        
        xn = (ct + ax[0]**2*(1-ct))*xx           + (ax[0]*ax[1]*(1-ct) - ax[2]*st)*yy + (ax[0]*ax[2]*(1-ct) + ax[1]*st)*zz
        yn = (ax[0]*ax[1]*(1-ct) + ax[2]*st)*xx + (ct + ax[1]**2*(1-ct))*yy           + (ax[1]*ax[2]*(1-ct) + ax[0]*st)*zz
        zn = (ax[0]*ax[2]*(1-ct) - ax[1]*st)*xx + (ax[2]*ax[1]*(1-ct) + ax[0]*st)*yy + (ct + ax[2]**2*(1-ct))*zz
        
        return [xn+self.placement.location[0],yn+self.placement.location[1],zn+self.placement.location[2]]   

    def propagate(self,previous,inray) :
        outrays = []
        
        # compute intersections
        self.intersection(inray)

        of = (inray.p1-self.placement.location)
        #print self.radius, pl.norm(of)

        if pl.norm(of) < self.radius:        
            outray = Ray(inray.p1,inray.d, inray.material, inray.wavelength, inray.color)

            # compute out going ray
            return outray
        else:
            outray = Ray(inray.p1,inray.d, inray.material, inray.wavelength, np.array(inray.color)*.2)

            # compute out going ray
            return outray
            return None

    def intersection(self,inray) :
#        lam = pl.linalg.dot(self.placement.location,inray.p0)/pl.linalg.dot(self.placement.orientation,self.placement.location)
        lam = abs(pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d))
        
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
#    dp  = pl.dot(ray.d,sn)
#    gam = ((nr*dp)**2-(n1/n2)+1)**0.5 - nr*dp
#    # new direction
#    d2  = ray.d*nr+gam*sn
#    d2 = d2/pl.linalg.norm(d2)
 #   r = Ray(ray.p1,d2)
    nr = n1/n2    
    ct1 = -pl.linalg.dot(sn,ray.d)
    ct2 = pl.sqrt(1-(nr**2)*(1-ct1**2))
    if ct1 < 0 :
        ct2 = -ct2
    d2 = nr*ray.d+(nr*ct1-ct2)*sn
    r = Ray(ray.p1,d2, material2, ray.wavelength, ray.color)
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
    d2 = ray.d-2*pl.dot(ray.d,sn)*sn
    r = Ray(ray.p1,d2, ray.material, ray.wavelength, ray.color)
    return r


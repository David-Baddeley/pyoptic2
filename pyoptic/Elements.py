import pylab as pl

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
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        
        return [xx,yy,zz]        
        
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
        lam = pl.linalg.dot(self.placement.location-inray.p0,self.placement.orientation)/pl.linalg.dot(self.placement.orientation,inray.d)
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
        cv = self.placement.location+self.placement.orientation*self.radcurv        
        zz = -pl.sign(self.radcurv)*pl.sqrt(self.radcurv**2-(xx-cv[0])**2-(yy-cv[1])**2)+cv[2]
        
        return [xx,yy,zz]

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
    def __init__(self,volume) :
        print 'CylindricalSurface.__init__>'
        self.volume = volume

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
    def __init__(self,volume) :
        print 'ThisLens.__init__>'

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
    r = Ray(ray.p1,d2, material2)
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
    r = Ray(ray.p1,d2)
    return r


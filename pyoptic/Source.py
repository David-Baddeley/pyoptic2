import pylab as pl
import numpy as np

from Elements import *
from Ray import *
from Placement import *

class Source(Volume,list) :
    def __init__(self,name,placement, wavelength=635.) :
        print 'Source.__init__'
        self.name = name
        self.placement = placement
        self.dimension = pl.array([0.05,0.05,0.01])
        self.material = Material(Material.refract,1.0)
        self.wavelength=wavelength

    def __str__(self) :
        s  = 'Source                   : '+self.name+'\n'
        s += 'Source.placement         : \n'+str(self.placement)
        return s

    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        
        return [xx,yy,zz]        
    
    def exampleRays(self,d) :
        r0  = Ray(self.placement.location,[0,0,1], self.material)
        ryp = Ray(self.placement.location+pl.array([0,d,0]),[0,0,1], self.material)
        rypp = Ray(self.placement.location+pl.array([0,d/2,0]),[0,0,1], self.material)
        ryn = Ray(self.placement.location+pl.array([0,-d,0]),[0,0,1], self.material)
        rynn = Ray(self.placement.location+pl.array([0,-d/2,0]),[0,0,1], self.material)
        rxp = Ray(self.placement.location+pl.array([d,0,0]),[0,0,1], self.material)
        rxn = Ray(self.placement.location+pl.array([-d,0,0]),[0,0,1], self.material)
        self.append(r0)
        self.append(ryp)
        self.append(rypp)
        self.append(rynn)
        self.append(ryn)
        self.append(rxp)
        self.append(rxn)

#        r0p  = Ray(self.placement.location,[ 0.        ,  0.04993762,  0.99875234])
#        rypp = Ray(self.placement.location+pl.array([0,d,0]),[ 0.        ,  0.04993762,  0.99875234])
#        rynp = Ray(self.placement.location+pl.array([0,-d,0]),[ 0.        ,  0.04993762,  0.99875234])
#        rxpp = Ray(self.placement.location+pl.array([d,0,0]),[ 0.        ,  0.04993762,  0.99875234])
#        rxnp = Ray(self.placement.location+pl.array([-d,0,0]),[ 0.        ,  0.04993762,  0.99875234])
#        self.append(r0p)
#        self.append(rypp)
#        self.append(rynp)
#        self.append(rxpp)
#        self.append(rxnp)

class PointSource(Source) :
    def __init__(self,name,placement, NA=1.49/1.51, color=(1.0, 0, 0), wavelength=635.) :
        print 'Source.__init__'
        self.name = name
        self.placement = placement
        self.dimension = pl.array([0.05,0.05,0.01])
        self.material = Material(Material.refract,1.0)
        self.NA = NA
        self.color = color
        self.wavelength=wavelength

    def __str__(self) :
        s  = 'Source                   : '+self.name+'\n'
        s += 'Source.placement         : \n'+str(self.placement)
        return s

    def surface(self) :
        x = pl.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = pl.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = pl.meshgrid(x,y)
        zz = 1/self.placement.orientation[2]*(pl.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        
        return [xx,yy,zz]
    
    #def _mray(self, theta, phi):
    
    @property
    def chief_rays(self):
        '''returns the principal ray + marginal rays at full NA along each of the source axes.
        
        suitable for visualizing the beam extent
        '''
        return self._generate_rays(5, 2)
    
    @property
    def pupil_rays(self):
        '''returns a  set of rays which densely sample the pupil for abberation estimation'''
        return  self._generate_rays(50,10)
    
    def _generate_rays(self, nph, nth, jit=False):
        d2 = self.placement.orientation
        #d1 = pl.cross(d2, pl.array([1,0,0]))
        d1 = pl.cross(d2, pl.array([1, 1, 1]) - d2)
        d1 = d1 / pl.norm(d1)
        d0 = pl.cross(d1, d2)
        d0 = d0 / pl.norm(d0)
        
        rays = []
    
        def mray(theta, phi):
            r = np.exp(1j * phi)
            #print r.real, r.imag, np.cos(theta)
            dn = np.sin(theta) * r.real * d0 + np.sin(theta) * r.imag * d1 + np.cos(theta) * d2
            #print pl.norm(dn)
            #dn = dn/pl.norm(dn)
        
            return Ray(self.placement.location, dn, self.material, color=self.color, wavelength=self.wavelength)
    
        r0 = Ray(self.placement.location, d2, self.material, color=self.color, wavelength=self.wavelength)
        #        ryp = Ray(self.placement.location,y1*d1 + z1*d2, self.material, color=self.color)
        #        rypp = Ray(self.placement.location,y2*d1 + z2*d2, self.material, color=self.color)
        #        ryn = Ray(self.placement.location,-y1*d1 + z1*d2, self.material, color=self.color)
        #        rynn = Ray(self.placement.location,-y2*d1 + z2*d2, self.material, color=self.color)
        #        rxp = Ray(self.placement.location,y1*d0 + z1*d2, self.material, color=self.color)
        #        rxn = Ray(self.placement.location,-y1*d0 + z1*d2, self.material, color=self.color)
        rays.append(r0)
        #        self.append(ryp)
        #        self.append(rypp)
        #        self.append(rynn)
        #        self.append(ryn)
        #        self.append(rxp)
        #        self.append(rxn)
        for th in np.linspace(0, 1, nth)[1:]:
            for phi in (np.linspace(0, 2 * np.pi, np.maximum(nph * th, 3))[:-1] + 3 * np.pi / 4):
                if jit:
                    phi = phi + (np.random.rand(1) - .5) * 2 * np.pi / (nph * th)
                    thm = th + (np.random.rand(1) - .5) / nth
                else:
                    thm = th
                rays.append(mray(np.arcsin(self.NA) * thm, phi))
                
        return rays
    
    def exampleRays(self,d, nph = 9, nth=2, jit=True):
        self.extend(self._generate_rays(nph, nth, jit))


class FanSource(PointSource):
    def exampleRays(self, d, phi=0, nth=2, jit=True):
        y1 = self.NA
        z1 = pl.sqrt(1 - y1 ** 2)
        y2 = y1 / pl.sqrt(2)
        z2 = pl.sqrt(1 - y2 ** 2)
        d2 = self.placement.orientation
        #d1 = pl.cross(d2, pl.array([1,0,0]))
        d1 = pl.cross(d2, pl.array([1, 1, 1]) - d2)
        d1 = d1 / pl.norm(d1)
        d0 = pl.cross(d1, d2)
        d0 = d0 / pl.norm(d0)
        
        def mray(theta, phi):
            r = np.exp(1j * phi)
            #print r.real, r.imag, np.cos(theta)
            dn = np.sin(theta) * r.real * d0 + np.sin(theta) * r.imag * d1 + np.cos(theta) * d2
            #print pl.norm(dn)
            #dn = dn/pl.norm(dn)
            
            return Ray(self.placement.location, dn, self.material, color=self.color, wavelength=self.wavelength)
        
        r0 = Ray(self.placement.location, d2, self.material, color=self.color, wavelength=self.wavelength)
        #        ryp = Ray(self.placement.location,y1*d1 + z1*d2, self.material, color=self.color)
        #        rypp = Ray(self.placement.location,y2*d1 + z2*d2, self.material, color=self.color)
        #        ryn = Ray(self.placement.location,-y1*d1 + z1*d2, self.material, color=self.color)
        #        rynn = Ray(self.placement.location,-y2*d1 + z2*d2, self.material, color=self.color)
        #        rxp = Ray(self.placement.location,y1*d0 + z1*d2, self.material, color=self.color)
        #        rxn = Ray(self.placement.location,-y1*d0 + z1*d2, self.material, color=self.color)
        self.append(r0)
        #        self.append(ryp)
        #        self.append(rypp)
        #        self.append(rynn)
        #        self.append(ryn)
        #        self.append(rxp)
        #        self.append(rxn)
        for th in np.linspace(-1, 1, nth):
            #for phi in (np.linspace(0, 2 * np.pi, np.maximum(nph * th, 3))[:-1] + 3 * np.pi / 4):
            if jit:
                phi = phi + (np.random.rand(1) - .5) * 2 * np.pi / (nph * th)
                thm = th + (np.random.rand(1) - .5) / nth
            else:
                thm = th
            self.append(mray(np.arcsin(self.NA) * thm, phi))
            
            #for phi in np.linspace(0, 2*np.pi, nph):
            #    self.append(mray(np.arcsin(self.NA)/2, phi))
        

def SourceTest() :
    s = Source("test",Placement([0,0,0],[0,0,1]))
    s.exampleRays(0.04)
    return s

import numpy as np

from .elements import *
from .rays import *
from .placement import *
from .material import *

class Source(Element, list) :
    def __init__(self,name,placement, wavelength=635.) :
        
        self.name = name
        self.placement = placement
        self.dimension = np.array([0.05,0.05,0.01])
        self.material = Material(Material.REFRACT, 1.0)
        self.wavelength=wavelength

    def __str__(self) :
        s  = 'Source                   : '+self.name+'\n'
        s += 'Source.placement         : \n'+str(self.placement)
        return s

    def surface(self) :
        x = np.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = np.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = np.meshgrid(x,y)
        zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        
        return [xx,yy,zz]        
    
    def example_rays(self, d) :
        r0  = RayBundle(self.placement.location, [0, 0, 1], self.material)
        ryp = RayBundle(self.placement.location + np.array([0, d, 0]), [0, 0, 1], self.material)
        rypp = RayBundle(self.placement.location + np.array([0, d / 2, 0]), [0, 0, 1], self.material)
        ryn = RayBundle(self.placement.location + np.array([0, -d, 0]), [0, 0, 1], self.material)
        rynn = RayBundle(self.placement.location + np.array([0, -d / 2, 0]), [0, 0, 1], self.material)
        rxp = RayBundle(self.placement.location + np.array([d, 0, 0]), [0, 0, 1], self.material)
        rxn = RayBundle(self.placement.location + np.array([-d, 0, 0]), [0, 0, 1], self.material)
        self.append(r0)
        self.append(ryp)
        self.append(rypp)
        self.append(rynn)
        self.append(ryn)
        self.append(rxp)
        self.append(rxn)



class PointSource(Source) :
    def __init__(self,name,placement, NA=1.49/1.51, color=(1.0, 0, 0), wavelength=635., **kwargs) :
        #print 'Source.__init__'
        self.name = name
        self.placement = placement
        self.dimension = np.array([0.05,0.05,0.01])
        self.material = Material(Material.REFRACT, 1.0)
        self.NA = NA
        self.color = color
        self.wavelength=wavelength

    def __str__(self) :
        s  = 'Source                   : '+self.name+'\n'
        s += 'Source.placement         : \n'+str(self.placement)
        return s
    
    def copy(self, **kwargs):
        args = dict(self.__dict__)
        
        args.update(kwargs)
        
        return self.__class__(**args)

    def surface(self) :
        x = np.arange(-self.dimension[0],self.dimension[0]+1e-8,self.dimension[0]/5)
        y = np.arange(-self.dimension[1],self.dimension[1]+1e-8,self.dimension[1]/5)
        xx,yy = np.meshgrid(x,y)
        zz = 1/self.placement.orientation[2]*(np.linalg.dot(self.placement.orientation,self.placement.location)-self.placement.orientation[0]*xx-self.placement.orientation[1]*yy)
        
        return [xx,yy,zz]
    
    
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
    
    @property
    def principle_ray(self):
        return RayBundle(self.placement.location, self.placement.orientation, self.material, color=self.color, wavelength=self.wavelength)
    
    def _generate_rays(self, nph, nth, jit=False):
        d2 = self.placement.orientation
        #d1 = np.cross(d2, np.array([1,0,0]))
        d1 = np.cross(d2, np.array([1, 1, 1]) - d2)
        d1 = d1 / np.linalg.norm(d1)
        d0 = np.cross(d1, d2)
        d0 = d0 / np.linalg.norm(d0)
        
        rays = []
    
        def mray(theta, phi):
            r = np.exp(1j * phi)

            dn = (np.sin(theta) * r.real) * d0 + (np.sin(theta) * r.imag) * d1 + np.cos(theta) * d2
        
            return RayBundle(np.zeros_like(dn) + self.placement.location, dn, self.material, color=self.color, wavelength=self.wavelength)
    
        
        ths = [0]
        phis = [0]
        
        
        for th in np.linspace(0, 1, nth)[1:]:
            for phi in (np.linspace(0, 2 * np.pi, np.maximum(nph * th, 3))[:-1] + 3 * np.pi / 4):
                if jit:
                    phi = phi + (np.random.rand(1) - .5) * 2 * np.pi / (nph * th)
                    thm = th + (np.random.rand(1) - .5) / nth
                else:
                    thm = th
                    
                ths.append(np.arcsin(self.NA) * thm)
                phis.append(phi)
                
        
        rays = [mray(np.array(ths)[:,None], np.array(phis)[:,None]), ]
                
        return rays
    
    def example_rays(self, d, nph = 9, nth=2, jit=True):
        self.extend(self._generate_rays(nph, nth, jit))
        
class ColimatedSource(PointSource):
    def __init__(self, name, placement, diameter=3.0, color=(1.0,0,0), wavelength=635.):
        PointSource.__init__(self, name, placement, NA=0, color=color, wavelength=wavelength)
        
        self.diameter = diameter

    def _generate_rays(self, nph, nth, jit=False):
        d2 = self.placement.orientation
        #d1 = np.cross(d2, np.array([1,0,0]))
        d1 = np.cross(d2, np.array([1, 1, 1]) - d2)
        d1 = d1 / np.linalg.norm(d1)
        d0 = np.cross(d1, d2)
        d0 = d0 / np.linalg.norm(d0)
    
        rays = []
    
        def mray(theta, phi):
            r = np.exp(1j * phi)

            dn = (theta * r.real) * d0 + (theta * r.imag) * d1

        
            return RayBundle(dn + self.placement.location, d2*np.ones_like(dn), self.material, color=self.color,
                             wavelength=self.wavelength)
    
    
        ths = [0]
        phis = [0]
    
        for th in np.linspace(0, 1, nth)[1:]:
            for phi in (np.linspace(0, 2 * np.pi, np.maximum(nph * th, 3))[:-1] + 3 * np.pi / 4):
                if jit:
                    phi = phi + (np.random.rand(1) - .5) * 2 * np.pi / (nph * th)
                    thm = th + (np.random.rand(1) - .5) / nth
                else:
                    thm = th
            
                ths.append(0.5*self.diameter * thm)
                phis.append(phi)
    
        #rays = [mray(t, p) for t, p in zip(ths, phis)]
    
        rays = [mray(np.array(ths)[:, None], np.array(phis)[:, None]), ]
    
        return rays


class FanSource(PointSource):
    #FIXME - does this still work?
    def example_rays(self, d, phi=0, nth=2, jit=True):
        y1 = self.NA
        z1 = np.sqrt(1 - y1 ** 2)
        y2 = y1 / np.sqrt(2)
        z2 = np.sqrt(1 - y2 ** 2)
        d2 = self.placement.orientation
        #d1 = np.cross(d2, np.array([1,0,0]))
        d1 = np.cross(d2, np.array([1, 1, 1]) - d2)
        d1 = d1 / np.linalg.norm(d1)
        d0 = np.cross(d1, d2)
        d0 = d0 / np.linalg.norm(d0)
        
        def mray(theta, phi):
            r = np.exp(1j * phi)
            #print r.real, r.imag, np.cos(theta)
            dn = np.sin(theta) * r.real * d0 + np.sin(theta) * r.imag * d1 + np.cos(theta) * d2
            #print np.linalg.norm(dn)
            #dn = dn/np.linalg.norm(dn)
            
            return RayBundle(self.placement.location, dn, self.material, color=self.color, wavelength=self.wavelength)
        
        r0 = RayBundle(self.placement.location, d2, self.material, color=self.color, wavelength=self.wavelength)
        
        self.append(r0)

        for th in np.linspace(-1, 1, nth):
            #for phi in (np.linspace(0, 2 * np.pi, np.maximum(nph * th, 3))[:-1] + 3 * np.pi / 4):
            if jit:
                phi = phi + (np.random.rand(1) - .5) * 2 * np.pi / (nph * th)
                thm = th + (np.random.rand(1) - .5) / nth
            else:
                thm = th
            self.append(mray(np.arcsin(self.NA) * thm, phi))
            


def SourceTest() :
    s = Source("test",Placement([0,0,0],[0,0,1]))
    s.example_rays(0.04)
    return s

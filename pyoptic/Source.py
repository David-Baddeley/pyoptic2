import pylab as pl

from Elements import *
from Ray import *
from Placement import *

class Source(Volume,list) :
    def __init__(self,name,placement) :
        print 'Source.__init__'
        self.name = name
        self.placement = placement
        self.dimension = pl.array([0.05,0.05,0.01])
        self.material = Material(Material.refract,1.0)

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
        

def SourceTest() :
    s = Source("test",Placement([0,0,0],[0,0,1]))
    s.exampleRays(0.04)
    return s

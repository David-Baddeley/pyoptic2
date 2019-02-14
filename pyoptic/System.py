import pylab as py

from Placement import *
from  Source import *
from Material  import *
from Elements  import *
#from Display3D import *

# Complete optical system
class System(list) : 

    def __init__(self) :
        print "System:__init__>"

        # maximum size as determined from the elements
        self.size = [0,0,0]

    def loadSystem(self) :
        print "System:laodSystem>"

    def checkSystem(self) :
        print "System:checkSystem>"
        
    def propagate(self, source=None) :
        #print "System:rayTracing>"

        raytree = []
        
        if source is None:
            source = self[0]
            
        if isinstance(self[0],Source):
            startAt = 1
        else:
            startAt = 0

        # iterate over rays from source
        i = 0 
        ri = iter(source)
        try :                
            while True :
                r = ri.next()       
                #print 'System.proagate> ray='+str(i)
#                raybranch = []
#                raybranch.append(r)
#                for j in range(1,len(self)) :
#                    #print 'System.propagate> element=',j
#                    rp = self[j].propagate(self[j-1],raybranch[j-1])                    
#                    raybranch.append(rp)
                raybranch = self.propagateRay(r, startAt=startAt, source=source)
                raytree.append(raybranch)
                i += 1
        except StopIteration :
            pass
            

        return raytree
        
    def propagateRay(self, r, startAt=1, source=None):
        raybranch = []
        raybranch.append(r)
        for j in range(startAt,len(self)) :
            #print 'System.propagate> element=',j
            if j == 0:
                rp = self[j].propagate(source, raybranch[j-startAt])
            else:
                rp = self[j].propagate(self[j-1],raybranch[j-startAt])
                
            if rp == None:
                break                   
            raybranch.append(rp)
            
        return raybranch
    
    def prepend(self, elements):
        try:
            for e in elements[::-1]:
                self.insert(0, e)
        except AttributeError:
            #single element
            self.insert(0, elements)
            
    def add(self, elements):
        try:
            self.extend(elements)
        except TypeError:
            self.append(elements)
                

                
    def __str__(self) :
        s = ''
        for e in self :
            s +=str(e)
        return s

def SystemTest() :
    print "System:SystemTest>"

    # source 
    ppos = [0,0,0]
    pdir = [0,0,1]
    s0 = Source("light source",Placement(ppos,pdir))
    s0.exampleRays(0.04)

    # curved surface 
    ppos = [0,0,0.05]
    pdir = [0,0,1]
    pdim = [0.05,0.05,0.01]
    s1 = SphericalSurface("spherical 1", Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.4),0.08)
    
    # plane surface
    ppos = [0,0,0.07]
    pdir = [0,0,1]
    pdim = [0.05,0.05,0.01]
    s2 = PlaneSurface("plane 1",Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.0))
    
    # plane surface
    ppos = [0,0,0.40]
    pdir = [0,0,1]
    pdim = [0.05,0.05,0.01]
    s3 = PlaneSurface("plane 2",Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.4))

    # curved surface
    ppos = [0,0,0.42]
    pdir = [0,0,1]
    pdim = [0.05,0.05,0.01]
    s4 = SphericalSurface("spherical 2", Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.0),-0.08)
    
    # plane surface
    ppos = [0,0,0.50]
    pdir = pl.array([0,-1,1])
    pdir = pdir/pl.linalg.norm(pdir)
    pdim = [0.05,0.05,0.01]
    s5 = PlaneSurface("plane 3",Volume.circ,pdim,Placement(ppos,pdir),Material(Material.mirror,1.0))

    # plane surface
    ppos = [0,0.15,0.50]
    pdir = [0,1,0]
    pdim = [0.05,0.05,0.01]
    s6 = PlaneSurface("plane 4",Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.0))

    # system
    s = System()
    s.append(s0)
    s.append(s1)
    s.append(s2)
    s.append(s3)
    s.append(s4)
    s.append(s5)
    s.append(s6)

    return s

    # ray trace through system
#    rl = s.propagate()
#
#    while r != None :
#        rl.append(r)
#        r1 = s[1].propagate(s[0],r)
#        print r
#        print r1[0]
#        r = s[0].nextRay()
#
#    d = Display3D(s,rl)
#    d.Draw()

    return s


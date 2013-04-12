from System import *
from Elements import *
from Source import *

from enthought.tvtk.tools import mlab
#from mayavi import mlab
from enthought.tvtk.tools import visual

class Display3D :
    def __init__(self,s,r) :
        self.s = s
        self.r = r
        self.f = mlab.figure()

        self.e3d = []

    def Draw(self) :        
        # loop over optical elements
        #self.f.scene.disable_render = True
        for e in self.s :
            x,y,z = e.surface()
            #print e.name
            #print z
            edo = mlab.Surf(x,y,z)
            edo.trait_set(representation='wireframe')
            self.f.add(edo)
            self.e3d.append(edo)
            
        # loop over rays
        #print self.r
#        for r in pl.flatten(self.r) :
#            if r.p1 != None :
#                print 'line'
#                print r.p0
#                print r.p1                
#                # ray display object
#
#                rdo = mlab.Line3([r.p0,r.p1],radius=0.1,representation='wireframe')
##                rdo.representation = 'wireframe'
#                self.f.add(rdo)
#                self.e3d.append(rdo)
        for rb in self.r:
            rbp = [r.p0 for r in rb] 
            if (not rb[-1].p1 == None):
                rbp.append(rb[-1].p1)
                
            if (len(rbp) > 1):
                rdo = mlab.Line3(rbp,radius=0.05,representation='wireframe')
    #           rdo.representation = 'wireframe'
                self.f.add(rdo)
                self.e3d.append(rdo)
                
        #self.f.scene.disable_render = False

def Display3DTest() :
    s = SystemTest()
    r = s.propagate()
    d = Display3D(s,r)
    d.Draw()
    return d

    

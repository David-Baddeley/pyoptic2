from System import *
from Elements import *
#from Source import *

import numpy as np

try:
    from enthought.tvtk.tools import mlab as ml2
    from mayavi import mlab
    from enthought.tvtk.tools import visual
    from enthought.tvtk.api import tvtk
    from enthought.tvtk.common import configure_input_data
except ImportError:
    from tvtk.tools import mlab as ml2
    from tvtk.tools import visual
    from tvtk.api import tvtk
    from tvtk.common import configure_input_data
    from mayavi import mlab

class Display3D :
    def __init__(self,s,r) :
        self.s = s
        self.r = r
        self.f = mlab.figure()
        
        #if 

        self.e3d = []
        
    def _drawLine(self, points, color=[1,1,1]):
        npts = len(points) - 1
        lines = np.zeros((npts, 2), 'l')
        lines[:,0] = np.arange(0, npts-0.5, 1, 'l')
        lines[:,1] = np.arange(1, npts+0.5, 1, 'l')
        d = tvtk.PolyData(points=points, lines=lines)
        m = tvtk.PolyDataMapper()
        configure_input_data(m, d)
        a = tvtk.Actor(mapper=m)
        a.property.color = color
        
        self.f.scene.add_actor(a)
        return a

    def Draw(self) :        
        # loop over optical elements
        self.f.scene.disable_render = True
        for e in self.s :
            x,y,z = e.surface()
            #print e.name
            #print z
            edo = mlab.mesh(x,y,z, color=(.9,1,1), opacity=0.7)
            #edo.trait_set(representation='wireframe')
            #edo.use_tubes = False
            #self.f.add(edo)
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
                #print rbp
                #rdo = ml2.Line3(rbp,radius=0.05, color=r.color)
    #           rdo.representation = 'wireframe'
                #rdo.use_tubes = False
                #self.f.scene.add(rdo)
                rdo = self._drawLine(rbp, color=r.color)
                self.e3d.append(rdo)
                
        self.f.scene.background = (0.38, 0.4, 0.55)
        self.f.scene.disable_render = False

def Display3DTest() :
    s = SystemTest()
    r = s.propagate()
    d = Display3D(s,r)
    d.Draw()
    return d

    

from .system import *
from .elements import *
#from Source import *

import numpy as np

try:
    #from enthought.tvtk.tools import mlab as ml2
    from mayavi import mlab
    #from enthought.tvtk.tools import visual
    from enthought.tvtk.api import tvtk
    from enthought.tvtk.common import configure_input_data
except ImportError:
    #from tvtk.tools import mlab as ml2
    #from tvtk.tools import visual
    from tvtk.api import tvtk
    from tvtk.common import configure_input_data
    from mayavi import mlab


class STLPart(object):
    def __init__(self, filename, colour=None, position=[0, 0, 0], orientation=None, rotation=None):
        self.filename = filename
        self.reader = tvtk.STLReader(file_name=filename)
        self.mapper = tvtk.PolyDataMapper()
        configure_input_data(self.mapper, self.reader.output)
        self.reader.update()
        self.actor = tvtk.Actor(mapper=self.mapper)
        
        if colour:
            self.actor.property.color = colour
        
        self.actor.position = position
        if rotation is None:
            if orientation is None:
                self.actor.orientation = [0,0,0]
            else:
                a2 = -180*np.arctan(orientation[0]/orientation[1])/np.pi
                a1 = 180*np.arcsin(orientation[2])/np.pi
                self.actor.orientation = [a1, 0, a2]
        else:
            self.actor.orientation = rotation
    
    def as_scad(self, origin=[0, 0, 0]):
        trans = '[%f, %f, %f]' % tuple(np.array(self.actor.position) - np.array(origin))
        rot = '[%f,%f,%f]' % tuple(self.actor.orientation)
        return 'translate(%s)rotate(%s)import("%s");\n' % (trans, rot, self.filename)
    
    def add_to_scene(self, scene):
        scene.add_actor(self.actor)
    
    @property
    def position(self):
        return self.actor.postion
    
    @position.setter
    def set_pos(self, value):
        self.actor.position = value
    
    @property
    def orientation(self):
        return self.actor.orientation
    
    @orientation.setter
    def set_or(self, value):
        self.actor.orientation = value

class Display3D :
    def __init__(self,s,r, stl_parts=[]):
        self.s = s
        self.r = r
        self.f = mlab.figure()
        
        self.stl_parts = stl_parts
        
        #if 

        self.e3d = []
        
    def _drawLines(self, points, color=[1,1,1]):
        nrays = points.shape[0]
        npts = points.shape[2]
        nsegs = npts - 1
        lines = np.zeros((nsegs*nrays, 2), 'l')
        pts_cat = np.zeros([npts*nrays, 3], 'f')
        for j in range(nrays):
            pts_cat[(j*npts):((j+1)*npts), :] = points[j,:,:].T
            lines[(j*nsegs):((j+1)*nsegs),0] = (j*npts) + np.arange(0, nsegs-0.5, 1, 'l')
            lines[(j * nsegs):((j + 1) * nsegs),1] = (j * npts) + np.arange(1, nsegs + 0.5, 1, 'l')
            #lines[:,1] = np.arange(1, nsegs+0.5, 1, 'l')
            
        d = tvtk.PolyData(points=pts_cat, lines=lines)
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
            ra = np.concatenate([np.atleast_3d(ri.p0) for ri in rb], 2)
            #ra = np.array([ri.p0 for ri in ray])
            if not rb[-1].p1 == None:
                ra = np.concatenate([ra, np.atleast_3d(rb[-1].p1)], 2)
            
            # rbp = [r.p0 for r in rb]
            # if (not rb[-1].p1 == None):
            #     rbp.append(rb[-1].p1)
                
            if (ra.shape[0] > 1):
                #print rbp
                #rdo = ml2.Line3(rbp,radius=0.05, color=r.color)
    #           rdo.representation = 'wireframe'
                #rdo.use_tubes = False
                #self.f.scene.add(rdo)
                rdo = self._drawLines(ra, color=rb[0].color)
                self.e3d.append(rdo)
                
        for p in self.stl_parts:
            p.add_to_scene(self.f.scene)
                
        self.f.scene.background = (0.38, 0.4, 0.55)
        self.f.scene.disable_render = False

def Display3DTest() :
    s = SystemTest()
    r = s.propagate()
    d = Display3D(s,r)
    d.Draw()
    return d

    

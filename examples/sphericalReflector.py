from System   import *
from Elements import *
from Source   import *
import Display3D

def sphericalReflector() :    
    # source 
    ppos = [0,0,0]
    pdir = [0,0,1]
    s0 = Source("light source",Placement(ppos,pdir))
    s0.example_rays(0.04)

    # spherical reflector 
    ppos = [0,0,0.20]
    pdir = [0,0,1]    
    pdim = [0.05,0.05,0.01]
    sr = SphericalSurface("spherical", Volume.circ,pdim,Placement(ppos,pdir),Material(Material.mirror,1.0),-0.10)
    
    # plane stop
    ppos = [0,0,0.1]
    pdir = [0,0,-1]    
    pdim = [0.05,0.05,0.01]                          
    ss = PlaneSurface("stop",Volume.circ,pdim,Placement(ppos,pdir),Material(Material.refract,1.0))    

    s = System()
    s.append(s0)
    s.append(sr)
    s.append(ss)

    r = s.propagate()
    d = Display3D.Display3D(s,r)
    d.Draw()    

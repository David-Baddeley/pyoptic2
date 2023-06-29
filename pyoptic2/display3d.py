import numpy as np
import pythreejs as p3j


class Display3D(object):
    def __init__(self,s,r, parts=[]):
        self.s = s
        self.r = r
        
        self._parts = parts
        
        #if 

        self.e3d = []
        
    def _drawLines(self, points, color=[1,1,1]):
        nrays = points.shape[0]
        
        color = f'rgb({int(255*color[0])}, {int(255*color[1])}, {int(255*color[2])})'
        o = []
        
        for j in range(nrays):
            o.append(p3j.Line(geometry=p3j.BufferGeometry(attributes={'position': p3j.BufferAttribute(points[j,:,:].T)}),
                              material=p3j.LineBasicMaterial(color=color)))
            
        return o

    def _material(self, color, transparent=False, opacity=1.0, specular='#afafff', shininess=30):
        material = p3j.MeshPhongMaterial()
        #material = p3.CustomMaterial("standard")
        material.color = color
        material.clipping = True
        material.side = "DoubleSide"
        material.polygonOffset = True
        material.polygonOffsetFactor = 1
        material.polygonOffsetUnits = 1
        material.transparent = transparent
        material.opacity = opacity
        material.specular = specular
        material.shininess = shininess
        material.alpha = opacity
        #material.update("metalness", 0.3)
        #material.update("roughness", 0.8)
        return material
    
    def _calc_normals(self, verts):
        v0 = verts[0::3, :] - verts[1::3, :] + np.array([0.1, 0, 0])
        v1 = verts[0::3, :] - verts[2::3, :] + np.array([0, 0.1, 0])
        n = np.cross(v0, v1)
        n = n / np.sqrt(np.sum(n**2, axis=1))[:, np.newaxis]

        n = np.repeat(n.astype('f4'), 3, axis=0)

        return n
    
    def Draw(self, shape=(1000, 600)) :        
        # loop over optical elements

        elements = []

        bbox = np.array([-1, 1, -1, 1, -1, 1]) # x0, x1, y0, y1, z0, z1

        for e in self.s :
            verts = e.surface_triangles()
            #print(verts.shape)

            if (verts.shape[0] == 0):
                continue

            bbox[0] = min(bbox[0], np.min(verts[:,0]))
            bbox[1] = max(bbox[1], np.max(verts[:,0]))
            bbox[2] = min(bbox[2], np.min(verts[:,1]))
            bbox[3] = max(bbox[3], np.max(verts[:,1]))
            bbox[4] = min(bbox[4], np.min(verts[:,2]))
            bbox[5] = max(bbox[5], np.max(verts[:,2]))

            
            if e.material.type == e.material.REFLECT:
                m = self._material('silver', specular='#efefff', shininess=90)
            else:
                m = self._material('lightblue', opacity=0.8, shininess=20)
                                                      
            
            mesh = p3j.Mesh(geometry=p3j.BufferGeometry(attributes={'position': p3j.BufferAttribute(verts), 
                                                                   'index': p3j.BufferAttribute(np.arange(verts.shape[0], dtype='uint32')),
                                                                   'normal': p3j.BufferAttribute(self._calc_normals(verts))
                                                                   }),
                            material=m)

            elements.append(mesh)

        #print(elements)

        # loop over rays
        for rb in self.r:
            ra = np.concatenate([np.atleast_3d(ri.p0) for ri in rb], 2)
            if not rb[-1].p1 == None:
                ra = np.concatenate([ra, np.atleast_3d(rb[-1].p1)], 2)
                
            if (ra.shape[0] > 1):
                bbox[0] = min(bbox[0], np.min(ra[:,:,0]))
                bbox[1] = max(bbox[1], np.max(ra[:,:,0]))
                bbox[2] = min(bbox[2], np.min(ra[:,:,1]))
                bbox[3] = max(bbox[3], np.max(ra[:,:,1]))
                bbox[4] = min(bbox[4], np.min(ra[:,:,2]))
                bbox[5] = max(bbox[5], np.max(ra[:,:,2]))

                rdo = self._drawLines(ra, color=rb[0].color)
                elements.extend(rdo)

        
        # loop over parts (things which are responsible for drawing themselves, typically optomechanics with a STL or STEP description)
        # these should either already be a pythreejs object, or have a to_threejs() method which returns one
        for p in self._parts:
            if hasattr(p, 'to_threejs'):
                elements.append(p.to_threejs())
            else:
                elements.append(p)
        
        # create camera, light, and scene
        key_light = p3j.DirectionalLight(color='white', position=[3, 5, 1], intensity=0.5)

        box_centre = (bbox[0::2] + bbox[1::2]) / 2

        c = p3j.CombinedCamera(position=tuple(box_centre + np.array([0, 0, 100])), up=[1, 0, 0], children=[key_light], 
                            width=shape[0], height=shape[1], mode='orthographic')
        c.up=[1, 0, 0]
        c.right=[0, 1, 0]
        c.position = tuple(box_centre + np.array([0, 0, 100]))
        
        c.target = tuple(box_centre)

        
        s = p3j.Scene(children=elements + [c, p3j.AmbientLight(color='#777777')], background='#eeeeee')
        renderer = p3j.Renderer(camera=c,
                            scene=s,
                            alpha=True,
                            clearOpacity=0,
                            controls=[p3j.OrbitControls(controlling=c)],
                            width=shape[0], height=shape[1])
        
        return renderer
                        
        

def Display3DTest() :
    from pyoptic2.system import SystemTest
    s = SystemTest()
    r = s.propagate()
    d = Display3D(s,r)
    d.Draw()
    return d

    

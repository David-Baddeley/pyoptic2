"""
Classes defining optomechanical components (represented by a 3D bounding box).
"""

import numpy as np

class OpticHolder(object):
    _AXES = {'x': np.array([1, 0, 0]), 'y': np.array([0, 1, 0]), 'z': np.array([0, 0, 1])}
    bbox = [-20, 20, 0, 10, -20, 20]
    hole_radius = 2.5
    
    _face_triangs = [
        [0, 2, 4],
        [4, 2, 6],
        [5, 4, 6],
        [7, 6, 5],
        [2, 3, 6],
        [6, 3, 7],
        [4, 1, 0],
        [4, 5, 1],
        [0, 3, 2],
        [1, 0, 3],
        [6, 7, 3],
        [6, 5, 7],
        [0, 4, 1],
        [4, 5, 1],
    ]
    
    def __init__(self, optic, up=[0, 0, 1], bbox=None, offset=0, stl=None, bbox_optic=None,
                 bbox_cutout=None, hole_pos=None, hole_radius=None, **kwargs):
        """ An Optomechanical holder for an optical element
        
        parameters:
        ===========
        
        optic : the element which the optomechanical part holds, defines the position and orientation of the part.
                Can also be a placement.Placement object
        up : the parts up direction (a 3-tuple)
        bbox : the bounding box used for display purposes (and as a fallback for generating cutouts). Bounding boxes are
               [xmin, xmax, ymin, ymax, zmin, zmax] and should be specified relative to the position of the optical
               element (i.e. the element thould be at (0,0,0) within the bounding box. The co-ordinate system is such that
               z is along the "up direction" y is along the element orientation (usually the direction of propagation), and
               x is orthogonal to the previous 2 directions.
        offset : a 3-tuple or a single scalar representing an offset of the optic from the (0,0,0) position of the bounding box
                 used where the optic could be referenced to different positions within the holder (e.g. lens tube or cage plates)
        stl : Not implemented (will hopefully allow bounding box to be derived automatically from an STL of the part at some
              point in the future.
        bbox_optic : a bounding box for the optic. Currently used for drawing mirrors on mirror holders where the
                     representation of the element is simply a flat plane, but the actual physical optic has thickness
        bbox_cutout : the bounding box to use when generating a cutout (see cutout method below), if different from the
                      overall bounding box
        hole_position : the position of the mounting hole (used for cutouts)
        hole_radius : the clearance radius for the mounting hole
        
        """
        self._optic = optic
        
        self._stl = stl
        if stl is not None:
            bbox = self._stl_bbox(stl)
        elif bbox is None:
            bbox = self.bbox
        
        if np.isscalar(offset):
            offset = np.array([0, offset, 0])
        
        self._bbox = self._offset_bbox(bbox, offset)
        
        self._bbox_optic = self._offset_bbox(
            bbox_optic if bbox_optic is not None else getattr(self, 'bbox_optic', None), offset)
        self._bbox_cutout = self._offset_bbox(
            bbox_cutout if bbox_cutout is not None else getattr(self, 'bbox_cutout', None), offset)
        
        self._hole_pos = hole_pos if hole_pos is not None else getattr(self, 'hole_pos', None)
        if self._hole_pos is not None:
            self._hole_pos = np.array(self._hole_pos) + offset[:2]
        
        self._hole_radius = hole_radius if hole_radius else self.hole_radius
        
        self._up = np.array(up)
        self._front = np.array(kwargs.get('orientation', optic.orientation))
        self._right = np.cross(up, self._front)
        self._right = self._right / np.linalg.norm(self._right)
        
        self._rmat = np.array([self._right, self._front, self._up])
    
    def _offset_bbox(self, bbox, offset):
        if bbox is None:
            return None
        
        x_min, x_max, y_min, y_max, z_min, z_max = bbox
        if len(offset) == 3:
            return np.array([x_min + offset[0], x_max + offset[0],
                             y_min + offset[1], y_max + offset[1],
                             z_min + offset[2], z_max + offset[2]])
        
        else:
            raise RuntimeError('Offset should be a scalar or a 3-tuple')
    
    def _bbox_corners(self, bbox):
        x_min, x_max, y_min, y_max, z_min, z_max = bbox
        return np.array([[x_min, y_min, z_min],
                         [x_min, y_max, z_min],
                         [x_min, y_min, z_max],
                         [x_min, y_max, z_max],
                         [x_max, y_min, z_min],
                         [x_max, y_max, z_min],
                         [x_max, y_min, z_max],
                         [x_max, y_max, z_max],
                         ])
    
    def _corners(self, bbox):
        return np.dot(self._bbox_corners(bbox), self._rmat) + self._optic.location[None, :]
    
    @property
    def corners(self):
        return self._corners(self._bbox)
    
    def _project(self, points, proj='xy'):
        x = (points * self._AXES[proj[0]][None, :]).sum(1)
        y = (points * self._AXES[proj[1]][None, :]).sum(1)
        
        return x, y
    
    def _projected_footprint(self, points, proj='xy'):
        from scipy.spatial import ConvexHull
        
        x, y = self._project(points, proj)
        
        pts = np.array([x, y]).T
        hull = ConvexHull(pts)
        vertices = np.hstack([hull.vertices, hull.vertices[0]])
        
        return pts[vertices, :].T
    
    def projected_footprint(self, proj='xy'):
        return self._projected_footprint(self.corners, proj=proj)
    
    def plot(self, proj='xy'):
        import matplotlib.pyplot as plt
        x, y = self.projected_footprint(proj=proj)
        plt.plot(y, x, 'k--')
        
        if self._bbox_optic is not None:
            #support drawing mirrors
            x, y = self._projected_footprint(self._corners(self._bbox_optic), proj=proj)
            plt.plot(y, x, 'k-')
        
        if self._bbox_cutout is not None:
            x, y = self._projected_footprint(self._corners(self._bbox_cutout), proj=proj)
            plt.plot(y, x, 'r:')
        
        if self._hole_pos is not None:
            hp = np.dot(np.hstack([self._hole_pos, 0])[None, :], self._rmat) + self._optic.location[None, :]
            hx, hy = self._project(hp, proj=proj)
            plt.plot(hy, hx, '+r')
    
    def plot3d(self, color=(0.7, 0.5, 0.5), opacity=0.4):
        from mayavi import mlab
        
        x, y, z = self.corners.T
        mlab.triangular_mesh(x, y, z, self._face_triangs, opacity=opacity, color=color)
    
    def cutout(self, plane='xy', radius=2):
        """
        Generate OpenSCAD code for creating cutouts in a plate to be used for positioning an optical element.
        
        These cutouts consist of the a rectangular cutout, with clearance overcuts in the corners, and a through hole
        for mounting
    
        relies on the following OpenSCAD module definition:
        
        module cutout(corners, radii, z_offset, hole_pos, hole_radius, corner_dir=[1,0]){
          union(){
            translate([0,0,z_offset])linear_extrude(50, center=false)offset(-radii)offset(radii)union(){
              polygon(corners);
              offs = [1,1,-1,-1,1];
              
                for (i = [0:(len(corners) - 1)]) {
                  p = corners[i];
                    //echo(sign(corners[i+1][0] - p[0]));
                  pc = p + radii*corner_dir*offs[i];
                  translate([pc[0], pc[1], 0]) circle(radii);};
              };
          
            if (hole_pos != false) translate([hole_pos[0], hole_pos[1], 0]) cylinder(r=hole_radius, h=100, center=true);
          };
        };
        
        """
        if self._hole_pos is None:
            hole_pos = 'false'
        else:
            hp = np.dot(np.hstack([self._hole_pos, 0])[None, :], self._rmat) + self._optic.location[None, :]
            hx, hy = self._project(hp, proj=plane)
            #print(hx, hy)
            hole_pos = '[%3.2f, %3.2f]' % (float(hx), float(hy))
        
        cut_bbox = self._bbox_cutout if self._bbox_cutout is not None else self._bbox #fall back on bbox if we haven't explicitly specified a cutting bbox
        #TODO - allow this to be non-rectangular
        
        pts = self._projected_footprint(self._corners(cut_bbox), proj=plane)
        corners = '[' + ', '.join(['[%3.2f, %3.2f]' % (float(p[0]), float(p[1])) for p in pts.T]) + ']'
        
        z_offset = cut_bbox[4]
        
        midp = pts.mean(1)
        mid_dist = pts - midp[:, None]
        
        right = np.array(self._project(self._right))
        right = right / np.linalg.norm(right)
        #print len(points.T)
        #FIXME - this isn't working correctly
        if np.mean((mid_dist - right[:, None] * np.array([1, 1, -1, -1, 1])[None, :len(pts.T)]) ** 2) < np.mean(
                        mid_dist ** 2):
            right = -right
        
        return 'cutout(%s, %3.2f, %3.2f, %s, %3.2f, [%3.2f, %3.2f]);\n' % (
        corners, radius, z_offset, hole_pos, self._hole_radius, right[0], right[1])


############################
# A few Thorlabs parts

class CP02T(OpticHolder):
    bbox = [-20.3, 20.3, 0, 12.7, -20.3, 20.3]
    bbox_cutout = [-20.8, 20.8, 0, 13.7, -20.3, 20.3]
    hole_pos = [0, 6.4]
    offsets = {'centre': 12.7 / 2, 'front': 2, 'back': 12.7 - 2}
    
    def __init__(self, lens, align='centre', **kwargs):
        offset = kwargs.pop('offset', 0) - self.offsets[align]
        OpticHolder.__init__(self, lens, offset=offset, **kwargs)


class KMSS(OpticHolder):
    bbox = [-12.7, 12.7, 6, 26, -15.9, 12.7]
    bbox_optic = [-12.5, 12.5, 0, 6, -12.5, 12.5]
    bbox_cutout = [-12.7, 12.7, 14.4, 22.4, -15.9, 12.7]
    hole_pos = [0, 18.8]


class GVS211(OpticHolder):
    bbox = [-14, 32, -31, 15, -16.2, 69.8]
    bbox_cutout = [-14, 33, -32, 15, -16.2, 69.8]
    hole_pos = [-(15 - 23), -(-14 + 23)]


class VA100(OpticHolder):
    bbox = [-25.4, 25.4, -6.9, 5, -17.8, 16]
    slit_bb = [-15, -1.5, -.05, .05, -5, 5]
    hole_pos = [0, 3 - (9.9 - 4.8)]
    
    def plot(self, proj='xy'):
        import matplotlib.pyplot as plt
        OpticHolder.plot(self, proj)
        x, y = self._projected_footprint(self._corners(self.slit_bb), proj=proj)
        plt.plot(y, x, 'k-')
        x, y = self._projected_footprint(self._corners(-np.array(self.slit_bb)), proj=proj)
        plt.plot(y, x, 'k-')


class KCB1EC(OpticHolder):
    bbox = [-24.1, 24.1, -24.1, 24.1, -24.1, 24.1]
    hole_pos = [24.1 - 17.8, -24.1 + 17.8]
    
    def __init__(self, optic, **kwargs):
        from pyoptic2.elements import rotation_matrix
        up = np.array(kwargs.get('up', [0, 0, 1]))
        orientation = np.array(kwargs.pop('orientation', optic.orientation))
        orientation = np.dot(rotation_matrix(up, pi / 4.), orientation)
        OpticHolder.__init__(self, optic, orientation=orientation, **kwargs)


class DFM1(KCB1EC): #override KCB1EC to get the 45 degrees bit right
    bbox = [-25.4, 25.4, -25.4, 25.4, -30, 30]
    hole_pos = [0, 0]
    hole_radius = 3.5


class CS2100M(OpticHolder):
    #sCMOS camera
    bbox = [-30, 30, -9.5, 53, -25, 25]


class TTL200(OpticHolder):
    bbox = [-18, 18, -18, 18, -14, 14]
    
    def _cutout(self, plane='xy', tol=0.5):
        templ = """
        union(){
            translate([0,0,-5])cylinder(r=39, h=10);
            cylinder(r=36.2, h=22.5);
            translate([0,0,-5])cylinder(r=34, h=60);
        }
        """


class TEXT(OpticHolder):
    #T-Slot extrusion
    def __init__(self, position, length=100, width=40, **kwargs):
        bbox = [-width / 2, width / 2, 0, length, -width / 2, width / 2]
        OpticHolder.__init__(self, optic=position, bbox=bbox, **kwargs)
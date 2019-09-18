import numpy as np

class Placement(object):
    """ Class to describe the location and orientation of a volume """

    def __init__(self, location, orientation) :
        self.location    = np.squeeze(np.array(location, 'f'))
        self.orientation = np.squeeze(np.array(orientation, 'f'))
        self.orientation /= np.linalg.norm(self.orientation)

    def __str__(self) :
        s  = 'Placement\n' 
        s += 'Placement.location       : '+str(self.location) + '\n'
        s += 'Placement.orientation    : '+str(self.orientation) +'\n'
        return s
    
    def offset(self, delta=0, orientation=None):
        return OffsetPlacement(self, delta=delta, orientation=orientation)
    
    
class OffsetPlacement(Placement):
    def __init__(self, old_placement, delta, orientation=None):
        self.old_placement = old_placement
        
        self.delta = np.array(delta, 'f')
        
        if orientation is None:
            self._orientation = None
        else:
            self._orientation = np.squeeze(np.array(orientation, 'f'))
            self._orientation /= np.linalg.norm(self._orientation)
        
    @property
    def location(self):
        if np.ndim(self.delta) > 0:
            #vector delta
            return self.old_placement.location + self.delta
        else:
            #scalar delta - go along orientation
            return self.old_placement.location + self.delta*self.orientation
    
    @property
    def orientation(self):
        if self._orientation is None:
            return self.old_placement.orientation
        else:
            return self._orientation
        
    @orientation.setter
    def orientation(self, value):
        self._orientation = np.array(value, 'f')
    
    
    

P = Placement
OP = OffsetPlacement
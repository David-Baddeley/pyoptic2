import numpy as np

class Placement :
    """ Class to describe the location and orientation of a volume """

    def __init__(self, location, orientation) :
        self.location    = np.array(location, 'f')
        self.orientation = np.array(orientation, 'f')
        self.orientation /= np.linalg.norm(self.orientation)

    def __str__(self) :
        s  = 'Placement\n' 
        s += 'Placement.location       : '+str(self.location) + '\n'
        s += 'Placement.orientation    : '+str(self.orientation) +'\n'
        return s
        

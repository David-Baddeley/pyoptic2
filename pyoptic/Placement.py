import pylab as pl

class Placement :
    """ Class to describe the location and orientation of a volume """

    def __init__(self, location, orientation) :
        self.location    = pl.array(location)
        self.orientation = pl.array(orientation)

    def __str__(self) :
        s  = 'Placement\n' 
        s += 'Placement.location       : '+str(self.location) + '\n'
        s += 'Placement.orientation    : '+str(self.orientation) +'\n'
        return s
        

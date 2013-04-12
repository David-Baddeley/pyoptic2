import pylab as pl 

class Material(object):

    mirror  = 1
    refract = 2

    def __init__(self, type, data = None) :

        self.type = type
        self.n_ = data
    
    def __str__(self) :
        s =  'Material                 : '+str(self.type)+'\n'
        if self.type == self.refract :
            s = 'Material                 : '+str(self.n)
        
        return s
        
    def n(self, wavelength):
        return self.n_
        
Air = Material(Material.refract, 1.0)

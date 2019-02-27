class Material(object):

    MIRROR  = 1
    REFRACT = 2

    def __init__(self, type, data = None) :

        self.type = type
        self.n_ = data
    
    def __str__(self) :
        s =  'Material                 : '+str(self.type)+'\n'
        if self.type == self.REFRACT :
            s = 'Material                 : '+str(self.n)
        
        return s
        
    def n(self, wavelength):
        return self.n_
    
class CauchyGlass(Material):
    """
    Glass model based on the Cauchy model for refractive index
    """
    def __init__(self, type, coeffs=[1.0,]):
        Material.__init__(self, type)
        
        self.coeffs = coeffs
        self._ns = {} #cache refractive index computations
        
    def n(self, wavelength):
        #cache our refractive index to avoid un-necessary calculation
        try:
            n = self._ns[wavelength]
        except KeyError:
            n = 0.0
            for i, c in enumerate(self.coeffs):
                n += c/(wavelength**(2*i))
            
            self._ns[wavelength] = n
    
        return n


class SellmeierGlass(Material):
    """
    Glass model based on the Cauchy model for refractive index
    """
    
    def __init__(self, type, B =[0,0,0], C = [0,0,0]):
        import numpy as np
        Material.__init__(self, type)
        
        self.B = np.array(B)
        self.C = np.array(C)
        self._ns = {} #cache refractive index computations

    def _sellmeier(self, lamb, B, C):
        import numpy as np
        L2 = (lamb / 1e3) ** 2
        return np.sqrt( 1 + ((B*L2)/(L2-C)).sum())
    
    def n(self, wavelength):
        #cache our refractive index to avoid un-necessary calculation
        try:
            n = self._ns[wavelength]
        except KeyError:
            n = self._sellmeier(wavelength, self.B, self.C)
            
            self._ns[wavelength] = n
        
        return n


class AbbeGlass(CauchyGlass):
    """
    Glass model if given only nd, vd.
    
    Fits a Cauchy model.
    """
    _lambda_d = 587.56
    _lambda_F = 486.13
    _lambda_C = 656.27
    def __init__(self, type, nd=1.0, vd=0):
        C = (nd - 1.0)/(vd*(1.0/self._lambda_F**2 - 1.0/self._lambda_C**2))
        B = (nd - C/self._lambda_d**2)
        
        CauchyGlass.__init__(self, type, [B,C])
        
    
    
        
        
Air = Material(Material.REFRACT, 1.0)

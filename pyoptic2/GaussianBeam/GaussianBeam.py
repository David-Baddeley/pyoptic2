# $ Id: $

import numpy as np
import pylab as pl

class GaussianBeam(object):
    """
    Zero-order Gaussian beam
    """    
    def __init__(self, w0z0, k):
        # TODO: handle array input (or not?)-should constructors return arrays?
        """
        GaussianBeam(w0, k)
        GaussianBeam((w0, z0), k)
        Constructor
        """
        self._path = [] # a list of paraxial elements
        #TODO: use a Path class to constrain operations and element types
        self._z0 = 0.0
        if isinstance(w0z0, tuple):
            if len(w0z0) > 0:
                self._w0 = w0z0[0]
            if len(w0z0) > 1:
                self._z0 = w0z0[1]
            if len(w0z0) > 2:
                raise TypeError('Expecting a scalar or 2-element tuple.')
        else:
            self._w0 = w0z0
            
        try:
            self._w0 = float(self._w0)
            self._z0 = float(self._z0)
        except:
            raise TypeError('Expecting a scalar or 2-element tuple.')
            
        self._k = k
        self._zc = confocalDistance(self._w0, k)
        
    def w(self, z, k):
        q = complexBeamParameter(self._zc, z - self._z0)
        invq = 1.0/q;
        w = np.sqrt(-2.0/self._k/invq.imag)
        return w
    
    def q(self, z):
        """
        q(z)
        Complex beam parameter at position Z
        """
        # assume that elements in _path are ordered by Z
        # TODO: make sure that elements are ordered by Z
        if len(self._path) > 0:
            raise Error('Not implemented')
        else:
            q = ComplexBeamParameter(self._zc, z - self._z0) 
        return q
    
    def __rmul__(self, el):
        qin = self.q(el.z)
        qout = el.abcd*qin
        w0, z0 = beamWaist(complex(qout), self._k)
        # TODO: this assumes that element has zero thickness -> add thickness property to ParaxialElement?
        return GaussianBeam((w0, el.z + z0), self._k)
        
    def field(self, r, z):
        q = complexBeamParameter(self._zc, z - self._z0)
        A = fieldAmplitude(q, self._k, r)
        P = fieldPhase(q, self._k, r)
        # TODO: improve efficiency by calculating A and P from private members
        return A*np.exp(P*1j)
        
class ComplexBeamParameter(object):
    """
    Zero-order Gaussian beam complex beam parameter
    """
    def __init__(self, *args):
        # TODO: handle array input (or not?)-should constructors return arrays?
        """
        ComplexBeamParameter(q)
        ComplexBeamParameter(zc, z)
        """
        nargin = len(args)
        if nargin == 1:
            self._q = complex(args[0])
        elif nargin == 2:
            self._q = complexBeamParameter(args[0], args[1])
        else:
            raise Error('Invalid number of input arguments.')
        self._invq = 1.0/self._q
            
    def __rmul__(self, abcd):
        abcd = np.matrix(abcd, dtype=float)
        if (abcd.shape != (2,2)):
            raise TypeError('Invalid ABCD matrix size.')
        q2 = (abcd[0,0]*self._q + abcd[0,1])/(abcd[1,0]*self._q + abcd[1,1])
        return ComplexBeamParameter(q2.imag, q2.real)
        
    def __complex__(self):
        return self._q
        
    def __repr__(self):
        return repr(self._q)
        
    def __eq__(self, other):
        return self._q == complex(other)
        
    def w(self, k):
        w = np.sqrt(2.0/(-k*self._invq.imag))
        return w
        
class ParaxialElement(object):
    """
    An optical element which can be modeled by an ABCD matrix
    """
    def __init__(self, abcd, z):
        """
        ParaxialElement(ABCD, Z) element at position Z with ABCD matrix 
        """
        #TODO: perhaps define this more strictly
        self.abcd = np.matrix(abcd, dtype=float);
        if (self.abcd.shape != (2,2)):
            raise TypeError('Invalid ABCD matrix size.')
        self.z = float(z)
        
class ThinLens(ParaxialElement):
    """
    A thin lens
    """
    def __init__(self, f, z0):
        """
        ThinLens(F, Z0) thin lens at position Z0, focal length F
        """
        self._f = float(f)
        ParaxialElement.__init__(self, [[1.0, 0.0],[-1.0/self._f, 1]], z0)
        
    def transformBeamWaist(self, w0z0, k):
        """
        W1, Z1 = transformBeamWaist(W0, K)
        W1, Z1 = transformBeamWaist((W0, Z0), K)        
        Input to output Gaussian beam waist transformation for input beam with
        beam waist W0 at position Z0 and wavenumber K
        """
        _w0 = w0z0
        _z0 = 0
        if isinstance(w0z0, tuple):
            if len(w0z0) > 0:
                _w0 = w0z0[0]
            if len(w0z0) > 1:
                _z0 = w0z0[1]
            if len(w0z0) > 2:
                raise TypeError('Expecting a scalar or 2-element tuple.')
        zc = np.array(confocalDistance(_w0, k), dtype=float)/self._f
        din = np.array(self.z - _z0, dtype=float)/self._f
        M = 1.0/np.sqrt((din - 1)**2 + zc**2)
        dout = 1 + (din - 1)/((din - 1)**2 + zc**2);
        return (_w0*M, self.z + dout*self._f)
        

def confocalDistance(w0, k):
    """
    Confocal distance ZC [d] from beam waist radius W0 [d] and wave number K [1/d]
    
    See also: complexBeamParameter, beamRadius
    """
    return (0.5*k)*w0**2

def complexBeamParameter(zc, z):
    """
    Complex beam parameter Q [d] from confocal distance ZC [d] and position Z
    [d] along propagation direction.
    
    See also: confocalDistance
    """
    z = np.array(z, dtype=float)
    tmp1 = pl.zeros(z.shape, np.complex)
    tmp1.real = z
    zc = np.array(zc, dtype=float)
    tmp2 = pl.zeros(zc.shape, np.complex)
    tmp2.imag = zc
    return tmp1 + tmp2
    
def radiusOfCurvature(q):
    """
    Radius of curvature R [d] from complex beam parameter Q [d]
    
    See also: complexBeamParameter
    """
    invq = 1.0/q;
    idx = invq.real==0;
    if (idx.ndim == 0):
        if (idx == True):
            R = np.inf
        else:
            R = 1.0/invq.real
    else:
        R = pl.zeros(q.shape, np.double)
        R[idx] = np.inf
        idx = ~idx
        R[idx] = 1.0/invq[idx].real;
    return R

def beamRadius(q, k):
    """
    Beam radius W [d] from complex beam parameter Q [d] and wave number K [1/d]
    
    See also: complexBeamParameter
    """
    lam = 2.0*np.pi/k;
    invq = 1.0/q;
    w = np.sqrt(lam/(-np.pi*invq.imag))
    return w

def beamWaistRadius(q, k):
    """
    Beam waist radius W0 [d] from complex beam parameter Q [d] and wave number K
    [1/d].
    
    See also: complexBeamParameter, beamWaistPosition
    """
    zc = q.imag
    w0 = np.sqrt(2*zc/k)
    return w0

def beamWaistPosition(q):
    """
    Beam waist position Z0 [d] from complex beam parameter Q [d] and wave number K
    [1/d].
    
    See also: complexBeamParameter, beamWaistRadius
    """
    z = q.real
    return -z

def beamWaist(q, k):
    """
    Beam waist radius and position as tuple (W0, Z0) [d] from complex beam
    parameter Q [d] and wave number K [1/d].
    
    See also: complexBeamParameter, beamWaistRadius, beamWaistPosition
    """
    zc = q.imag
    w0 = np.sqrt(2*zc/k)
    z = q.real
    return (w0, -z)

def fieldAmplitude(q, k, r):
    """
    Scalar field amplitude A relative to beam waist center from complex beam
    parameter Q [d], wave number K [1/d], at radial distance R [d] from the
    propagation axis.
    
    See P.F. Goldsmith, "Quasioptical Systems", Section 2.1, p.15
    
    See also: complexBeamParameter, beamRadius, fieldPhase
    """
    w = beamRadius(q, k)
    w0 = beamWaistRadius(q, k)
    A = w0/w*np.exp(-r**2/w**2)
    return A

def fieldPhase(q, k, r):
    """
    Field phase P [rad] referred to beam waist, at radial distance R [d] from
    the propagation axis.

    See P.F. Goldsmith, "Quasioptical Systems", Section 2.1, p.15
    
    See also: complexBeamParameter, radiusOfCurvature, fieldAmplitude
    """
    R = radiusOfCurvature(q)
    P = -0.5*k*r**2/R
    z = q.real
    zc = q.imag
    P0 = -k*z + np.arctan(z/zc)
    # TODO: sign of result assumes jwt time convention
    return np.mod(P + P0, 2*np.pi)
    
def test(dryTest=True):
    dryTest = bool(dryTest)
    pi = np.pi
    print "Running scalar tests ..."
    w0 = 8.0
    lam = 3.0
    zc = confocalDistance(w0, 2*pi/lam)
    assert zc == pi*w0**2/lam
    z = 0.0;
    q = complexBeamParameter(zc, z)
    assert q.real==z
    assert q.imag==zc
    R = radiusOfCurvature(q)
    assert R == np.inf
    w = beamRadius(q, 2*pi/lam)
    assert abs(w-w0) < 1.0e-15
    print "Pass"
    print "Running array tests ..."
    w0 = 8.0
    lam = 3.0
    zc = confocalDistance(w0, 2*pi/lam)
    assert zc == pi*w0**2/lam
    z = np.array([0, zc])
    q = complexBeamParameter(zc, z)
    assert q.real[0]==z[0]
    assert q.real[1]==zc
    assert q.imag[0]==zc
    assert q.imag[1]==zc
    R = radiusOfCurvature(q)
    assert R[0] == np.inf
    assert R[1] == 2*zc
    w = beamRadius(q, 2*pi/lam)
    assert abs(w[0]-w0) < 1.0e-15
    assert abs(w[1]-np.sqrt(2.0)*w0) < 1.0e-12
    assert (beamWaistRadius(q[1], 2*pi/lam) - w0) < 1.0e-15
    if not dryTest: pl.figure(); pl.plot(z, w)
    print "Running field tests ..."
    r = np.array([0, w0, 2*w0])
    A = fieldAmplitude(q[0], 2*pi/lam, r)
    assert abs(A[0]-1) < 1.0e-15
    assert abs(A[1]-np.exp(-1.0)) < 1.0e-15
    assert abs(A[2]-np.exp(-4.0)) < 1.0e-15
    P = fieldPhase(q[0], 2*pi/lam, r)
    assert abs(P[0]) < 1.0e-15
    assert abs(P[1]) < 1.0e-15
    # TODO: add some characteristic points along Z
    print "Testing GaussianBeam class"
    gb = GaussianBeam(beamWaist(q[0], 2*pi/lam), 2*pi/lam)
    (R,Z)= pl.meshgrid(np.arange(-24,24), np.arange(0,100))
    if not dryTest: pl.figure(); pl.imshow(abs(gb.field(R, Z)))
    print "Testing ComplexBeamParameter class"
    d = 10
    q = gb.q(0)
    abcd = np.matrix([[1, d],[0, 1]], dtype=float)
    qo = abcd*q
    assert qo == gb.q(d)
    print "Testing ParaxialElement class"
    el = ParaxialElement(abcd, 0)
    gb2 = el*gb
    print gb2.q(0)
    print gb.q(d)
    assert gb2.q(0)==gb.q(d)
    print "Pass"

# if __name__ == "__main__":

    

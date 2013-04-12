# $ Id: $
"""
The formulas used in this package are based on P. Goldsmith, Quasioptical Systems

Physical units of input and output parameters are based on distance units 'd',
denoted as [d]

 notes:
 -* Express properties in distance or reciprocal distance units
 -* Gaussian beam could keep a 'path' property that describes
     how it propagated (list of ABCD matrices) --> from current position
     always possible to calculate back in any direction
 -* keep basic formulas accessible as normal functions (e.g. zc = confocal distance)

 properties:
  beamWaistRadius (can change by propagation)
  position (changes by propagation)
  path = list of ABCD matrices
  waveLength or waveNumber

 calculated properties:
  confocalDistance (depends on beam waist)
  complexBeamParameter (depends on beam waist & distance from beam waist)

 methods:
  propagate(path)
  tfield(x, y), where x, y are orthogonal to beam axis
  tamp(x, y)
  tphase(x, y)
"""

#__all__ = ['GaussianBeam']

from GaussianBeam import *

# TODO exclude import of locally imported modules like numpy, pylab
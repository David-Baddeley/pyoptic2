#!/usr/bin/env python
import unittest
import numpy as np

from GaussianBeam import *

class GaussianBeamTest(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def testFunScalarArgs(self):
        w0 = 8.0
        lam = 3.0
        k = 2*np.pi/lam
        zc = confocalDistance(w0, k)
        self.assertEqual(zc, np.pi*w0**2/lam)
        z = 0.0;
        q = complexBeamParameter(zc, z)
        self.assertEqual(q.real, z)
        self.assertEqual(q.imag, zc)
        R = radiusOfCurvature(q)
        self.assertEqual(R, np.inf)
        w = beamRadius(q, k)
        self.assertAlmostEqual(w, w0)
    
    def testFunArrayArgs(self):
        w0 = 8.0
        lam = 3.0
        k = 2*np.pi/lam
        zc = confocalDistance(w0, k)
        z = np.array([0, zc])
        q = complexBeamParameter(zc, z)
        self.assertEqual(q.real[0], z[0])
        self.assertEqual(q.real[1], zc)
        self.assertEqual(q.imag[0], zc)
        self.assertEqual(q.imag[1], zc)
        R = radiusOfCurvature(q)
        self.assertEqual(R[0], np.inf)
        self.assertEqual(R[1], 2*zc)
        w = beamRadius(q, k)
        self.assertAlmostEqual(w[0], w0)
        self.assertAlmostEqual(w[1], np.sqrt(2.0)*w0)
        self.assertAlmostEqual(beamWaistRadius(q[1], k), w0)
        #if not dryTest: pl.figure(); pl.plot(z, w)
        #print "Running field tests ..."
        r = np.array([0, w0, 2*w0])
        A = fieldAmplitude(q[0], k, r)
        self.assertAlmostEqual(abs(A[0]), 1)
        self.assertAlmostEqual(abs(A[1]), np.exp(-1.0))
        self.assertAlmostEqual(abs(A[2]), np.exp(-4.0))
        P = fieldPhase(q[0], k, r)
        self.assertAlmostEqual(abs(P[0]), 0)
        self.assertAlmostEqual(abs(P[1]), 0)
        # TODO: add some characteristic points along Z
    
    def testClasses(self):
        w0 = 8.0
        lam = 3.0
        k = 2*np.pi/lam
        zc = confocalDistance(w0, k)
        z = np.array([0, zc])
        q = complexBeamParameter(zc, z)        
        #print "Testing GaussianBeam class"
        gb = GaussianBeam(beamWaist(q[0], k), k)
        (R,Z)= pl.meshgrid(np.arange(-24,24), np.arange(0,100))
        #if not dryTest: pl.figure(); pl.imshow(abs(gb.field(R, Z)))
        d = 10.0
        q = gb.q(0)
        abcd = np.matrix([[1, d],[0, 1]], dtype=float)
        qo = abcd*q
        self.assertEqual(qo, gb.q(d))
        #print "Testing ParaxialElement class"
        el = ParaxialElement(abcd, 0)
        gb2 = el*gb
        #print gb2.q(0)
        #print gb.q(d)
        self.assertEqual(gb2.q(0), gb.q(d))
    
    def testTelescope(self):
        import matplotlib
        matplotlib.use('AGG')
        import matplotlib.mlab as ml
        import pylab as pl
        import time        
        w0 = 8.0
        k = 2*np.pi/3.0
        gb = GaussianBeam(w0, k)
        lens = ThinLens(150, 150)
        gb2 = lens*gb
        self.assertAlmostEqual(gb2._z0, gb._z0 + 2*150.0)
        lens2 = ThinLens(300, 600)
        gb3 = lens2*gb2
        self.assertAlmostEqual(gb3._z0, gb2._z0 + 2*300.0)
        self.assertAlmostEqual(gb._w0, gb3._w0/2.0)
        z = np.arange(0, 150)
        z2 = np.arange(150, 600)
        z3 = np.arange(600, 900)
        pl.plot(z, gb.w(z, k), z2, gb2.w(z2, k), z3, gb3.w(z3, k))
        pl.grid()
        pl.xlabel('z')
        pl.ylabel('w')
        pl.savefig('testTelescope1.png')
        time.sleep(0.1)
        pl.close('all')        
    
class ParaxialElementTest(unittest.TestCase):
    def setUp(self):
        self.gb = GaussianBeam((8.0, -50.0), 2*np.pi/3.0)
        
    def testThinLens(self):
        import matplotlib
        matplotlib.use('AGG')
        import matplotlib.mlab as ml
        import pylab as pl
        import time
        lens = ThinLens(150.0, 100.0)
        gbOut = lens*self.gb
        self.assertAlmostEqual(gbOut._z0, 250.0)
        (w1, z1) = lens.transformBeamWaist((self.gb._w0, self.gb._z0), self.gb._k)
        self.assertAlmostEqual(gbOut._w0, w1)
        
        f = 150.0
        lens = ThinLens(f, 0)
        w0, z0 = ml.meshgrid([4.0, 8.0, 16.0], np.arange(-3*f, f))
        (w1, z1) = lens.transformBeamWaist((w0, z0), self.gb._k)
        #h = pl.figure()
        pl.plot(z0/f, z1/f)
        pl.grid()
        pl.xlabel('d_{in} [f]')
        pl.ylabel('d_{out} [f]')
        pl.savefig('testThinLens1.png')
        time.sleep(0.1)
        pl.close('all')
    
if __name__ == '__main__':
    #unittest.main()
    suiteGB = unittest.TestLoader().loadTestsFromTestCase(GaussianBeamTest)
    suitePE = unittest.TestLoader().loadTestsFromTestCase(ParaxialElementTest)
    suite = unittest.TestSuite([suiteGB, suitePE])
    unittest.TextTestRunner(descriptions=2, verbosity=2).run(suite)

    


import pylab as pl
import Intensity 

def FraunhoferTest(z) :
    s = 0.0005
    i1 = Intensity.Intensity2D(256,-0.005,0.005,256,-0.005,0.005)
    i1.makeRectangularFlatTop(0,0,0.001,0.001)
#    i1.makeEllipticalFlatTop(0,0,0.001,0.001)
#    i1.makeGaussian(0,0,s,s)
    i1.plot(1)
    i1.calculate()
    i2 = Fraunhofer(i1,z)
    i2.plot(2)


#    za = []
#    sa = []

#    for z in pl.arange(0.1,200.2,10) :
#        i2 = Fraunhofer(i1,z)
    #    i2.plot(2)
#        i2.calculate()
#        za.append(z)
#        sa.append(i2.xrms)

#    i2.plot(1)
#    pl.figure(2)
#    pl.plot(za,sa)
#    pl.axhline(s)

def Fraunhofer(i, z) :
    print "Propagation:Fraunhofer"
    ft = pl.fftshift(pl.fftn(pl.fftshift(i.i)))
    dx = i.wl*z/(i.nx*i.dx)
    dy = i.wl*z/(i.ny*i.dy)
    po = pl.exp(1j*2*pl.pi/i.wl*i.dx*i.dx)/(1j*i.wl*z)
    p = pl.arange(0,i.nx)-(i.nx+0.5)/2.0
    q = pl.arange(0,i.ny)-(i.ny+0.5)/2.0
    [pp,qq] = pl.meshgrid(p,q)
    pm = pl.exp(1j*pl.pi/(i.wl*z)*((pp*dx)**2+(qq*dy)**2))
    i2 = Intensity.Intensity2D(i.nx,-i.nx*dx/2,i.nx*dy/2,i.ny,-i.ny*dy/2,i.ny*dy/2)
    i2.i = po*pm*ft
    return i2
    
    print "Propagation:Fraunhofer>",dx,dy,i.nx*dx,i.ny*dy
    
def FraunhoferRectangularAperture(xw,yw,z) :
    print "Propagation:FraunhoferRectangularAperture"

def FraunhoferEllipticalAperture(xw,yw,z) :
    print "Propagation:FraunhoferEllipticalAperture"

def FresnelSingleTest(z) :
    i1 = Intensity.Intensity2D(128,-1,1,128,-1,1)
    
def FresnelSingle(i,z) :
    print "FresnelSingle"
    
def FresnelDouble(i,z) :
    print "FresnelDouble" 



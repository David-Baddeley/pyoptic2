# README #

This is a fork to keep track of my changes and bugfixes to the google code project pyoptic (code.google.com/p/pyoptic) which now seems to be unmaintained. Pyoptic is a python-based raytracing/optical CAD package which works from a programmatic definition of the optical system. It was originally envisaged as a 'toy' tracer for educational purposes, but is arguably capable of being used in some (limited) production uses. I use it for testing the arrangement of optical components in large systems such as microscopes.

### Features ###

* Spherical surfaces
* Mirrors
* 3D Visualisation
* Thin lenses (new in fork)
* Dispersion (new in fork)
* Import of Zemax lens files (e.g. for Thorlabs lenses) - reads geometry and glass information (new in fork)

There are also a number of bug fixes and visualisation tweaks that are new in the fork.

### How do I get set up? ###

Pyoptic requires a number of dependencies:

* Python 2.7
* Numpy, Scipy, Matplotlib, Matplotlib
* Mayavi2

A python distribution such as Enthought Canopy, WinPython, Python(x,y), or Anaconda ought to provide all of these.

After downloading and installing python & dependencies, pyoptic can be installed by executing:
python setup.py install in the pyoptic base directory.
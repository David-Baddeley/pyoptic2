# README #

This is a fork to keep track of my changes and bugfixes to the google code project pyoptic (code.google.com/p/pyoptic) which now seems to be unmaintained. Pyoptic is a python-based raytracing/optical CAD package which works from a programmatic definition of the optical system. It was originally envisaged as a 'toy' tracer for educational purposes, but is arguably capable of being used in some (limited) production uses. I use it for testing the arrangement of optical components in large systems such as microscopes.

By now, things have divereged a lot from the original pyoptic, motivating the change to pyoptic2.

pyoptic2 is released under the [GNU GPL v3 license](https://choosealicense.com/licenses/gpl-3.0/)

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

* Python >= 3.6
* Numpy, Scipy, Matplotlib
* Mayavi2 (optional, for 3D visualisation)

A python distribution such as Enthought Canopy, WinPython, Python(x,y), or Anaconda ought to provide all of these.

After downloading and installing python & dependencies, pyoptic can be installed by executing:
`python setup.py install` in the pyoptic base directory.

### Getting started ###

The examples directory has a couple of ipython notebooks with comments:

*pyoptic_testsv2.ipynb* : Sets up a simple, folded 4-f system using catalog lenses. This is a good place to start.
*olympus_25x_objective.ipynb*: Simulates a microscope objective (based on patent literature). This is a bit more involved, 
and less well documented, but is an example of how the package can be used for something practical. It shows how to construct
systems from individual surfaces.

### Optional dependencies ###

These packages are needed to automatically grab Zemax files from Thorlabs:

* Requests, Beautifulsoup4, Lxml

### Comparison to other optical CAD packages/ What is missing ###

There are a few things that you might find in professional optical CAD that pyoptic does not yet do

* Aspherical surfaces (this would be resonably simple to add)
* Polarization
* Fresnel coefficients (or any other form of intensity tracking)
* Non-sequential ray tracing
* Both reflection and transmission at the same surface (this can be spoofed in order to model, e.g. dichroic mirrors 
  separating excitation and detection paths in a microscope by creating a copy of the system up until the split) 
  
On the other side of the coin, there are a few things which we are particularly well equipped to do. Because 
system definition is programatic, it is very easy to generate parametric systems, and, by leveraging the 
numpy/scipy/matplotlib stack, to generate plots of characteristics which might be specific to your system.

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 11:57:53 2013

@author: David Baddeley
"""
import lzw
import urllib
import os
from bs4 import BeautifulSoup
import numpy as np

import read_agf

def decompZar(filename):
    """zemax .zar files are a series of concatenated, LZW compressed text-format
    files.
    
    This function attempts to locate and decompress each of the sections, returning
    a dictionary of strings where the keys are the origninal(uncompressed) file names.    
    """
    f = open(filename, 'rb')
    raw = f.read()
    f.close()
    
    #delimiter seems to be first 4 bytes
    delim = raw[:4]
    
    #sections seem to be separated by at least 4 nulls, often more
    #the length comparison deals with multiple consecutive nulls
    sections = [s for s in raw.split('\0\0\0\0') if len(s) > 0]
    
    componentFiles = {}
    
    for i, s in enumerate(sections):
        #compressed data is preceeded by a plain-text file name
        if s.endswith('.LZW'):
            #this is a filename
            compFName = s[:-4]
            
            compData = ''.join([c for c in lzw.decompress(sections[i+1].lstrip('\0'))])
            
            componentFiles[compFName] = compData
        elif len(s) < 256:
            try: 
                #print s
                s = (s.lstrip('\0') + '\0').decode('utf16')
            except:
                s = ''
                
            #print s.encode('ascii', 'replace')
            if s.endswith('.LZW'):
                compFName = s[:-4]
        
                compData = ''.join([c for c in lzw.decompress(sections[i+1].lstrip('\0'))])
                
                if (len(compData) % 2):
                    #make sure we can decode as utf16
                    compData+='\0'
                compData = compData.decode('utf16').encode('ascii', 'replace')
                
                componentFiles[compFName] = compData
            
    return componentFiles
 
   
class DefaultGlass(object):
    #for compatibility with pyoptic
    REFLECT = 1
    REFRACT = 2
    
    def __init__(self, n=1.0, typ=2) :
        self.type = typ
        self.data = {'name':'air'}
        self._n = n
    
    def __str__(self) :
        return   'Material                 : ' + self.data['name']+'\n'

        
    def n(self, wavelength):
        return self._n
        
AIR = DefaultGlass()
    
class Glass(DefaultGlass):
    def __init__(self, data, typ=2) :

        self.type = typ
        self.data = data
        self.ns = {}
    
    def n(self, wavelength):
        #cache our refractive index to avoid un-necessary calculation
        try:
            n = self.ns[wavelength]
        except KeyError:
            n = read_agf.get_index(self.data, wavelength)
            self.ns[wavelength] = n
        
        return n
    
class Surface(object):
    def __init__(self, lines, glasses):
        self.glassName=None
        self.glass = AIR
        self.parseLines(lines, glasses)
        
        
    
    def parseLines(self, lines, glasses):
        self.lines = []
        while (len(lines) > 0) and (lines[0].startswith('  ') or lines[0] == ''):
            l = lines.pop(0)
            self.lines.append(l)
            l = l.lstrip()
            code = l[:4]
            data = l[5:]
            if code == 'TYPE':
                self.type = data
            elif code == 'CURV':
                self.curv = data.split()
                c = float(self.curv[0])
                if not c == 0:
                    self.radius = 1./float(self.curv[0])
                else:
                    self.radius = None
            elif code == 'GLAS':
                self.glas = data.split()
                self.glassName = self.glas[0]
                self.glass = Glass(glasses[self.glassName])
            elif code == 'DISZ':
                self.disz = float(data)
            elif code == 'DIAM':
                self.diam = data.split()
                self.diameter = float(self.diam[0])
                

    
class ZMX(object):
    def __init__(self, data, glasses):
        self.surfaces = []
        self._parseZMX(data, glasses)
        
        
    def _parseZMX(self,data, glasses):
        lines = data.splitlines()
        self.lines = list(lines)
        while len(lines) > 0:
            l = lines.pop(0)
            if l.startswith('SURF'):
                self.surfaces.append(Surface(lines, glasses))  
                
    
    def toPyOptic(self, position, direction=[0,0,1], fb=None, flip = False, f = None, orientation=[0,1,0]):
        import pyoptic2 as pyo

        d = np.array(direction)
        pos = np.array(position)

        plc = pyo.Placement(pos, d)
        
        return self.to_pyoptic(plc, fb, flip, f, orientation)
    
    def to_pyoptic(self, placement, fb=None, flip = False, f = None, orientation=[0,1,0]):
        import pyoptic2 as pyo
    
        #we're only interested in real surfaces
        surfs = self.surfaces[1:-1]
        
        l = sum([s.disz for s in surfs[:-1]])
        #print(l)
        
        if fb is None and f is None:
            #calculate length of lens
            
            #midpoint at l/2
            z0 = -l/2
        elif fb is None:
            z0 = float(f) - sum([s.disz for s in surfs])
        else:
            z0 = -float(fb) - l + float(f)
        
        outSurfs = []
        
        
        
        if flip:
            z0 = -z0 - l
            i_s = range(len(surfs))
            for i in i_s[::-1]:
                #print(i)
                s = surfs[i]
                glass = self.surfaces[i].glass
                if s.radius is None:
                    #planar
                    outSurfs.append(pyo.PlaneSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement= pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=glass))
                elif s.type == 'TOROIDAL':
                    #cylindrical lens
                    outSurfs.append(pyo.CylindricalSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement = pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=glass,
                                                     curvature_radius=-s.radius,
                                                     curvature_axis=orientation))
                else:
                    outSurfs.append(pyo.SphericalSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement = pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=glass,
                                                     curvature_radius=-s.radius))
                #print((z0, s.disz))
                z0 += self.surfaces[i].disz         
        else:
            for i, s in enumerate(surfs):
                if s.radius is None:
                    #planar
                    outSurfs.append(pyo.PlaneSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement= pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=s.glass))
                elif s.type == 'TOROIDAL':
                    outSurfs.append(pyo.CylindricalSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement= pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=s.glass,
                                                     curvature_radius=s.radius,
                                                     curvature_axis=orientation))
                else:
                    outSurfs.append(pyo.SphericalSurface(name='ZMX_%d' % i,
                                                     shape=pyo.SHAPE_CIRC,
                                                     dimension=np.ones(3) * float(s.diam[0]),
                                                     placement= pyo.OffsetPlacement(placement, (z0 - 0)),
                                                     material=s.glass,
                                                     curvature_radius=s.radius))
                #print((z0, s.disz))
                z0 += s.disz
                
            
        return outSurfs
    
    def __call__(self, *args, **kwargs):
        return self.to_pyoptic(*args, **kwargs)
        
#keep glass database global (lets us read files with broken glass info if we have read good files first)
glasses = {}
                
def readZar(filename):
    components = decompZar(filename)
    
    #glasses = {}
    
    #find all the glass catalogs    
    glass_cats = [n for n in components.keys() if n.endswith('.AGF')]
    for gc in glass_cats:
        glasses.update(read_agf.read_cat_from_string(components[gc]))
        
    #print((glasses.keys()))
        
    #find all components (probably only one)
    zmxs = [ZMX(components[n], glasses) for n in components.keys() if (n.endswith('.ZMX') or n.endswith('.zmx'))]
    
    return zmxs, glasses


def read_cached_zar(filename):
    import shelve
    #rom PYME.misc import zemax
    
    shv = shelve.open('lensdb2.shv')
    try:
        lens = shv[filename]
    except KeyError:
        lens = readZar(filename)[0][0]
        shv[filename] = lens
    
    shv.close()
    
    return lens
    
def load_thorlabs_zar(partnumber, cached=True):
    # Where on ThorLab's website do we look for parts?
    root = 'https://www.thorlabs.com'
    extension = '/thorproduct.cfm?partnumber='

    # Where do we save thorlabs cached files?
    home = os.path.expanduser("~")
    if not os.path.exists(home + '/.pyoptic'):
        os.makedirs(home + '/.pyoptic')
    if not os.path.exists(home + '/.pyoptic/thorlabs'):
        os.makedirs(home + '/.pyoptic/thorlabs')
    save_dir = home + '/.pyoptic/thorlabs/'

    # Load up the part
    address = root + extension + partnumber
    response = urllib.urlopen(address)
    data = response.read()
    soup = BeautifulSoup(data)

    # Grab the Zemax file
    found = False
    for link in soup('a'):
        if link.get('alt', '') == 'Zemax':
            download_link = root + link['href']
            saved_file_name = save_dir + partnumber + '-Zemax.zar'
            urllib.urlretrieve(download_link, saved_file_name)
            found = True
    
    if not found:
        raise ValueError('Part not found.')

    if cached:
        lens = read_cached_zar(saved_file_name)
    else:
        lens = readZar(saved_file_name)[0][0]

    return lens
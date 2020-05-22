# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 11:57:53 2013

@author: David Baddeley
"""
from . import lzw
import numpy as np

from . import read_agf

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

        zvs = -np.cumsum([s.disz for s in surfs[::-1]])[::-1]
        
        if fb is None and f is None:
            #calculate length of lens
            
            #midpoint at l/2
            z0 = -l/2
            zvs = zvs - zvs.mean()
        elif fb is None:
            z0 = float(f) - sum([s.disz for s in surfs])
            zvs = float(f) + zvs
        else:
            z0 = -float(fb) - l + float(f)
            zvs = float(f) + zvs - zvs[-1] - float(fb)
        
        outSurfs = []
        
        #print('z0: %f, l: %f' % (z0, l))
        
        #print(zvs)
        
        if flip:
            z0 = -z0 - l
            #print('z0: %f' % z0)
            i_s = range(len(surfs))
            for i in i_s[::-1]:
                #print(i)
                s = surfs[i]
                
                kwargs = dict(name=('ZMX_%d' % i),
                              shape=pyo.SHAPE_CIRC,
                              dimension=np.ones(3) * float(s.diam[0]),
                              placement= pyo.OffsetPlacement(placement, -zvs[i]),
                              material=self.surfaces[i].glass, material2=s.glass)
                
                if s.radius is None:
                    #planar
                    outSurfs.append(pyo.PlaneSurface(**kwargs))
                elif s.type == 'TOROIDAL':
                    #cylindrical lens
                    outSurfs.append(pyo.CylindricalSurface(curvature_radius=-s.radius,
                                                     curvature_axis=orientation, **kwargs))
                else:
                    outSurfs.append(pyo.SphericalSurface(curvature_radius=-s.radius, **kwargs))
                #print((z0, s.disz))
                z0 += self.surfaces[i].disz
            #print('z0: %f' % z0)
        else:
            z0 = -z0 -l
            #print('z0: %f' % z0)
            for i, s in enumerate(surfs):
                kwargs = dict(name='ZMX_%d' % i,
                             shape=pyo.SHAPE_CIRC,
                             dimension=np.ones(3) * float(s.diam[0]),
                             placement= pyo.OffsetPlacement(placement, zvs[i]),
                             material=s.glass, material2=self.surfaces[i+1].glass)
                
                if s.radius is None:
                    #planar
                    outSurfs.append(pyo.PlaneSurface(**kwargs))
                elif s.type == 'TOROIDAL':
                    outSurfs.append(pyo.CylindricalSurface(curvature_radius=s.radius,
                                                     curvature_axis=orientation, **kwargs))
                else:
                    outSurfs.append(pyo.SphericalSurface(curvature_radius=s.radius, **kwargs))
                    
                #print((z0, s.disz))
                z0 += s.disz

            #print('z0: %f' % z0)
                
                
            
        return pyo.ElementGroup(outSurfs, placement)
    
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
        glasses.update(shv.get('glasses', {}))
        lens = readZar(filename)[0][0]
        shv[filename] = lens
        
        #cache glass info
        shv['glasses'] = glasses
    
    shv.close()
    
    return lens
    
def load_thorlabs_zar(partnumber, cached=True):
    """
    Goes to the ThorLabs website and grabs the .zar file for the specified part number.

    Parameters
    ---------- 
        partnumber : str
            ThorLabs part number.
        cached : bool
            Use cached .zar files.

    Returns
    -------
        lens : pyoptic2.util.zemax.ZMX
            Zemax lens.
    """

    import os
    import requests
    from bs4 import BeautifulSoup

    # Where on ThorLab's website do we look for parts?
    root = 'http://www.thorlabs.com'
    extension = '/thorproduct.cfm?partnumber='

    # Where do we save thorlabs cached files?
    home = os.path.expanduser("~")
    if not os.path.exists(home + '/.pyoptic'):
        os.makedirs(home + '/.pyoptic')
    if not os.path.exists(home + '/.pyoptic/thorlabs'):
        os.makedirs(home + '/.pyoptic/thorlabs')
    save_dir = home + '/.pyoptic/thorlabs/'

    save_file_name = save_dir + partnumber + '-Zemax.zar'

    found = True

    # Check if we've already downloaded the file or if we've disabled the cache
    if (not os.path.isfile(save_file_name)) or (not cached):

        found = False

        # Load up the part
        address = root + extension + partnumber
        print(address)
        
        headers = requests.utils.default_headers()
        #print(headers)
        headers['Connection'] = 'close'
        headers['User-Agent'] =  "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/77.0.3865.120 Safari/537.36"
        response = requests.get(address, headers=headers)
        data = response.content
        soup = BeautifulSoup(data, 'lxml')

        #if response.status_code != 200:
        #    #print the error message
        #    print(data)

        # Grab the Zemax file
        for link in soup('a'):
            if link.get('alt', '') in ['Zemax', 'Zemax (ZAR)']:
                download_link = root + link['href']
                # Need to add a header to get requests to work
                #headers = {'Connection': 'close'}
                response = requests.get(download_link, allow_redirects=True, stream=True, headers=headers)
                if response.status_code != 200:
                    print(response.content)
                    
                response.raise_for_status()
                with open(save_file_name, 'wb') as fp:
                    for block in response.iter_content(1024):
                        fp.write(block)
                found = True
                
        if not found:
            print(data)
        
    if not found:
        raise ValueError('Part %s not found.' % partnumber)

    # Use the cache
    if cached:
        lens = read_cached_zar(save_file_name)
    else:
        lens = readZar(save_file_name)[0][0]

    return lens
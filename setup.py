#! /usr/bin/env python

__version__ = '0.9.0'

from setuptools import setup, find_packages
from os import path

if __name__ == '__main__':

    if sys.version_info[0] == 2:
        sys.exit("Sorry, Python 2 is not supported")

    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as file_id:
        long_description = file_id.read()
    short_description = ('An experimental python-based raytracing/optical CAD package which '
                         'works from a programmatic definition of the optical system.')

    setup(name='pyoptic2',
          version=__version__,
          license="GPLv3",
          description=short_description,
          long_description=long_description,
          long_description_content_type='text/markdown',
          author="David Baddeley",
          author_email='david.baddeley.nz@gmail.com',
          packages=find_packages(),
          keywords=['Optics', 'Ray Tracing'],
          url='https://github.com/David-Baddeley/pyoptic2',
          project_urls={"Bug Tracker": "https://github.com/David-Baddeley/pyoptic2/issues",
                        "Source Code": "https://github.com/David-Baddeley/pyoptic2g"},
          classifiers=['Development Status :: 5 - Production/Stable',
                       'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                       'Programming Language :: Python :: 3.6'],
          install_requires=['numpy',
                            'scipy',
                            'matplotlib',
                           ],
          extras_require={'dev': ['twine'],
                          '3d': ['mayavi'], # for 3D visualisation of optical layouts. Will likely also require either wx or QT (but we leave that for mayavi itself to specify)
                          'fetch': ['lxml', 'requests', 'beautifulsoup4'], # to allow automatic download of zars for thorlabs lenses
                          'all': ['mayavi',
                                   'lxml',
                                   'requests',
                                   'beautifulsoup4']},},
          python_requires=">=3.6"
    )

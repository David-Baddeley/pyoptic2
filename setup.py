#!/usr/bin/python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyoptic2',parent_package,top_path)
    config.add_subpackage('pyoptic2')
    
    #config.make_svn_version_py()  # installs __svn_version__.py
    #config.make_config_py()
    #config.get_version('PYME/version.py')
    return config

if __name__ == '__main__':
    import setuptools
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

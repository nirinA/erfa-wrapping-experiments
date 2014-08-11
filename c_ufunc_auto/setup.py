import numpy
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
import erfa_wrapp

def configuration():
    config = Configuration('erfa',
                           parent_package='',
                           top_path=None)
    config.add_extension('_erfa',
                         sources = ['erfa.c',
                                    '_erfamodule.c'],
                         include_dirs = ['.'],
                         library_dirs = [],
                         libraries = ['m'])
    return config

setup(configuration=configuration)

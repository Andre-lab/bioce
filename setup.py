"""
Setup script for distutils currently working with UNIX based systemes

"""
# setup.py
import distutils
import sys
import os
from distutils.core import setup, Extension

os.environ["CC"] = "clang"
os.environ["CXX"] = "clang"

is_64bits = sys.maxsize > 2**32
enable_openmp = False
if sys.platform == 'darwin':
    if not is_64bits:
        # Disable OpenMP
        enable_openmp = False
    else:
        # Newer versions of Darwin don't support openmp
        try:
            darwin_ver = int(os.uname()[2].split('.')[0])
            if darwin_ver >= 12:
                enable_openmp = False
        except:
            print("PROBLEM determining Darwin version")

extensions = [Extension('_vbwSC', ["vbw_sc.i", "VBW_sc.cpp"],
                             swig_opts=["-c++","-I/usr/include/python2.7"],
                             extra_compile_args= ["-fpic", "-fopenmp",
                                              "-O3", "-std=c++11" ],
                             libraries =[ "gsl", "gslcblas", "m"]
                             )]

setup(name = "bioce",
      version = "1.0",
      description='Bayesian inference of conformational ensembles',
      author='Wojtek, Potrzebowski, Jill Trewhella and Ingemar Andre',
      author_email='wojciech.potrzebowski@esss.se',
      url='https://www.andrelab.org',
      ext_package='vbwSC',
      ext_modules=extensions,
      py_modules=['psis','psisloo','stan_models','stan_utility',
                  'statistics','variationalBayesian','fullBayesian','prepareBayesian']
      )

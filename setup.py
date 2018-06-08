# setup.py
import distutils
from distutils.core import setup, Extension

setup(name = "bioce",
      version = "1.0",
      description='Bayesian inference of conformational ensembles',
      author='Wojtek, Potrzebowski, Jill Trewhella and Ingemar Andre',
      author_email='wojciech.potrzebowski@esss.se',
      url='https://www.andrelab.org',
      ext_package='vbwSC',
      ext_modules=[Extension('_vbwSC', ["vbw_sc.i", "VBW_sc.cpp"],
                             swig_opts=["-I/usr/include/python2.7"],
                             extra_postargs= ["-shared", "-fpic", "-fopenmp",
                                              "-O3", "-std=c+11" ],
                             libraries =[ "lgsl", "lgslcblas", "lm"]
                             )],
      py_modules=['psis','psisloo','stan_models','stan_utility',
                  'statistics','variationalBayesian','fullBayesian']
      )

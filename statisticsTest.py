"""
Tests for statistics modules
1) chi2
2) waic
3) psis loo?
"""
from __future__ import print_function

__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"


from __future__ import division, print_function


import unittest
import time

import numpy as np


def load_data(filename="98929.txt"):
    data = np.loadtxt(filename, dtype=np.float64)

if __name__ == '__main__':
    unittest.main()
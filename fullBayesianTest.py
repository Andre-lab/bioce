"""
Unit tests for full Bayesian analysis including:
1) file reading
2) and combining output
"""
from __future__ import division, print_function

import unittest
import time

import numpy as np


def load_data(filename="98929.txt"):
    data = np.loadtxt(filename, dtype=np.float64)

def test_read_file_safe(filename, dtype="float64"):

def test_combine_curve(simulated, weights, scale):

if __name__ == '__main__':
    unittest.main()
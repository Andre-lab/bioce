#!/bin/bash
CONDA_DIR=$HOME/anaconda3/envs/bioce/

swig -python -c++ -o vbw_sc_wrap.cpp vbw_sc.i
clang -O3 -fPIC -c VBW_sc.cpp -fpic -fopenmp -std=c++11 -I$CONDA_DIR/include
clang -O3 -fPIC -c vbw_sc_wrap.cpp -fpic -fopenmp -std=c++11 -I$CONDA_DIR/include/python3.7m -I$CONDA_DIR/include
clang -bundle -undefined dynamic_lookup VBW_sc.o vbw_sc_wrap.o -L$CONDA_DIR/lib -fopenmp -lgsl -lgslcblas -lm -o _vbwSC.so

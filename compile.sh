swig -python -c++ -o vbw_sc_wrap.cpp vbw_sc.i
clang -O3 -fPIC -c VBW_sc.cpp -fpic -fopenmp -std=c++11 -I/Users/wojtek/anaconda3/envs/bioce/include
clang -O3 -fPIC -c vbw_sc_wrap.cpp -fpic -fopenmp -std=c++11 -I/Users/wojtek/anaconda3/envs/bioce/include/python3.7m -I/Users/wojtek/anaconda3/envs/bioce/include
clang -bundle -undefined dynamic_lookup VBW_sc.o vbw_sc_wrap.o -L/Users/wojtek/anaconda3/envs/bioce/lib -fopenmp -lgsl -lgslcblas -lm -o _vbwSC.so

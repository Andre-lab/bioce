/* vbw.i */
%module vbwSC
%include "std_string.i"
%{
  #define SWIG_FILE_WITH_INIT
  #include "VBW_sc.hh"
%}
 %include "VBW_sc.hh" 

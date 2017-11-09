/* vbw.i */
%module vbwSC
%include "std_string.i"
%{
  #define SWIG_FILE_WITH_INIT
  #include "VBW_sc.hh"
 /* Put header files here or function declarations like below */
 /*extern void run_vbw(const int &again, const int &k, const std::string &mdfile,
        const int &N, const std::string &presaxsfile, const std::string &saxsfile, const std::string &saxserrfile,
        const std::string &outfile, const int &nprocs, const double &w_cut);*/
%}

 %include "VBW_sc.hh" 
 /*extern void run_vbw(const int &again, const int &k, const std::string &mdfile,
        const int &N, const std::string &presaxsfile, const std::string &saxsfile, const std::string &saxserrfile,
        const std::string &outfile, const int &nprocs, const double &w_cut);*/
 

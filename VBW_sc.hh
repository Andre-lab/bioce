/*******************************************************************************
Bayesian weight inference for single scattering curve
Author: Wojciech Potrzebowski
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include "gsl_siman.h"

typedef struct {
        size_t size;
        double *alphas;
        double shiftEnergy;
        void *saxsExpPtr;
        void *saxsErrPtr;
        void *saxsEnsPtr;
        void *saxsPrePtr;
        void *saxsMixPtr;
        double saxsScale;
        void *csExpPtr;
        void *csErrPtr;
	    void *csRmsPtr;
        void *csPrePtr;
        void *csMixPtr;
        int numberProcs;
        void *alphaPre;
        bool rosettaPrior;
        bool saxsOn;
        bool chemicalShiftsOn;
        } block;

block * block_alloc(size_t n);

void block_free(block * t);

void block_copy(void *inp, void *outp);

void * block_copy_construct(void *xp);

void block_destroy(void *xp);

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k);
double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N);
double L_function(void *xp);
double L_distance(void *xp, void *yp);
void L_print (void *xp);
void L_take_step(const gsl_rng * r, void *xp, double step_size);
void run_vbw(const int &again, const int &k, const std::string &pre_weight_file, const std::string &structure_energies_file,
        const int &N, const int &n, const std::string &presaxsfile, const int &Ncurves, const std::string &curvesfile,
        const std::string &outfile, const int &nprocs, const double &w_cut, const int &skip_vbw,
        const std::string &precsfile, const std::string &rmscsfile, const std::string &chemical_shifts_fite);

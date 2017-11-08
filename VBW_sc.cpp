/*******************************************************************************
Bayesian weight inference for single scattering curve
Author: Wojciech Potrzebowski
*******************************************************************************/

#include "VBW_sc.hh"

using namespace std;
const double pi = M_PI;

block * block_alloc(size_t n) {
        block * t = (block *) malloc(sizeof(block));
        t->alphas = (double *) malloc ((n+1) * sizeof(double));
        t->size = n;
        return t;
     	}

void block_free(block * t) {
        free(t->alphas);
        free(t);
     	}

void block_copy(void *inp, void *outp) {
       	int i;
       	block * in = (block *) inp;
       	block * out = (block *) outp;

       	for(i=0; i< in->size; i++){
               out->alphas[i] = in->alphas[i];
       	}
       	out->shiftEnergy = in->shiftEnergy;
       	out->size = in->size;
	    out->saxsExpPtr = in->saxsExpPtr;
	    out->saxsErrPtr = in->saxsErrPtr;
	    out->saxsEnsPtr = in->saxsEnsPtr;
	    out->saxsPrePtr = in->saxsPrePtr;
	    out->saxsMixPtr = in->saxsMixPtr;
        out->saxsScale = in->saxsScale;

        out->csExpPtr = in->csExpPtr;
        out->csErrPtr = in->csErrPtr;
	    //TODO: Seems not to be used
        out->csRmsPtr = in->csRmsPtr;
        out->csPrePtr = in->csPrePtr;
        out->csMixPtr = in->csMixPtr;

	    out->numberProcs = in->numberProcs;
	    out->rosettaPrior = in->rosettaPrior;
	    out->alphaPre = in->alphaPre;
	    out->saxsOn = in->saxsOn;
	    out->chemicalShiftsOn = in->chemicalShiftsOn;
	}

void * block_copy_construct(void *xp) {
	block * x = (block *) xp;
	block * y = block_alloc(x->size);
	block_copy(x, y);
	return y;
}

void block_destroy(void *xp){
	block_free( (block *) xp);
}

///////////////////////////////Simulated annealing handling finished////////////

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k) {

	double jsd=0.0, s1=0.0, s2=0.0;
	for (int i=0; i<k; i++) {
		if ( gsl_vector_get(w_a,i) == 0.0 || gsl_vector_get(w_b,i) == 0.0) continue;
		s1 +=  gsl_vector_get(w_a,i)*log2(2*gsl_vector_get(w_a,i)/(gsl_vector_get(w_a,i)+gsl_vector_get(w_b,i)));
		s2 +=  gsl_vector_get(w_b,i)*log2(2*gsl_vector_get(w_b,i)/(gsl_vector_get(w_a,i)+gsl_vector_get(w_b,i)));
	}
	jsd =  0.5*(s1+s2);
	return jsd;
}


double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N)
{
	double tempa = 0.0, tempb = 0.0;
	for( int i = 0; i< N; i++) {
		tempa += gsl_vector_get(saxs_ens,i)*gsl_vector_get(saxs_exp,i)/gsl_vector_get(err_saxs,i);
		tempb += pow(gsl_vector_get(saxs_ens,i),2.0)/gsl_vector_get(err_saxs,i);
	}
	return tempa/tempb;
}

///////////////////Simulated annealing functions////////////////////////////////
double L_function(void *xp)  
{
  //timeval t1, t2;
  //double elapsedTime;
  //gettimeofday(&t1, NULL);

  block *x = (block *) xp;

  bool saxsOn = x->saxsOn;
  bool chemicalShiftsOn = x->chemicalShiftsOn;
  bool rosettaPrior = x->rosettaPrior;
  int nprocs = x->numberProcs;
  gsl_vector *alpha_pre = (gsl_vector *) (x->alphaPre);

  //TODO: Add these variables to structure block

  size_t L = x->size;

  int rep = 0;
  double alpha_zero = 0.0;
  double energy_zero = 0.0;

  double log_gamma_2 = gsl_sf_lngamma(0.5);
  double Lfunc=0.0;
  double fit_saxs=0.0, fit_saxs_mix = 0.0;

  for (int i = 0; i < L; i++) {
	  alpha_zero+=x->alphas[i];
  }

  if (rosettaPrior) {
        //beta = 1.0/kBT in kcal
        double beta = 1.717472947;
        double shiftEnergyExp = exp(-x->shiftEnergy*beta);

        for (int i = 0; i < L; i++) {
	        Lfunc+=gsl_sf_lngamma(shiftEnergyExp*gsl_vector_get(alpha_pre,i))
	        - gsl_sf_lngamma( x->alphas[i] );

	        Lfunc+=(x->alphas[i]-shiftEnergyExp*gsl_vector_get(alpha_pre,i))
            *(gsl_sf_psi(x->alphas[i])-gsl_sf_psi(alpha_zero));

	        energy_zero+=shiftEnergyExp*gsl_vector_get(alpha_pre,i);
	    }
	    Lfunc+= gsl_sf_lngamma(alpha_zero)-gsl_sf_lngamma(energy_zero);
  }
  else {
        Lfunc+= gsl_sf_lngamma(alpha_zero)-gsl_sf_lngamma(L/2);
        for (int i = 0; i < L; i++) {
            Lfunc+=log_gamma_2 - gsl_sf_lngamma( x->alphas[i] );
            Lfunc+=((x->alphas[i]-0.5)*(gsl_sf_psi(x->alphas[i])-gsl_sf_psi(alpha_zero)));
        }
  }

  //Likelihood functions
  if (saxsOn) {

    gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
    gsl_vector *saxs_exp = (gsl_vector *) (x->saxsExpPtr);
    gsl_vector *err_saxs = (gsl_vector *) (x->saxsErrPtr);
    gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);

    double *mix_saxs = (double *) (x->saxsMixPtr);
    double saxs_scale = x->saxsScale;

    size_t N = saxs_exp->size;

    gsl_vector *weightsL = gsl_vector_alloc(N);
    double alpha_ens[N];

    for( int i = 0; i< N; i++) {
	    alpha_ens[i] = 0.0;
	    for (int k = 0; k < L; k++) {
		    alpha_ens[i]+=gsl_matrix_get(saxs_pre,i,k)*x->alphas[k];
	    }
	    gsl_vector_set(weightsL, i, alpha_ens[i]/alpha_zero);
    }
    double saxs_scale_current = SaxsScaleMean(weightsL,saxs_exp,err_saxs,N);

    for( int i = 0; i< N; i++) {
	    fit_saxs += ( pow(saxs_scale_current*alpha_ens[i]/alpha_zero -
	    gsl_vector_get(saxs_exp,i), 2) / pow(gsl_vector_get(err_saxs,i),2) );
    }

    double smix, deltamix;
    int i_ind,j_ind;
   
    //gettimeofday(&t1, NULL);
    #pragma omp parallel for \
    default(none) shared(L,x,mix_saxs, alpha_zero)\
    private (i_ind, j_ind, smix, deltamix) \
    num_threads(nprocs) \
    schedule(dynamic,16) \
    reduction(+:fit_saxs_mix)

    for(i_ind = 0; i_ind < L; i_ind++) {
  	    for ( j_ind = i_ind; j_ind < L; j_ind++) {
		    smix = mix_saxs[L*i_ind+j_ind];
		    deltamix = (i_ind!=j_ind) ? -2*x->alphas[i_ind]*x->alphas[j_ind] :
            x->alphas[i_ind]*(alpha_zero - x->alphas[i_ind]);
            fit_saxs_mix += deltamix * smix;
		}
    }
    fit_saxs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
    Lfunc+=0.5*(fit_saxs+fit_saxs_mix);
   
    gsl_vector_free(weightsL);
    //gsl_vector_free(saxs_ens);
    //gsl_vector_free(saxs_exp);
    //gsl_vector_free(err_saxs);
    //gsl_matrix_free(saxs_pre);
    //free(mix_saxs);

  }
  //This may be semi-optimal but is good for modularization
  if (chemicalShiftsOn) {
    //cout<<"Chemical shifts on"<<std::endl;
    gsl_vector *cs_exp = (gsl_vector *) (x->csExpPtr);
    gsl_vector *cs_err = (gsl_vector *) (x->csErrPtr);
    gsl_vector *cs_rms = (gsl_vector *) (x->csRmsPtr);
    gsl_matrix *cs_pre = (gsl_matrix *) (x->csPrePtr);
    double *mix_cs = (double *) (x->csMixPtr);
    size_t n = cs_exp->size;
    double cs_alpha_ens[n];
    double fit_cs=0.0, fit_cs_mix = 0.0;
    double cs_mix, deltamix;
    int i_ind,j_ind;

    for( int i = 0; i< n; i++) {
        cs_alpha_ens[i] = 0.0;
        for (int k = 0; k < L; k++) {
            cs_alpha_ens[i]+=gsl_matrix_get(cs_pre,i,k)*x->alphas[k];
        }
        fit_cs += ( pow(cs_alpha_ens[i]/alpha_zero - gsl_vector_get(cs_exp,i), 2)
        / ( pow(gsl_vector_get(cs_err,i),2) + pow(gsl_vector_get(cs_rms,i),2) ) );
    }

    #pragma omp parallel for \
    default(none) shared(L,x,mix_cs, alpha_zero)\
    private (i_ind, j_ind, cs_mix, deltamix) \
    num_threads(nprocs) \
    schedule(dynamic,16) \
    reduction(+:fit_cs_mix)

    for(i_ind = 0; i_ind < L; i_ind++) {
  	    for ( j_ind = i_ind; j_ind < L; j_ind++) {
		    cs_mix = mix_cs[L*i_ind+j_ind];
            deltamix = (i_ind!=j_ind) ? -2*x->alphas[i_ind]*x->alphas[j_ind] :
            x->alphas[i_ind]*(alpha_zero - x->alphas[i_ind]);
  		    fit_cs_mix += deltamix * cs_mix;
       }
    }
    fit_cs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
    Lfunc+=0.5*(fit_cs+fit_cs_mix);
    //gsl_vector_free(cs_exp);
    //gsl_vector_free(cs_err);
    //gsl_vector_free(cs_rms);
    //gsl_matrix_free(cs_pre);
    //free(mix_cs);
  }

  //gsl_vector_free(alpha_pre);

  return Lfunc;
}

double L_distance(void *xp, void *yp)
{
  block *x = (block *) xp;
  block *y = (block *) yp;
  double vector_distance = 0.0;
  for (int i=0; i<x->size; i++) {
	vector_distance+=fabs(x->alphas[i]-y->alphas[i]);
  }
  if (x->rosettaPrior) {
    vector_distance+=fabs(x->shiftEnergy-y->shiftEnergy);
  }
  return sqrt(vector_distance);
}

//No printing is done by default
void L_print (void *xp)
{
  block *x = (block *) xp;
  double alpha_zero = 0.0;
  double weight; 
  for(int i=0; i < x->size; i++){
  	alpha_zero += x->alphas[i];
  }
  for(int i=0; i < x->size; i++){
      weight =  x->alphas[i]/alpha_zero;
     //Add vector save here
  }
}

void L_take_step(const gsl_rng * r, void *xp, double step_size)
{
    block * x = (block *) xp;
  	//The index of which alpha should be modified
  	int i = (int) round(gsl_rng_uniform(r)*x->size);

    //TODO: Rosetta prior can be a part of structure block and therefore can be conditioned from here
    //TODO: The sampled energy pushes a lot towards priors
	double s = x->shiftEnergy+gsl_ran_gaussian_ziggurat(r,step_size);
	//x->shiftEnergy = GSL_MAX(0.001, s);
	x->shiftEnergy = s;

    double u = x->alphas[i]+gsl_ran_gaussian_ziggurat(r,step_size);
	x->alphas[i] = GSL_MAX(0.001, u);
}


double ModelEvidenceEnergy(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs,
                double saxs_scale, int N)
{

	double exp_scale = 1/10e8;
	double norm_term = 1.0; 

	double fit_saxs = 0.0;

   	for( int i = 0; i< N; i++) { 
	
	norm_term *= gsl_vector_get(err_saxs,i)*2.5066282746310002;

	fit_saxs +=
	(pow( saxs_scale*gsl_vector_get(saxs_ens,i) - gsl_vector_get(saxs_exp,i),2)/
	pow(gsl_vector_get(err_saxs,i),2)); 

	}
	return exp(-fit_saxs*exp_scale);
}


double mc_integrate(gsl_matrix *saxs_pre, gsl_vector *saxs_exp,
                    gsl_vector *err_saxs, int k, int N) {

  double energy_final;
  double *alphas;
  double *samples;
  double *alpha_ens;
  int Ntrials = 1000000;

  alphas = (double * ) malloc( k * sizeof( double ));
  samples = (double * ) malloc( k * sizeof( double ));
  //alpha_ens = (double * ) malloc( k * sizeof( double ));
  const gsl_rng_type *Krng;
	gsl_rng *r;
	gsl_rng_env_setup();
	Krng = gsl_rng_default;
	r = gsl_rng_alloc(Krng);
    gsl_rng_set(r,time(NULL));

  gsl_vector *weights = gsl_vector_alloc(k);
  gsl_vector *saxs_ens = gsl_vector_alloc(N);
  for (int i = 0; i<k; i++) alphas[i] = 0.5;
  gsl_ran_dirichlet(r, k, alphas, samples);
  //for( int j = 0; j< k; j++) {
  //}
  double energy_trial=0.0;
  double saxs_scale = 0.0;
  for (int i=0; i<Ntrials; i++) {
    gsl_ran_dirichlet(r, k, alphas, samples);
    for( int j = 0; j< k; j++) {
	    gsl_vector_set(weights, j, samples[j]);
    }
    gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, weights, 0.0, saxs_ens);
    saxs_scale = SaxsScaleMean(saxs_ens,saxs_exp,err_saxs,N);
    energy_trial+=ModelEvidenceEnergy(saxs_ens,saxs_exp,err_saxs,saxs_scale,N);
  }
  cout<<"Energy final "<<energy_trial<<std::endl;
  energy_final=energy_trial/Ntrials;
  gsl_rng_free (r);
  return energy_final;

}

void calculate_alpha_priors(gsl_vector* rosetta_engeries,
                            gsl_vector* alpha_priors,  int L ) {

    double energy_sum = 0.0;
    double kBT = 0.0019872041*293;

    for( int j = 0; j< L; j++) {
        energy_sum += exp(-gsl_vector_get(rosetta_engeries,j)/kBT);
    }

    for( int j = 0; j< L; j++) {
        //double alpha = (exp(-(c_ref + gsl_vector_get(rosetta_engeries,j))/kBT));
        double alpha = exp(-(gsl_vector_get(rosetta_engeries,j))/kBT);
        gsl_vector_set(alpha_priors,j,alpha);
    }
}

double calculate_chi2( gsl_vector* saxs_ens_current, double saxs_scale_current, gsl_vector* saxs_exp, gsl_vector* err_saxs,  int N) {
    double chi2 = 0.0;
	for (int i=0; i < N; i++) {
	    chi2+=pow((gsl_vector_get(saxs_exp, i) - saxs_scale_current*gsl_vector_get(saxs_ens_current, i)),2)
	    / pow(gsl_vector_get(err_saxs, i),2);
	}
	return chi2;
}

double calculate_chi2_crysol( gsl_vector* saxs_ens_current, double saxs_scale_current, gsl_vector* saxs_exp, gsl_vector* err_saxs,  int N) {
    double chi2 = 0.0;
    double saxs_scale_crysol, mixed_term, square_calc;
    double err;
    for (int i=0; i < N; i++) {
        err = pow(gsl_vector_get(err_saxs, i),2);
        mixed_term += gsl_vector_get(saxs_exp, i)*gsl_vector_get(saxs_ens_current, i)/err;
        square_calc += gsl_vector_get(saxs_ens_current, i)*gsl_vector_get(saxs_ens_current, i)/err;
    }
    saxs_scale_crysol = mixed_term/square_calc;

	for (int i=0; i < N; i++) {
	    chi2+=pow((gsl_vector_get(saxs_exp, i) - saxs_scale_crysol*gsl_vector_get(saxs_ens_current, i)),2)
	    / pow(gsl_vector_get(err_saxs, i),2);
	}
	return chi2/N;
}
/*Overall algorithm
1. Read experimental data and parameter priors
2. Run simulated anealing to minimize function 
3. Iteratively remove structures with weights lower than wcut
*/
void run_vbw(const int &again, const int &k, const std::string &pre_weight_file,
        const std::string &structure_energy_file, const int &N, const int &n,
        const std::string &presaxsfile, const int &Ncurves,
        const std::string &curvesfile, const std::string &outfile,
        const int &nprocs, const double &w_cut, const int &skip_vbw,
        const std::string &precsfile, const std::string &rmscsfile,
        const std::string &chemical_shifts_file)
{
	//////////////////// Init section /////////////////////////////////////
	double saxs_scale_current;
	double wdelta = 0.0001;	

	gsl_siman_params_t params;
	int N_TRIES; //Seems to be inactive?
    int ITERS_FIXED_T;
    double STEP_SIZE;
    double K;
    double T_INITIAL;
    double MU_T;
    double T_MIN;

    //TODO: Samples, set to maximum 5000, which is also the maximum number of iterations.
	int samples = 500;

    ofstream output(outfile,  std::ofstream::out | std::ofstream::trunc);

    //Number of models in single iteration
    int L = k;
	double alpha_zero;
	double energy_current, energy_min;
	double *saxs_mix; 
	double *cs_mix;
	float acceptance_rate = 1.0;
    bool rosettaPrior = false;
    bool chemicalShiftsOn = false;
    bool saxsOn = true;
 	if (std::strcmp(chemical_shifts_file.c_str(),"None")!=0) {
        chemicalShiftsOn = true;
    }
    if (std::strcmp(structure_energy_file.c_str(),"None")!=0) {
	    rosettaPrior = true;
	}
 	saxs_mix = (double * ) malloc( k * k * sizeof( double ));

	gsl_matrix *saxs_pre = gsl_matrix_alloc(N,k);
	gsl_matrix *saxs_file_matrix = gsl_matrix_alloc(N,3);

	gsl_vector *saxs_exp = gsl_vector_alloc(N),
		*err_saxs = gsl_vector_alloc(N),
		*saxs_qvector = gsl_vector_alloc(N),
		*alpha_pre = gsl_vector_alloc(k),
		*rosetta_engeries = gsl_vector_alloc(k),
		*w_pre = gsl_vector_alloc(k),
		*w_ens_current = gsl_vector_alloc(k),
		*alpha_ens_current = gsl_vector_alloc(k),
		*tostart = gsl_vector_alloc(k+2),
		*saxs_ens_current = gsl_vector_alloc(N),
		*memory = gsl_vector_alloc(k+4),
		*bayesian_weight1 = gsl_vector_alloc(k),
        *bayesian_weight1_current = gsl_vector_alloc(k),
        *cs_exp = gsl_vector_alloc(n),
        *cs_err = gsl_vector_alloc(n),
		*cs_rms = gsl_vector_alloc(n),
		*cs_ens_current = gsl_vector_alloc(n);

    cs_mix = (double * ) malloc( k * k * sizeof( double ));
    gsl_matrix *cs_pre = gsl_matrix_alloc(n,k);
    gsl_matrix *cs_file_matrix = gsl_matrix_alloc(n,2);
	
	gsl_vector_set_zero(bayesian_weight1);
	gsl_matrix *weight_samples = gsl_matrix_alloc(samples,k);;

	//Marks indexes that don't pass threshold filter
	bool removed_indexes[k];
	for (int i = 0; i < k; i++) removed_indexes[i]=false;

	// Read in data from files //
	FILE * inFile = fopen(presaxsfile.c_str(),"r"); gsl_matrix_fscanf(inFile,saxs_pre);fclose(inFile);
	//Read prior weights
	inFile = fopen(pre_weight_file.c_str(),"r"); gsl_vector_fscanf(inFile,w_pre); fclose(inFile);
	//Read prior alphas
	if (rosettaPrior) {
            std::cout<<"Loading structural priors"<<std::endl;
	    inFile = fopen(structure_energy_file.c_str(),"r");
	    gsl_vector_fscanf(inFile,rosetta_engeries);
	    fclose(inFile);
	}

    //Loading chemical shift files
    if (chemicalShiftsOn) {
	std::cout<<"Loading chemical shifts"<<std::endl;
        inFile = fopen(precsfile.c_str(),"r");
        gsl_matrix_fscanf(inFile,cs_pre); fclose(inFile);

        inFile = fopen(rmscsfile.c_str(),"r");
        gsl_vector_fscanf(inFile,cs_rms); fclose(inFile);

	    FILE *inFile = fopen(chemical_shifts_file.c_str(),"r");
            gsl_matrix_fscanf(inFile,cs_file_matrix);

	    for (int i = 0;  i< n; i++) {
            gsl_vector_set(cs_exp,i,gsl_matrix_get(cs_file_matrix,i,0));
       	    gsl_vector_set(cs_err,i,gsl_matrix_get(cs_file_matrix,i,1));
        }
	    fclose(inFile);
	}

	//Read scattering file
    FILE *inSAXSdat = fopen(curvesfile.c_str(),"r");
    gsl_matrix_fscanf(inSAXSdat,saxs_file_matrix);
    for (int i = 0;  i< N; i++) {
        gsl_vector_set(saxs_qvector,i,gsl_matrix_get(saxs_file_matrix,i,0));
        gsl_vector_set(saxs_exp,i,gsl_matrix_get(saxs_file_matrix,i,1));
       	gsl_vector_set(err_saxs,i,gsl_matrix_get(saxs_file_matrix,i,2));
    }
    fclose(inSAXSdat);
	
	cout<<"Files reading finished"<<std::endl;
	// initialize random number generators //
	const gsl_rng_type *Krng; 
	gsl_rng *r; 
	gsl_rng_env_setup(); 
	Krng = gsl_rng_default;
	r = gsl_rng_alloc(Krng); 
	gsl_rng_set(r,time(NULL)); 

    //Setting initial reference energy equal to number_of_models/2
	double energy_sum = 0.0;
    double kBT = 0.0019872041*293;
    double shiftEnergyInternal;

    //Jenso-Shannon Divergance
    double jsd1_sum = 0.0;

    //Skipping VBW and going directly to model evidence integration
	if (!skip_vbw) {

    //Setting priors based on energies
	if (rosettaPrior) {
        for( int j = 0; j< L; j++) {
            energy_sum += exp(-gsl_vector_get(rosetta_engeries,j)/kBT);
        }
	    shiftEnergyInternal = -kBT*log(0.5*L/energy_sum);
    }
    calculate_alpha_priors(rosetta_engeries, alpha_pre, k);

    block *simAnBlock = block_alloc(k);
	//Initiallize alphas based on rosetta energies
	simAnBlock->shiftEnergy = shiftEnergyInternal;
	simAnBlock->saxsExpPtr = saxs_exp;
	simAnBlock->saxsErrPtr = err_saxs;
	simAnBlock->saxsPrePtr = saxs_pre;
	simAnBlock->alphaPre = alpha_pre;
	simAnBlock->csExpPtr = cs_exp;
    simAnBlock->csErrPtr = cs_err;
	simAnBlock->csRmsPtr = cs_rms;
    simAnBlock->csPrePtr = cs_pre;
	simAnBlock->numberProcs = nprocs;
	simAnBlock->rosettaPrior = rosettaPrior;
    simAnBlock->chemicalShiftsOn = chemicalShiftsOn;
	simAnBlock->saxsOn = saxsOn;

    //Initialize alphas with flat prior values
	for (int i = 0; i < k; i++) {
        simAnBlock->alphas[i] = gsl_vector_get(w_pre,i);
	}

	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_pre, 0.0, saxs_ens_current);
	if (chemicalShiftsOn) {
	    gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_pre, 0.0, cs_ens_current);
    }
	//TODO: Scaling is not required here?
	saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
	simAnBlock->saxsScale = saxs_scale_current;
	simAnBlock->saxsEnsPtr = saxs_ens_current;
	//simAnBlock->csEnsPtr = cs_ens_current;
	if(again == 1){ inFile = fopen("restart.dat","r"); gsl_vector_fscanf(inFile,tostart); fclose(inFile); }
	//timeval t1, t2;
	//double elapsedTime;
  	//gettimeofday(&t1, NULL);

	double smix;
	double csmix;
    #pragma omp parallel for reduction(+:smix) reduction(+:csmix) num_threads(nprocs)
	//#pragma omp parallel for reduction(+:smix) num_threads(nprocs)
	for( int i = 0; i< k; i++) {
        for (int j = 0; j < k; j++) {
		    smix = 0.0;
		    csmix = 0.0;
            for (int m = 0; m < N; m++) {
			    smix+=gsl_matrix_get(saxs_pre,m,i)*gsl_matrix_get(saxs_pre,m,j)
			    /pow(gsl_vector_get(err_saxs,m),2);
		    }
		    if (chemicalShiftsOn) {
		        for (int m = 0; m < n; m++) {
                    csmix+=gsl_matrix_get(cs_pre,m,i)*gsl_matrix_get(cs_pre,m,j)
                    /(pow(gsl_vector_get(cs_err,m),2)+pow(gsl_vector_get(cs_rms,m),2));
                }
                cs_mix[i*k+j] = csmix;
            }
            saxs_mix[i*k+j] = smix;
        }
	}
	/*gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  	elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  	cout << "Time: "<< elapsedTime << " ms."<<std::endl;*/

	simAnBlock->saxsMixPtr = saxs_mix;
	simAnBlock->csMixPtr = cs_mix;
	///////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////
	if(again == 1) { 
		inFile = fopen("restart.dat","r"); 
		gsl_vector_fscanf(inFile,tostart); 
		fclose(inFile); 
		for( int i = 0; i< k; i++) gsl_vector_set(alpha_ens_current,i,gsl_vector_get(tostart,i));
		energy_min = gsl_vector_get(tostart,k);
		simAnBlock->saxsScale = gsl_vector_get(tostart,k+1);
	}
	else {	
		////////////////////// First iteration ////////////////////////////////
		cout<<"Equilibration starting..."<<std::endl;
		
		N_TRIES = 1; //Seems to be inactive?
		ITERS_FIXED_T = 1;
		STEP_SIZE = 1;
		K = 1.0;
		T_INITIAL = 2.0; 
		MU_T = 1.000025;
       	T_MIN = 2.7776e-11;
		params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//Define params before equilibration and after for next rounds
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
		 	block_copy, block_copy_construct, block_destroy,                
            0, params, &acceptance_rate);

		alpha_zero = 0.0;
		for (int i=0; i < k; i++) {
			alpha_zero+=simAnBlock->alphas[i];
			gsl_vector_set(alpha_ens_current,i,simAnBlock->alphas[i]);
		}
		for (int i=0; i < k; i++) {
            gsl_vector_set(w_ens_current,i,gsl_vector_get(alpha_ens_current,i)/alpha_zero);
        }
		energy_min = L_function(simAnBlock);
		if (rosettaPrior) {
    	    shiftEnergyInternal = simAnBlock->shiftEnergy;
		    cout<<"Shift energy: "<<shiftEnergyInternal<<std::endl;
	    }
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current, 0.0, saxs_ens_current);
		saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
        if (chemicalShiftsOn) {
		    gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_ens_current, 0.0, cs_ens_current);
		}
		block_destroy(simAnBlock);
		free(saxs_mix); //TODO: CHeck if these doesn't corrupt memory
		//free(cs_mix); //TODO: CHeck if these doesn't corrupt memory
		/////////////////////////////////////////////////////////////////////
	
		//Store alphas after equilibration stage
		ofstream restart("restart.dat");
        	for(int j = 0; j < k; j++) { restart << gsl_vector_get(alpha_ens_current,j)<<" "; }
        	restart <<energy_min<<" "<<saxs_scale_current<<std::endl;
        	restart.close();
	}
			
	///////////////////Next iterations //////////////////////////////////
	cout<<"Simulated annealing started"<<std::endl;
	int overall_iteration = 0;
	int sampling_step;
	int last_updated;
	int l, m, newL;
    L = k;
	//Energy from first iteration
	while ( L > 1 ) {
		cout<<"Starting "<<overall_iteration+1<<" iteration with "<<L<<" models"<<std::endl;
		
		block *simAnBlock = block_alloc(L);
		gsl_matrix *saxs_pre_round = gsl_matrix_alloc(N,L);
		gsl_matrix *cs_pre_round = gsl_matrix_alloc(n,L);
    	gsl_vector *alpha_pre_round = gsl_vector_alloc(L);
		gsl_vector *rosetta_engeries_round = gsl_vector_alloc(L);
		double  *saxs_mix_round =  (double * ) malloc( L * L * sizeof( double ));
        double  *cs_mix_round =  (double * ) malloc( L * L * sizeof( double ));
		l = 0;
		for (int i = 0; i < k; i++) {
			if (removed_indexes[i]==false) {
			    //TODO: Here we need to call function that recalculates alphas_pre
				if (rosettaPrior) {
				    gsl_vector_set(rosetta_engeries_round,l,gsl_vector_get(rosetta_engeries,i));
				}
                if (chemicalShiftsOn) {
                    for (int j = 0; j < n; j++) {
					    gsl_matrix_set(cs_pre_round,j,l,gsl_matrix_get(cs_pre,j,i));
				    }
                }
				for (int j = 0; j < N; j++) {
					gsl_matrix_set(saxs_pre_round,j,l,gsl_matrix_get(saxs_pre,j,i));
				}
                simAnBlock->alphas[l] = gsl_vector_get(alpha_ens_current,i);
				l++;
			}
        }

        //The value is recalculated based on number of structures
        //Alternatively the value from the previous iteration can be taken.
        energy_sum = 0.0;

        for( int j = 0; j< L; j++) {
            energy_sum += exp(-gsl_vector_get(rosetta_engeries,j)/kBT);
        }
	    shiftEnergyInternal = -kBT*log(0.5*L/energy_sum);
        simAnBlock->shiftEnergy = shiftEnergyInternal;

        cout<<"Shift energy calculated"<<std::endl;

		#pragma omp parallel for reduction(+:smix) reduction(+:csmix) num_threads(nprocs)
		//#pragma omp parallel for reduction(+:smix) num_threads(nprocs)
        for( int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
		         smix = 0.0;
		         csmix = 0.0;
                 for (int m = 0; m < N; m++) {
		            smix+=gsl_matrix_get(saxs_pre_round,m,i)
		            *gsl_matrix_get(saxs_pre_round,m,j)/pow(gsl_vector_get(err_saxs,m),2);
                 }
                 if (chemicalShiftsOn) {
                     for (int m = 0; m < n; m++) {
                        csmix+=gsl_matrix_get(cs_pre_round,m,i)*
                        gsl_matrix_get(cs_pre_round,m,j)/
                        (pow(gsl_vector_get(cs_err,m),2)+pow(gsl_vector_get(cs_rms,m),2));
                    }
                    cs_mix_round[i*L+j]=csmix;
                 }
	        saxs_mix_round[i*L+j]=smix;
           	}
		}

        cout<<"Mix terms calcualted"<<std::endl;
        if (rosettaPrior) {
            cout<<"Alphas: ";
            calculate_alpha_priors(rosetta_engeries_round, alpha_pre_round, L);
            for ( int i = 0; i < L; i++) {
                const double beta = 1.717472947;
                double shiftEnergyExp = exp(-simAnBlock->shiftEnergy*beta);
                cout<<shiftEnergyExp*gsl_vector_get(alpha_pre_round,i)<<" ";
            }
            cout<<std::endl;
            cout<<"Shift energy: "<<shiftEnergyInternal<<std::endl;
        }
		//saxs_exp and err_saxs are independent of run
        simAnBlock->saxsExpPtr = saxs_exp;
        simAnBlock->saxsErrPtr = err_saxs;
		simAnBlock->saxsPrePtr = saxs_pre_round;
        simAnBlock->saxsMixPtr = saxs_mix_round;
		simAnBlock->saxsEnsPtr = saxs_ens_current;
		simAnBlock->saxsScale = saxs_scale_current;

        simAnBlock->csExpPtr = cs_exp;
        simAnBlock->csErrPtr = cs_err;
		simAnBlock->csRmsPtr = cs_rms;
        simAnBlock->csPrePtr = cs_pre_round;
        simAnBlock->csMixPtr = cs_mix_round;
        //simAnBlock->csEnsPtr = cs_ens_current;

        simAnBlock->alphaPre = alpha_pre_round;
        simAnBlock->rosettaPrior = rosettaPrior;
		simAnBlock->chemicalShiftsOn = chemicalShiftsOn;
		simAnBlock->saxsOn = saxsOn;
		simAnBlock->numberProcs = nprocs;
		
		////////////////////////Short equilibration period to find step size/////////////////////////
		/*N_TRIES = 1;
        ITERS_FIXED_T = 1000;
        K = 1.0;
        T_INITIAL = 1.0;
        MU_T = 1.00005;
        T_MIN = 1.0;
		//Itertate over different step size
		float dmin = 10;
		cout<<"Initial step estimation"<<std::endl;
		for (double s=0.01; s<2.1; s+=0.1) { 
            params = {N_TRIES, ITERS_FIXED_T, s, K, T_INITIAL, MU_T, T_MIN};
		
             //alphas are used from the previous simulation
             gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
                block_copy, block_copy_construct, block_destroy,
                0, params, &acceptance_rate);
			//cout<<"Acceptance Rate "<<acceptance_rate<<std::endl;
			if(fabs(acceptance_rate -0.5) < dmin) { 
				dmin = fabs(acceptance_rate -0.5);
				STEP_SIZE = s;		
			}	
		}
		cout<<"STEP_SIZE set to: "<<STEP_SIZE<<std::endl;*/
		///////////////////////////////////////////////////////////////////////////////////////////
		N_TRIES = 1.0;
        ITERS_FIXED_T = 1.0;
        STEP_SIZE = 1.0;
        K = 1.0;
        T_INITIAL = 1.0;
		MU_T = 1.00005;
		T_MIN = 1.3888e-11;
        params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//alphas are used from the previous simulation 
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
                  		block_copy, block_copy_construct, block_destroy, 
                  		0, params, &acceptance_rate);

		energy_current = L_function(simAnBlock);

        //If L_function doesn't improve after 10 iterations exit program
		newL = 0;
		m = 0; 
		alpha_zero = 0.0;
		for ( int i = 0; i < L; i++ ) alpha_zero +=simAnBlock->alphas[i];

		double new_alpha_zero = 0.0;
		for ( int i = 0; i < k; i++ ) {
			if ( removed_indexes[i]==false ) {
				double wib = simAnBlock->alphas[m]/alpha_zero;
				if ( wib < w_cut ) {
					gsl_vector_set( alpha_ens_current, i, 0.0 );
					gsl_vector_set( w_ens_current, i, 0.0);
					removed_indexes[i] = true;
				} else {
					new_alpha_zero += simAnBlock->alphas[m];
					gsl_vector_set( alpha_ens_current, i, simAnBlock->alphas[m] );
					newL++;
				}
				m++;
			}
		}
		
		//int wdelta_count = 0;
		for ( int i = 0; i < k; i++ ) {
            if (removed_indexes[i]==false) {
				gsl_vector_set( w_ens_current,i,gsl_vector_get(alpha_ens_current,i)/new_alpha_zero );
			}
		}
		//Stoping simulations if weights don't change for more than delta (0.001)
		//if (wdelta_count == newL) {cout<<"Simulations stopped because weights don't progress"<<std::endl; break;}

		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current, 0.0, saxs_ens_current);
	    saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
		if (chemicalShiftsOn) {
		    gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_ens_current, 0.0, cs_ens_current);
		}
		//Structural library size after discarding structures with weight lower than cuttof
		shiftEnergyInternal = simAnBlock->shiftEnergy;
		L = newL;
		
		block_destroy(simAnBlock);	
		overall_iteration++;

	    double chi2 = calculate_chi2_crysol(saxs_ens_current, saxs_scale_current, saxs_exp, err_saxs, N);

		sampling_step = overall_iteration-1;
		for (int jind=0; jind<k; jind++) {
            gsl_matrix_set(weight_samples,sampling_step,jind,gsl_vector_get(w_ens_current,jind));
        }

        double niter = 1.0/double(sampling_step+1);
        gsl_vector_add(bayesian_weight1,w_ens_current);
        gsl_vector_memcpy(bayesian_weight1_current,bayesian_weight1);
        gsl_vector_scale(bayesian_weight1_current,niter);


        //Calculating posterior expected divergence
        //TODO: Make a cluean-up with vector
        //double jsd1_sum = 0.0;
        double jsd1 = 0.0;
        //for (int s=0; s<sampling_step; s++) {
        //for (int j=0; j<k; j++) {
        // gsl_vector_set(bayesian_weight1,j,gsl_matrix_get(weight_samples,s,j));
        //}
        jsd1 = jensen_shannon_div(bayesian_weight1_current,w_ens_current,k);
        //jsd1_sum += sqrt(jsd1);
	jsd1_sum += jsd1;
        //}

		if (energy_current < energy_min) {
			energy_min = energy_current;
            last_updated = overall_iteration;
		} //Adding this one WOjtek
			for( int l = 0; l < k; l++) {
                gsl_vector_set(memory, l , gsl_vector_get(w_ens_current,l));
        	}

        	gsl_vector_set(memory, k, saxs_scale_current);
        	gsl_vector_set(memory, k+1, energy_current);
            gsl_vector_set(memory, k+2, chi2);
            gsl_vector_set(memory, k+3, jsd1);
        	//All weights scale factor, energy_currenr and chi2
        	for( int j = 0; j < k + 3; j++) output << gsl_vector_get(memory,j) << " ";
        	output <<gsl_vector_get(memory,k+3)<<endl;
			//output.close();
		//}


		free(saxs_mix_round);
		gsl_matrix_free(saxs_pre_round);
		free(cs_mix_round);
        gsl_matrix_free(cs_pre_round);
		gsl_vector_free(alpha_pre_round);
        gsl_vector_free(rosetta_engeries_round);

		if ((overall_iteration-last_updated)>10) {
            cout<<"Energy hasn't decreased for 10 iterations. Stopping simulations"<<std::endl;
            break;
        }
		
		if (overall_iteration == samples) {
			cout<<"Maximum number of iteration has been reached. Stopping simulation"<<std::endl;
			break;
		}

	}	
	///////////////////////////////////////////////////////////////////////	

    output<<"Jesen-Shannon Div "<<sqrt(jsd1_sum/double(sampling_step))<<std::endl;
    if (rosettaPrior) {
        output<<"BoltzmanShift "<<shiftEnergyInternal<<std::endl;
    }
    //}//Finish VBW section
    }
    //Model Evidence calculation
    if (skip_vbw == 1)  L = k;

    gsl_matrix *saxs_pre_selected = gsl_matrix_alloc(N,L);
    if (skip_vbw == 1)  {
        int l = 0;
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < N; j++) {
				gsl_matrix_set(saxs_pre_selected,j,l,gsl_matrix_get(saxs_pre,j,i));
			}
         	l++;
        }
    }
    else {
        //Copy saxs_pre_round for monte carlo integration
		int l = 0;
		for (int i = 0; i < k; i++) {
			if (removed_indexes[i]==false) {
				for (int j = 0; j < N; j++) {
					gsl_matrix_set(saxs_pre_selected,j,l,gsl_matrix_get(saxs_pre,j,i));
				}
         		l++;
			}
        }
    }

    saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
    double chi2 = calculate_chi2_crysol(saxs_ens_current, saxs_scale_current, saxs_exp, err_saxs, N);
    output<<"\nChi2 "<<chi2<<std::endl;

    //TODO: Check if saxs_ens_current keeps latest value and if needs to be scaled
    ofstream fitoutput(outfile+".fit", std::ofstream::out | std::ofstream::trunc);
    for (int i = 0; i < N; i++) {
        fitoutput.precision(8);
        fitoutput.showpoint;
	    fitoutput.width(14);
	    fitoutput<<gsl_vector_get(saxs_qvector,i);
	    fitoutput.width(14);
	    fitoutput<<gsl_vector_get(saxs_exp,i);
	    fitoutput.width(14);
	    fitoutput<<saxs_scale_current*gsl_vector_get(saxs_ens_current,i);
	    fitoutput.width(14);
	    fitoutput<<gsl_vector_get(err_saxs,i)<<std::endl;
    }

    double model_evd;
    model_evd = mc_integrate(saxs_pre_selected, saxs_exp, err_saxs, L, N);
    output<<"\nModel Evidence: "<<model_evd<<std::endl;
    output.close();
    fitoutput.close();
	gsl_rng_free (r);
}

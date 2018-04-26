stan_code = """
data {
  int n_measures;
  int n_structures;
  vector[n_measures] target_curve;
  vector[n_measures] target_errors;
  matrix[n_measures, n_structures] sim_curves;
  vector[n_structures] priors;
}

parameters {
  simplex[n_structures] weights;
  real scale;
}

model {
  vector[n_structures] alphas;
  vector[n_measures] pred_curve;
  alphas = priors;
  //scale ~ uniform(0,2);
  scale ~ lognormal(0.8,0.4);
  weights ~ dirichlet(alphas);
  pred_curve = sim_curves * weights * scale;
  target_curve ~ normal(pred_curve, target_errors);
}

"""
stan_code_CS = """
data {
  int n_measures;
  int m_measures;
  int n_structures;
  vector[n_measures] target_saxs;
  vector[n_measures] target_saxserr;
  matrix[n_measures, n_structures] sim_saxs;
  vector[m_measures] target_cs;
  vector[m_measures] target_cserr;
  vector[m_measures] sim_cserr;
  matrix[m_measures, n_structures] sim_css;
  vector[n_structures] priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0> scale;
}

model {
  vector[n_measures] pred_saxs;
  vector[m_measures] pred_css;
  vector[m_measures] cs_errors;
  vector[n_structures] alphas;
  scale ~ uniform(0,100);
  alphas = priors;
  weights ~ dirichlet(alphas);
  pred_saxs = sim_saxs * weights * scale;
  pred_css = sim_css * weights;
  target_saxs ~ normal(pred_saxs, target_saxserr);
  for (i in 1:m_measures) {
    cs_errors[i] = pow(pow(sim_cserr[i],2)+pow(target_cserr[i],2),0.5);
  }
  target_cs ~ normal(pred_css, cs_errors);
}
"""

stan_code_EP = """
data {
  int n_measures;
  int n_structures;
  vector[n_measures] target_curve;
  vector[n_measures] target_errors;
  matrix[n_measures, n_structures] sim_curves;
  vector[n_structures] energy_priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0> scale;
  real<lower=0> boltzmann_shift;
}

model {
  vector[n_measures] pred_curve;
  vector[n_structures] alphas;
  scale ~ uniform(0,100);
  boltzmann_shift ~ uniform(0,300);
  alphas = exp(-1.717472947*(boltzmann_shift+energy_priors));
  weights ~ dirichlet(alphas);
  pred_curve = sim_curves * weights * scale;
  target_curve ~ normal(pred_curve, target_errors);
}
"""

stan_code_EP_CS = """
data {
  int n_measures;
  int m_measures;
  int n_structures;
  vector[n_measures] target_saxs;
  vector[n_measures] target_saxserr;
  matrix[n_measures, n_structures] sim_saxs;
  vector[m_measures] target_cs;
  vector[m_measures] target_cserr;
  vector[m_measures] sim_cserr;
  matrix[m_measures, n_structures] sim_css;
  vector[n_structures] energy_priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0> scale;
  real<lower=0> boltzmann_shift;
}

model {
  vector[n_measures] pred_saxs;
  vector[m_measures] pred_css;
  vector[m_measures] cs_errors;
  vector[n_structures] alphas;
  boltzmann_shift ~ uniform(0,300);
  scale ~ uniform(0,100);
  alphas = exp(-1.717472947*(boltzmann_shift+energy_priors));
  weights ~ dirichlet(alphas);
  pred_saxs = sim_saxs * weights * scale;
  pred_css = sim_css * weights;
  for (i in 1:m_measures) {
    cs_errors[i] = pow(pow(sim_cserr[i],2)+pow(target_cserr[i],2),0.5);
  }
  target_cs ~ normal(pred_css, cs_errors);
  target_saxs ~ normal(pred_saxs, target_saxserr);
}
"""
posterior_predict_quanities = """
generated quantities {
    vector[n_measures] pred_curve;
    vector[n_measures] target_curve_tilde;
    pred_curve = sim_curves * weights * scale;
    for (i in 1:n_measures) {
        target_curve_tilde[i] = cauchy_rng(pred_curve[i], target_errors[i]);
    }
}
"""

psisloo_quanities = """
generated quantities {
    real loglikes[n_measures];
    vector[n_measures] pred_curve;
    pred_curve = sim_curves * weights * scale;
    for (i in 1:n_measures) {
            loglikes[i] = normal_lpdf(target_curve[i] | pred_curve[i], target_errors[i]);
    }
}
"""
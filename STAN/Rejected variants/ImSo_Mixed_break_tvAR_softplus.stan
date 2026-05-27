functions {
  #include "backshift_matrix.stan"
}

data {
  int<lower=0> N;
  matrix[N,21] y;
  int<lower=0> max_shift;
}

transformed data {
  vector[21] freq;
  row_vector[21] log_freq;
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
  log_freq = log(freq)';
  
  int numb_shift_coeff = 21*max_shift-max_shift*(max_shift+1)/2;
  
  
  // Segments: 1 breaks => 2 segements => 3 betas
  int n_para;
  int n_beta;
  int first_bp;
  
  n_beta = 3;
  n_para = 4;
  first_bp = n_beta + 1;
  
  // prioris
  vector[n_beta] mean_beta;
  vector[n_beta] sd_beta;
  
  mean_beta = [76.61, -6.55, -11.7]';
  sd_beta = [1, 0.2, 0.6]';
  
}

parameters {
  // Fix
  vector[n_beta] fe_betas;
  real<lower=50,upper=5000> fe_bp1;
  // RE Hyperpriors
  //vector<lower=0>[n_para] hyper_v;
  
  // tvAR coeffs
  vector[numb_shift_coeff] shift_coeff;
  vector<lower=0>[21] variances;
  
  // RE
  cholesky_factor_corr[n_para] REcorr_chol;
  vector<lower=0>[n_para] re_var;
  real<lower=0> re_c_concentration;
  matrix[N,n_para] z;
}
transformed parameters {
  // Fix
  vector[n_para] fe_phis;
  fe_phis[1:n_beta] = fe_betas;
  fe_phis[first_bp] = fe_bp1;
  
  // Fix + Random
  matrix[N,n_para] phis;
  matrix[N,n_beta] betas;
  vector<lower=0>[N] bp1;

  // Data (Co)variance
  matrix[21,21] backshift;
  vector<lower=0>[21] st_dev;
  cholesky_factor_cov[21] st_dev_mat;
  cholesky_factor_cov[21] difference_matrix;
  cholesky_factor_cov[21] cov_mat_chol;
  
  
  backshift = backshift_matrix(21,shift_coeff);

  
  st_dev = sqrt(variances);
  st_dev_mat = diag_matrix(st_dev);
  
  difference_matrix = diag_matrix(rep_vector(1.0, 21)) - backshift; // idenität nach trafo data schieben wenn es klappt
  cov_mat_chol = mdivide_left_tri_low(difference_matrix,st_dev_mat);
  
  // RE
  matrix[N,n_para] re;
  cholesky_factor_cov[n_para] REcov_mat_chol;
  vector<lower=0>[n_para] re_sd = sqrt(re_var);
  
  REcov_mat_chol = diag_pre_multiply(re_sd, REcorr_chol);
  re = (REcov_mat_chol * z')';
  
  
  // Expected value
  matrix[N,21] mu;
  
  phis = rep_matrix(to_row_vector(fe_phis), N) + re;
  betas = phis[,1:n_beta];
  bp1 = phis[,first_bp];
  
  vector[N] off12 = -log(bp1);
   
  matrix[N, 21] intercepts = betas[,1] * rep_row_vector(1.0, 21);
  matrix[N, 21] bases = betas[,2] * log_freq;
  matrix[N,21] addon = rep_matrix(0,N,21);
  mu = intercepts + bases;
   
   
  // beta3 * log(1+exp[10*(log(f) - bp)])/10
  matrix[N,21] eta;
  vector[N] beta_col= betas[, n_beta-1];// Spalte extrahieren

  
  // Outer-Sum: N × 21
  eta =  rep_matrix(off12, 21) + rep_matrix(log_freq, N);   // <-- OHNE Transpose

  
  // elementweise softplus mit Skalierung
  eta = log1p_exp(10 * eta) / 10;
  
  // zeilenweise Skalierung
  addon = diag_pre_multiply(beta_col, eta);


  
  mu += addon;
  
  
  // Residuen
  matrix[N, 21] resid;
  matrix[N, 21] st_resid;
  
  // Berechnung der Residuen
  resid = y - mu;

  // Standardisierte Residuen
  st_resid = (mdivide_left_tri_low(cov_mat_chol, resid'))';
  
}
model {
  
  // Fix
  fe_betas ~ multi_normal_cholesky(mean_beta, diag_matrix(sd_beta));
  fe_bp1 ~ normal(1456, 21);
  
  // RE
  to_vector(z) ~ normal(0,1);
  re_c_concentration ~ normal(0,1);
  re_var ~ normal(0,1);
  REcorr_chol ~ lkj_corr_cholesky(re_c_concentration);
  
  
  // AR coeff
  shift_coeff ~ normal(0,1);


  // Data covariance
  variances ~ exponential(0.1);
  
  // Data distrubution
  to_vector(st_resid) ~ normal(0,1);
}
generated quantities {
  // Log-Likelihood
  vector[N] log_lik;
  real log_cov_det = sum(log(diagonal(cov_mat_chol)));
  
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(st_resid[i] | 0, 1) - log_cov_det;
  }
  
  // Data Cov/Corr
  cov_matrix[21] cov_mat;
  //corr_matrix[21] corr_mat;
  cov_mat = multiply_lower_tri_self_transpose(cov_mat_chol);
  //corr_mat = multiply_lower_tri_self_transpose(corr_mat_cholesky);
  
  // RE Cov/Corr
  matrix[n_para,n_para] REcov;
  matrix[n_para,n_para] REcorr;
  
  REcov = multiply_lower_tri_self_transpose(REcov_mat_chol);
  REcorr = multiply_lower_tri_self_transpose(REcorr_chol);
  
  
  
} 


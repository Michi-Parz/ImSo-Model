functions {
  #include "backshift_matrix.stan"
}

data {
  int<lower=0> N;
  matrix[N,21] y;
  int<lower=0> max_shift;
  //matrix[21, 3] priorHelpMat;
}

transformed data {
  vector[21] freq;
  vector[21] log_freq;
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
  log_freq = log(freq);

  int numb_shift_coeff = 21*max_shift-max_shift*(max_shift+1)/2;

  // Segments: 1 break => 2 segements => 3 betas
  int n_para;
  int n_beta;
  int first_bp;

  n_beta = 3;
  n_para = 4;
  first_bp = n_beta + 1;

  // prioris
  vector[n_beta] mean_beta;
  vector[n_beta] sd_beta;

  mean_beta = [38.64, 4.88, -18.8]';
}

parameters {
  // Fix
  vector[n_beta] fe_betas;
  real<lower=0> b_point_diff;
  // tvAR coeffs
  vector[numb_shift_coeff] shift_coeff;
  vector<lower=0>[21] variances;
}

transformed parameters {
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

  // Breakpoints
  real b_point;

  b_point = 50 + b_point_diff;
  // Expected value
  vector[21] mu;
  
  matrix[21,n_beta] designMat;
  
  designMat[,1] = rep_vector(1.0, 21);
  designMat[,2] = log_freq;
  designMat[,3] = fmax(log_freq - log(b_point), 0);

  // vector[21] intercept = fe_betas[1]*rep_vector(1.0, 21);
  // vector[21] base = fe_betas[2]*log_freq;
  // 
  // vector[21] hinge1 = fe_betas[3]*fmax(log_freq - log(b_point), 0);
  // 
  // mu = intercept + base + hinge1;
  
  mu = designMat * fe_betas;

  
  // Precision
  vector<lower=0>[21] InvStDev = inv(st_dev);
  matrix[21,21] InvStDevMat = diag_matrix(InvStDev);
  
  matrix[21,21] precisionMatrix;
  matrix[21,21] precisionMatrixChol = InvStDevMat*difference_matrix;
  
  precisionMatrix = precisionMatrixChol'*precisionMatrixChol;
  
  
  // beta prior precision matrix
  matrix[n_beta,n_beta] betaPriorPrecision;
  betaPriorPrecision = designMat' * precisionMatrix * designMat;


}

model {

  // Coeff
  fe_betas ~ multi_normal_prec(mean_beta, betaPriorPrecision);

  # Der Bruchpunkt liegt bei 50+b_point_diff
  # A-Priori soll er bei 86 sein.
  # Deshalb 86-50=36
  b_point_diff ~ exponential(1/36.0);

  // AR coeff
  shift_coeff ~ normal(0,1);

  // Data covariance
  variances ~ exponential(1/900.0);
  // Data distrubution
  //to_vector(st_resid) ~ normal(0,1);

  for (i in 1:N) {
    y[i] ~ multi_normal_cholesky(mu, cov_mat_chol);
  }
}
generated quantities {

  // Residuen
  matrix[N, 21] resid;
  matrix[N, 21] st_resid;

  // Berechnung der Residuen
  resid = y - rep_matrix(mu', N);
  // Standardisierte Residuen
  st_resid = (mdivide_left_tri_low(cov_mat_chol, resid'))';

  // Log-Likelihood
  vector[N] log_lik;
  real log_cov_det = sum(log(diagonal(cov_mat_chol)));

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(st_resid[i] | 0, 1) - log_cov_det;
  }
  
  
  // Data cov
  matrix[21,21] cov_mat;
  cov_mat = multiply_lower_tri_self_transpose(cov_mat_chol);

   // Data corr
  matrix[21,21] corr_mat_cholesky;
  matrix[21,21] corr_mat;
  vector[21] unconditional_st_dev;
  unconditional_st_dev = sqrt(diagonal(cov_mat));
  corr_mat_cholesky = diag_pre_multiply(1 ./ unconditional_st_dev, cov_mat_chol);
  corr_mat = corr_mat_cholesky * corr_mat_cholesky';

} 


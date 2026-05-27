data {
  int<lower=0> N;
  matrix[N,21] y;
}
transformed data {
  vector[21] freq;
  vector[21] log_freq;
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
  for (j in 1:21) {
    log_freq[j] = log(freq[j]);
  }
  
  // prioris
  vector[3] mean_beta;
  vector[3] sd_beta;
  
  mean_beta = [73.6,-6,-13]';
  sd_beta = [0.8,0.2,0.6]';
  
  
}
parameters {
  // Fix
  vector[3] fe_betas;
  real<lower=50,upper=5000> fe_bp1;
  // RE Hyperpriors
  vector<lower=0>[4] hyper_v;
  // Kovarianz
  real<lower=0> c_concentration;
  vector<lower=0>[21] variances;
  cholesky_factor_corr[21] corr_mat_cholesky;
  // RE
  cholesky_factor_corr[4] REcorr_chol;
  vector<lower=0>[4] re_var;
  real<lower=0> re_c_concentration;
  matrix[N,4] z;
}
transformed parameters {
  // Fix
  vector[4] fe_phis;
  fe_phis[1:3] = fe_betas;
  fe_phis[4] = fe_bp1;
  
  // Fix ü Random
  matrix[N,4] phis;
  matrix[N,3] betas;
  vector<lower=0>[N] bp1;

  // Data variance
  vector<lower=0>[21] st_dev;
  cov_matrix[21] cov_mat;
  corr_matrix[21] corr_mat;
  cholesky_factor_cov[21] cov_mat_chol;

  // Cov Mat
  st_dev = sqrt(variances);

  cov_mat_chol = diag_matrix(st_dev)*corr_mat_cholesky;

  cov_mat = cov_mat_chol*cov_mat_chol';
  corr_mat = corr_mat_cholesky * corr_mat_cholesky';
  
  // RE
  matrix[N,4] re;

  cholesky_factor_cov[4] REcov_mat_chol;

  REcov_mat_chol = diag_matrix(sqrt(re_var))*REcorr_chol;
  
  re = z * REcov_mat_chol';
  
  cov_matrix[4] REcov;
  corr_matrix[4] REcorr;
  
  REcov = REcov_mat_chol*REcov_mat_chol';
  REcorr = REcorr_chol*REcorr_chol';
  
  // Expected value
  matrix[N,21] mu;
  
   for (i in 1:N) {
     phis[i] = fe_phis' + re[i];
      
      betas[i] = phis[i,1:3];
      bp1[i] = phis[i,4];
   }
   
   vector[N] off12 = -log(bp1);
   
   for (i in 1:N) {
  
    for (j in 1:21) {
      mu[i,j] = betas[i,1] + betas[i,2] * log_freq[j];
  
      if (freq[j] > bp1[i]) {
        mu[i,j] += betas[i,3] * (log_freq[j]+off12[i]);
        }
    }
  }
}
model {
  
  // Fix
  fe_bp1 ~ normal(1500, 59);
  fe_betas ~ multi_normal_cholesky(mean_beta, diag_matrix(sd_beta));
  
  // RE
  to_vector(z) ~ normal(0,1);
  hyper_v ~ gamma(2, 0.1);
  
  for (i in 1:4) {
    re_var[i] ~ exponential(1.0 / hyper_v[i]);
  }
  re_c_concentration ~ exponential(1);
  REcorr_chol ~  lkj_corr_cholesky(re_c_concentration);


  // Data covariance
  for (j in 1:21) {
    variances[j] ~ exponential(0.1);
  }
  c_concentration ~ exponential(1);
  corr_mat_cholesky ~ lkj_corr_cholesky(c_concentration);
  
  // Data distrubution
  for (i in 1:N) {
    y[i] ~ multi_normal_cholesky(mu[i], cov_mat_chol);
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = multi_normal_cholesky_lpdf(y[i] | mu[i], cov_mat_chol);
  }
  
  matrix[21, 21] CC_inv; // Inverse der Kovarianzmatrix
  matrix[N, 21] resid;
  matrix[N, 21] st_resid;
  
  // Berechnung der Residuen
  resid = y - mu;


  // Berechne die Inverse der Cholesky Kovarianzmatrix:
  CC_inv = inverse(cov_mat_chol);

  // Standardisierte Residuen
  st_resid = resid * CC_inv';
  
} 


data {
  int<lower=0> N;
  matrix[N,21] y;
}

transformed data {
  vector[21] freq;
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
}


parameters {
  // Quasi-Physical
  real intercept;
  //real slope1;
  real slope;//2;
  //real breakpoint;
  
  // Covariance
  //real<lower=0> c_concentration;
  real<lower=0> variances;
  //cholesky_factor_corr[21] corr_mat_cholesky;
}


transformed parameters {
  vector[21] mu;
  matrix[21,21] Sigma;
  
  
  // for (j in 1:21){
  //   if (breakpoint <= freq[j]){
  //     mu[j] = intercept;// + slope1*freq[j];
  //   }
  //   if (breakpoint > freq[j]){
  //     //mu[j] = intercept + breakpoint * (slope1 - slope2) + slope2 * freq[j];
  //     mu[j] = intercept + slope * (freq[j] - breakpoint);
  //   }
  // }
  
  mu = rep_vector(intercept, 21) + slope * log(freq);
  
  Sigma = diag_matrix(rep_vector(variances, 21));
  
  
  // Covariance
  // real<lower=0> st_dev;
  // cov_matrix[21] cov_mat;
  // corr_matrix[21] corr_mat;
  // cholesky_factor_cov[21] cov_mat_chol;
  // 
  // st_dev = sqrt(variances);
  // 
  // cov_mat_chol = st_dev*corr_mat_cholesky;
  // 
  // cov_mat = cov_mat_chol*cov_mat_chol';
  // corr_mat = corr_mat_cholesky * corr_mat_cholesky';
  
}


model {
  // c_concentration ~ exponential(1);
  // corr_mat_cholesky ~ lkj_corr_cholesky(c_concentration);
  

  for (i in 1:N) {
    y[i] ~ multi_normal(mu, Sigma);
  }
  
  
}


generated quantities {
  vector[N] log_lik;
  matrix[N,21] y_rep;         // Simulierte Daten
  
  for (i in 1:N) {
    log_lik[i] = multi_normal_lpdf(y[i] | mu, Sigma);
    y_rep[i] = to_row_vector(multi_normal_rng(mu, Sigma));
  }
  
}

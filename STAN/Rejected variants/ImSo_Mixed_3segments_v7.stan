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
}
parameters {
  // Fix
  real fe_intercept;
  real fe_slope1;
  real fe_d_slope2;
  real fe_d_slope3;
  real<lower=1,upper=5000> fe_bp1;
  real<lower=0> fe_bp_diff;

  // RE Hyperpriors
  real<lower=0> hyper_v_int;
  real<lower=0> hyper_v_slope1;
  real<lower=0> hyper_v_slope2;
  real<lower=0> hyper_v_slope3;
  real<lower=0> hyper_v_bp1;
  real<lower=0> hyper_v_bp_diff;
  // Kovarianz
  real<lower=0> c_concentration;
  vector<lower=0>[21] variances;
  cholesky_factor_corr[21] corr_mat_cholesky;
  // RE
  cholesky_factor_corr[6] REcorr_chol;
  vector<lower=0>[6] re_var;
  real<lower=0> re_c_concentration;
  matrix[N,6] re;
}
transformed parameters {
  real fe_slope2;
  real fe_slope3;

  vector[N] intercept;
  vector[N] slope1;

  vector[N] d_slope2;
  vector[N] d_slope3;

  vector[N] slope2;
  vector[N] slope3;
  vector[N] bp1;
  vector[N] bp_diff;
  vector[N] bp2;
  matrix[N,21] mu;

  
  // Data (CO)variance
  vector<lower=0>[21] st_dev;
  cov_matrix[21] cov_mat;
  corr_matrix[21] corr_mat;
  cholesky_factor_cov[21] cov_mat_chol;

  st_dev = sqrt(variances);

  cov_mat_chol = diag_matrix(st_dev)*corr_mat_cholesky;

  cov_mat = cov_mat_chol*cov_mat_chol';
  corr_mat = corr_mat_cholesky * corr_mat_cholesky';
  
  // RE cov
  cholesky_factor_cov[6] REcov_mat_chol;

  REcov_mat_chol = diag_matrix(sqrt(re_var))*REcorr_chol;
  
  cov_matrix[6] REcov;
  corr_matrix[6] REcorr;
  
  REcov = REcov_mat_chol*REcov_mat_chol';
  REcorr = REcorr_chol*REcorr_chol';
  

  // slopes
  fe_slope2 = fe_slope1 + fe_d_slope2;
  fe_slope3 = fe_slope2 + fe_d_slope3;

  // fix + re
  intercept = fe_intercept + re[,1];
  slope1 = fe_slope1 + re[,2];

  d_slope2 = fe_d_slope2 + re[,3];
  d_slope3 = fe_d_slope3 + re[,4];

  slope2 = slope1 + d_slope2;
  slope3 = slope2 + d_slope3;

  //slope2 = fe_slope2 + re_slope2;
  //slope3 = fe_slope3 + re_slope3;
  bp1 = fe_bp1 + re[,5];
  bp_diff = fe_bp_diff + re[,6];
  bp2 = bp1 + bp_diff;
  for (i in 1:N) {
  //real off12 = (slope1[i] - slope2[i]) * log(bp1[i]);
  //real off23 = (slope2[i] - slope3[i]) * log(bp2[i]);

  real off12 = -log(bp1[i]);
  real off23 = -log(bp2[i]);

  for (j in 1:21) {
    mu[i,j] = intercept[i] + slope1[i] * log_freq[j];

    if (freq[j] > bp1[i]) {
      mu[i,j] = mu[i,j] + d_slope2[i] * (log_freq[j]+off12);
      }
    if (freq[j] > bp2[i]) {
      mu[i,j] = mu[i,j] + d_slope3[i] * (log_freq[j]+off23);
      }
  }


}
}
model {
  fe_intercept ~ normal(92, 4);
  fe_slope1 ~ normal(-10,1);
  fe_d_slope2 ~ normal(5,1);
  fe_d_slope3 ~ normal(-14, 0.6);
  fe_bp1 ~ normal(100, 11);
  fe_bp_diff ~ normal(1300, 61);


  hyper_v_int ~ exponential(0.01);
  hyper_v_slope1 ~ exponential(0.01);
  hyper_v_slope2 ~ exponential(0.01);
  hyper_v_slope3 ~ exponential(0.01);
  hyper_v_bp1 ~ exponential(0.01);
  hyper_v_bp_diff ~ exponential(0.01);
  re_var[1] ~ exponential(1.0/hyper_v_int);
  re_var[2] ~ exponential(1.0/hyper_v_slope1);
  re_var[3] ~ exponential(1.0/hyper_v_slope2);
  re_var[4] ~ exponential(1.0/hyper_v_slope3);
  re_var[5] ~ exponential(1.0/hyper_v_bp1);
  re_var[6] ~ exponential(1.0/hyper_v_bp_diff);

  // RE
  re_c_concentration ~ exponential(1);
  REcorr_chol ~  lkj_corr_cholesky(re_c_concentration);

  // Data covariance
  for (j in 1:21) {
    variances[j] ~ exponential(0.1);
  }
  c_concentration ~ exponential(1);
  corr_mat_cholesky ~ lkj_corr_cholesky(c_concentration);
  // Data and Re

  for (i in 1:N) {

    re[i] ~ multi_normal_cholesky(rep_vector(0, 6), REcov_mat_chol);

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


#include "ImpactModel.stan"

data {
  int<lower=0> N;
  matrix[N,21] y;
  vector[N] l1;
  vector[N] l2;
  real rho;
}


transformed data {
  vector[21] freq;
  real c0;
  real cL;
  real rho0;
  real max_rl;
  real reference_f;
  real offset;
  
  // fix
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
  c0 = 340;
  cL = 3500;
  rho0 = 1.23;
  max_rl = 2;
  reference_f = 1000;
  offset = 155;
  
  //calcu
  real fc;
  vector[N] f11;
  matrix[N,21] rad_fac;
  matrix[N,21] rad_fac3;
  
  fc = fc_fun(c0, cL);
  
  
  for (i in 1:N) {
    f11[i] = f11_fun(c0, fc, l1[i], l2[i]);
    
    for (j in 1:21) {
      rad_fac3[i,j] = rad3(freq[j], fc, l1[i], l2[i], c0);
      rad_fac[i,j] = radiation_fac(freq[j], fc, f11[i], c0,
                                    l1[i], l2[i], max_rl, rad_fac3[i,j]);
    }
  }
  

}


parameters {
  // Physical
  real<lower=0,upper=1> thickness;
  real<lower=0> q;
  real<lower=50> reso_freq;
  
  // Covariance
  real<lower=0> c_concentration;
  vector<lower=0>[21] variances;
  cholesky_factor_corr[21] corr_mat_cholesky;
}




transformed parameters {
  // Physical
  real mass;
  matrix[N,21] Ln;
  matrix[N,21] L_delta;
  matrix[N,21] Ln_dash;
  
  mass = thickness * rho;
  
  for (i in 1:N){
    for (j in 1:21) {
      Ln[i,j] = impact_sound_L_n(freq[j], mass, reference_f, offset, fc, rad_fac[i,j]);
      L_delta[i,j] = impact_sound_L_delta(freq[j], q, reso_freq);
      
      Ln_dash[i,j] = Ln[i,j] - L_delta[i,j];
    }
  }
  
  // Covariance
  vector<lower=0>[21] st_dev;
  cov_matrix[21] cov_mat;
  corr_matrix[21] corr_mat;
  cholesky_factor_cov[21] cov_mat_chol;
  
  st_dev = sqrt(variances);
  
  cov_mat_chol = diag_matrix(st_dev)*corr_mat_cholesky;
  
  cov_mat = cov_mat_chol*cov_mat_chol';
  corr_mat = corr_mat_cholesky * corr_mat_cholesky';
  
}


model {
  c_concentration ~ exponential(1);
  corr_mat_cholesky ~ lkj_corr_cholesky(c_concentration);
  

  for (i in 1:N) {
    y[i] ~ multi_normal_cholesky(Ln_dash[i], cov_mat_chol);
  }
}

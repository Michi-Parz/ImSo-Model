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

// Auskommentieren wenn ohne Trafo
parameters {
  // // Physical: Fixed effects
  // real<lower=0,upper=1> fix_thickness;
  // real<lower=0> fix_q;
  // real<lower=50> fix_reso_freq;
  // // Physical: Random effects
  // vector<lower = -fix_thickness, upper = 1-fix_thickness>[N] rando_thickness;
  // vector<lower = -fix_q>[N] rando_q;
  // vector<lower=50-fix_reso_freq>[N] rando_reso_freq;
  
  // Physical: Fixed effects
  real<lower =0> fix_thickness;
  real fix_q;
  real<lower = 0> fix_reso_freq;
  // Hyperpara
  real<lower=0> hyp_t;
  real<lower=0> hyp_q;
  real<lower=0> hyp_f0;
  // Physical: Random effects
  vector<lower=-fix_thickness>[N] rando_thickness;
  vector[N] rando_q;
  vector<upper = 5000-fix_reso_freq>[N] rando_reso_freq;

  // Covariance
  //real<lower=0> c_concentration;
  real<lower=0> variances;
  //cholesky_factor_corr[21] corr_mat_cholesky;

}


// // Auskommentieren wenn mit Trafo
// parameters {
//   // Physical: Fixed effects
//   real fix_thick_trafo;
//   real fix_q_trafo;
//   real fix_rf_trafo;
//   // Physical: Random effects
//   vector[N] rando_thick_trafo;
//   vector[N] rando_q_trafo;
//   vector[N] rando_rf_trafo;
//   
//   // Covariance
//   real<lower=0> c_concentration;
//   vector<lower=0>[21] variances;
//   cholesky_factor_corr[21] corr_mat_cholesky;
// }



transformed parameters {
  
  // // AUskommentieren wenn ohne Trafo
  // // Fixed + Random
  // vector[N] thickness_trafo;
  // vector[N] q_trafo;
  // vector[N] rf_trafo;
  
  // Fixed + Random
  vector[N] thickness; // [0,1]
  vector[N] q; //>0
  vector[N] reso_freq; //[50,5000]
  
  // Physical
  vector[N] mass;
  matrix[N,21] Ln;
  matrix[N,21] L_delta;
  matrix[N,21] Ln_dash;

  
  for (i in 1:N){
    // // Auskommentieren wenn ohne Trafo
    // thickness_trafo[i] = fix_thick_trafo + rando_thick_trafo[i];
    // q_trafo[i] = fix_q_trafo + rando_q_trafo[i];
    // rf_trafo[i] = fix_rf_trafo + rando_rf_trafo[i];
    // thickness[i] = 1/(1+exp(-thickness_trafo[i]));
    // q[i] = exp(q_trafo[i]);
    // reso_freq[i] = 50 + exp(rf_trafo[i]);
    
    
    // Auskommentieren wenn Trafo
    thickness[i] = fix_thickness + rando_thickness[i];
    q[i] = fix_q + rando_q[i];
    reso_freq[i] = fix_reso_freq + rando_reso_freq[i];
    
    
    mass[i] = thickness[i] * rho;
    for (j in 1:21) {
      Ln[i,j] = impact_sound_L_n(freq[j], mass[i], reference_f, offset, fc, rad_fac[i,j]);
      L_delta[i,j] = impact_sound_L_delta(freq[j], q[i], reso_freq[i]);
      
      Ln_dash[i,j] = Ln[i,j] - L_delta[i,j];
    }
  }
  
  // Covariance
  // vector<lower=0>[21] st_dev;
  // cov_matrix[21] cov_mat;
  // corr_matrix[21] corr_mat;
  // cholesky_factor_cov[21] cov_mat_chol;
  // 
  // st_dev = sqrt(variances);
  // 
  // cov_mat_chol = diag_matrix(st_dev)*corr_mat_cholesky;
  // 
  // cov_mat = cov_mat_chol*cov_mat_chol';
  // corr_mat = corr_mat_cholesky * corr_mat_cholesky';
  
}


// Versuche die nichts gebracht haben:
// random effekt als Normalverteilt mit erwartungswert 0 angenommen.
// Ohne Fix effekt sondern nur die Parameter als Vektoren
// Transformierte Parameter damit eigentlicher Fit unbegrenzt ist
// Grenzen angepasst
// IID Fehlerannahme (Eine Varianz und keine Korrelationen)

model {
  // c_concentration ~ exponential(1);
  // corr_mat_cholesky ~ lkj_corr_cholesky(c_concentration);
  
  rando_thickness ~ multi_normal(rep_vector(0, N), diag_matrix(rep_vector(hyp_t, N)));
  rando_q ~ multi_normal(rep_vector(0, N), diag_matrix(rep_vector(hyp_q, N)));
  rando_reso_freq ~ multi_normal(rep_vector(0, N), diag_matrix(rep_vector(hyp_f0, N)));

  for (i in 1:N) {
    //y[i] ~ multi_normal_cholesky(Ln_dash[i], cov_mat_chol);
    y[i] ~ multi_normal(Ln_dash[i], diag_matrix(rep_vector(variances, 21)));
  }
  
  
}

data {
  int<lower=0> N;
  matrix[N,21] y;
}

transformed data {
  vector[21] freq;
  freq = [50,63,80,100,125, 160, 200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000]';
}


parameters {
  // Fix
  real fe_intercept;
  real fe_slope;
  
  // Random
  vector[N] re_intercept;
  vector[N] re_slope;
  
  // Hyper
  real<lower=0> hyp_i;
  real<lower=0> hyp_s;

  // Varianz
  real<lower=0> variances;
  
}


transformed parameters {
  vector[N] intercept;
  vector[N] slope;
  matrix[N,21] mu;
  matrix[21,21] Sigma;
  
  intercept = rep_vector(fe_intercept, N) + re_intercept;
  slope = rep_vector(fe_slope, N) + re_slope;
  
  for (i in 1:N) {
    mu[i] = to_row_vector(rep_vector(intercept[i], 21) + slope[i] * log(freq));
  }
  
  
  Sigma = diag_matrix(rep_vector(variances, 21));

}


model {
  
  re_intercept ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(hyp_i, N)));
  re_slope ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(hyp_s, N)));

  for (i in 1:N) {
    y[i] ~ multi_normal(mu[i], Sigma);
  }
  
  
}


generated quantities {
  vector[N] log_lik;
  //matrix[N,21] y_rep;         // Simulierte Daten
  
  for (i in 1:N) {
    log_lik[i] = multi_normal_lpdf(y[i] | mu[i], Sigma);
    //y_rep[i] = to_row_vector(multi_normal_rng(mu, Sigma));
  }
  
}

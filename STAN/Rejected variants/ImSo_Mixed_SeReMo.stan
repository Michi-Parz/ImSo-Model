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
  real fe_slope1;
  real fe_slope2;
  real<lower=1,upper=5000> fe_bp;
  
  // Hyper
  real<lower=0> var_int;
  real<lower=0> var_slope1;
  real<lower=0> var_slope2;
  real<lower=0> var_bp;
  
  // Random
  vector[N] re_intercept;
  vector[N] re_slope1;
  vector[N] re_slope2;
  vector<lower=1-fe_bp,upper=5000-fe_bp>[N] re_bp;

  // Varianz
  real<lower=0> variances;
  
}


transformed parameters {
  vector[N] intercept;
  vector[N] slope1;
  vector[N] slope2;
  vector[N] bp;
  matrix[N,21] mu;
  matrix[21,21] Sigma;
  
  intercept = rep_vector(fe_intercept, N) + re_intercept;
  slope1 = rep_vector(fe_slope1, N) + re_slope1;
  slope2 = rep_vector(fe_slope2, N) + re_slope2;
  bp = rep_vector(fe_bp, N) + re_bp;
  
  for (i in 1:N) {
    for (j in 1:21) {
      if (freq[j] <= bp[i]) {
        mu[i,j] = intercept[i] + slope1[i] * log(freq[j]);
      }
      if (freq[j] > bp[i]) {
        mu[i,j] = intercept[i] + slope2[i] * log(freq[j]) + (slope1[i]-slope2[i])*log(bp[i]);
      }
      
    }
  }
  
  
  Sigma = diag_matrix(rep_vector(variances, 21));

}


model {
  
  re_intercept ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_int, N)));
  re_slope1 ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_slope1, N)));
  re_slope2 ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_slope2, N)));
  re_bp ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_bp, N)));
                  

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

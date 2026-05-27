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
  real fe_slope2;
  real fe_slope3;
  real<lower=1,upper=5000> fe_bp1;
  real<lower=0> fe_bp_diff;
  
  // Hyper
  real<lower=0> var_int;
  real<lower=0> var_slope1;
  real<lower=0> var_slope2;
  real<lower=0> var_slope3;
  real<lower=0> var_bp1;
  real<lower=0> var_bp_diff;
  
  // Random
  vector[N] re_intercept;
  vector[N] re_slope1;
  vector[N] re_slope2;
  vector[N] re_slope3;
  vector<lower = 1-fe_bp1,upper=5000-fe_bp1>[N] re_bp1;
  vector<lower = -var_bp_diff>[N] re_bp_diff;

  // Varianz
  real<lower=0> variances;
  
}


transformed parameters {
  vector[N] intercept;
  vector[N] slope1;
  vector[N] slope2;
  vector[N] slope3;
  vector[N] bp1;
  vector[N] bp_diff;
  vector[N] bp2;
  matrix[N,21] mu;
  matrix[21,21] Sigma;
  
  intercept = fe_intercept + re_intercept; // rep_vector(...,N) lt ChatGPT unnötig deshalb weg
  slope1 = fe_slope1 + re_slope1;
  slope2 = fe_slope2 + re_slope2;
  slope3 = fe_slope3 + re_slope3;
  bp1 = fe_bp1 + re_bp1;
  bp_diff = fe_bp_diff + re_bp_diff;
  bp2 = bp1 + bp_diff;
  
  
  
  
  for (i in 1:N) {
  real off12 = (slope1[i] - slope2[i]) * log(bp1[i]);
  real off23 = (slope2[i] - slope3[i]) * log(bp2[i]);

  for (j in 1:21) {
    if (freq[j] <= bp1[i]) {
      mu[i, j] = intercept[i] + slope1[i] * log_freq[j];
    } else if (freq[j] <= bp2[i]) {
      mu[i, j] = intercept[i] + slope2[i] * log_freq[j] + off12;
    } else {
      mu[i, j] = intercept[i] + slope3[i] * log_freq[j] + off12 + off23;
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
  re_slope3 ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_slope3, N)));
  re_bp1 ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_bp1, N)));
  re_bp_diff ~ multi_normal(rep_vector(0, N),
                  diag_matrix(rep_vector(var_bp_diff, N)));
                  

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

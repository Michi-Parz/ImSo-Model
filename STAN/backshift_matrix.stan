
matrix backshift_matrix (int N, vector coefficients){
  int M = size(coefficients);
  
  matrix[N,N] back_mat;
  
  back_mat = rep_matrix(0, N,N);
  
  int i=1;
  int j=1;
  
  for (m in 1:M) {
    back_mat[i+j,i] = coefficients[m];
    i += 1;
    if (i+j>N) {
      i = 1;
      j += 1;
    }
  }
  
  return back_mat;
}


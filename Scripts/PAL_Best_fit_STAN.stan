// generated with brms 2.19.0
functions {
  /* zero-inflated negative binomial log-PDF of a single response
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_lpmf(int y, real mu, real phi,
                                       real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   mu: mean parameter of negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_logit_lpmf(int y, real mu,
                                             real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_lpmf(0 | mu, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_lpmf(y | mu, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: zero-inflation probability
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_lpmf(int y, real eta,
                                           real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_lpmf(1 | zi),
                         bernoulli_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  /* zero-inflated negative binomial log-PDF of a single response
   * log parameterization for the negative binomial part
   * logit parameterization of the zero-inflation part
   * Args:
   *   y: the response value
   *   eta: linear predictor for negative binomial distribution
   *   phi: shape parameter of negative binomial distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real zero_inflated_neg_binomial_log_logit_lpmf(int y, real eta,
                                                 real phi, real zi) {
    if (y == 0) {
      return log_sum_exp(bernoulli_logit_lpmf(1 | zi),
                         bernoulli_logit_lpmf(0 | zi) +
                         neg_binomial_2_log_lpmf(0 | eta, phi));
    } else {
      return bernoulli_logit_lpmf(0 | zi) +
             neg_binomial_2_log_lpmf(y | eta, phi);
    }
  }
  // zero_inflated negative binomial log-CCDF and log-CDF functions
  real zero_inflated_neg_binomial_lccdf(int y, real mu, real phi, real hu) {
    return bernoulli_lpmf(0 | hu) + neg_binomial_2_lccdf(y | mu, phi);
  }
  real zero_inflated_neg_binomial_lcdf(int y, real mu, real phi, real hu) {
    return log1m_exp(zero_inflated_neg_binomial_lccdf(y | mu, phi, hu));
  }
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K_B0;  // number of population-level effects
  matrix[N, K_B0] X_B0;  // population-level design matrix
  int<lower=1> K_Tmin;  // number of population-level effects
  matrix[N, K_Tmin] X_Tmin;  // population-level design matrix
  int<lower=1> K_Tmax;  // number of population-level effects
  matrix[N, K_Tmax] X_Tmax;  // population-level design matrix
  int<lower=1> K_alphaii0;  // number of population-level effects
  matrix[N, K_alphaii0] X_alphaii0;  // population-level design matrix
  int<lower=1> K_alphaij;  // number of population-level effects
  matrix[N, K_alphaij] X_alphaij;  // population-level design matrix
  int<lower=1> K_alphaiiT;  // number of population-level effects
  matrix[N, K_alphaiiT] X_alphaiiT;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  int<lower=1> K_zi;  // number of population-level effects
  matrix[N, K_zi] X_zi;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_zi = K_zi - 1;
  matrix[N, Kc_zi] Xc_zi;  // centered version of X_zi without an intercept
  vector[Kc_zi] means_X_zi;  // column means of X_zi before centering
  for (i in 2:K_zi) {
    means_X_zi[i - 1] = mean(X_zi[, i]);
    Xc_zi[, i - 1] = X_zi[, i] - means_X_zi[i - 1];
  }
}
parameters {
  vector<lower=1e-04,upper=200>[K_B0] b_B0;  // population-level effects
  vector<lower=10,upper=22.5>[K_Tmin] b_Tmin;  // population-level effects
  vector<lower=27.9,upper=35>[K_Tmax] b_Tmax;  // population-level effects
  vector<lower=1e-05,upper=20>[K_alphaii0] b_alphaii0;  // population-level effects
  vector<lower=1e-05,upper=20>[K_alphaij] b_alphaij;  // population-level effects
  vector<lower=-2,upper=2>[K_alphaiiT] b_alphaiiT;  // population-level effects
  real<lower=0> shape;  // shape parameter
  vector[Kc_zi] b_zi;  // population-level effects
  real Intercept_zi;  // temporary intercept for centered predictors
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b_B0 | 80, 15)
    - 1 * log_diff_exp(normal_lcdf(200 | 80, 15), normal_lcdf(1e-04 | 80, 15));
  lprior += normal_lpdf(b_Tmin | 20, 50)
    - 1 * log_diff_exp(normal_lcdf(22.5 | 20, 50), normal_lcdf(10 | 20, 50));
  lprior += normal_lpdf(b_Tmax | 30, 50)
    - 1 * log_diff_exp(normal_lcdf(35 | 30, 50), normal_lcdf(27.9 | 30, 50));
  lprior += normal_lpdf(b_alphaii0 | 0.00001, 0.5)
    - 1 * log_diff_exp(normal_lcdf(20 | 0.00001, 0.5), normal_lcdf(1e-05 | 0.00001, 0.5));
  lprior += normal_lpdf(b_alphaij | 0.00001, 1)
    - 1 * log_diff_exp(normal_lcdf(20 | 0.00001, 1), normal_lcdf(1e-05 | 0.00001, 1));
  lprior += normal_lpdf(b_alphaiiT | 0, 0.05)
    - 1 * log_diff_exp(normal_lcdf(2 | 0, 0.05), normal_lcdf(-2 | 0, 0.05));
  lprior += gamma_lpdf(shape | 0.01, 0.01);
  lprior += logistic_lpdf(Intercept_zi | 0, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_B0 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_Tmin = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_Tmax = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_alphaii0 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_alphaij = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_alphaiiT = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    // initialize linear predictor term
    vector[N] zi = rep_vector(0.0, N);
    nlp_B0 += X_B0 * b_B0;
    nlp_Tmin += X_Tmin * b_Tmin;
    nlp_Tmax += X_Tmax * b_Tmax;
    nlp_alphaii0 += X_alphaii0 * b_alphaii0;
    nlp_alphaij += X_alphaij * b_alphaij;
    nlp_alphaiiT += X_alphaiiT * b_alphaiiT;
    zi += Intercept_zi + Xc_zi * b_zi;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (C_1[n] * (nlp_B0[n] * C_2[n] * (C_2[n] - nlp_Tmin[n]) * sqrt(nlp_Tmax[n] - C_2[n])) * (1 / ((1 + (nlp_alphaii0[n] + nlp_alphaiiT[n] * C_3[n]) * C_4[n] + nlp_alphaij[n] * C_5[n]))));
    }
    for (n in 1:N) {
      target += zero_inflated_neg_binomial_logit_lpmf(Y[n] | mu[n], shape, zi[n]);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_zi_Intercept = Intercept_zi - dot_product(means_X_zi, b_zi);
}


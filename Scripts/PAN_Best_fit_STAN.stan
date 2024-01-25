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
  int<lower=1> K_alphaii;  // number of population-level effects
  matrix[N, K_alphaii] X_alphaii;  // population-level design matrix
  int<lower=1> K_alphaij0;  // number of population-level effects
  matrix[N, K_alphaij0] X_alphaij0;  // population-level design matrix
  int<lower=1> K_alphaijT;  // number of population-level effects
  matrix[N, K_alphaijT] X_alphaijT;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  int<lower=1> K_shape;  // number of population-level effects
  matrix[N, K_shape] X_shape;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_shape = K_shape - 1;
  matrix[N, Kc_shape] Xc_shape;  // centered version of X_shape without an intercept
  vector[Kc_shape] means_X_shape;  // column means of X_shape before centering
  for (i in 2:K_shape) {
    means_X_shape[i - 1] = mean(X_shape[, i]);
    Xc_shape[, i - 1] = X_shape[, i] - means_X_shape[i - 1];
  }
}
parameters {
  vector<lower=1e-04,upper=200>[K_B0] b_B0;  // population-level effects
  vector<lower=10,upper=22.5>[K_Tmin] b_Tmin;  // population-level effects
  vector<lower=28.7,upper=40>[K_Tmax] b_Tmax;  // population-level effects
  vector<lower=1e-05,upper=20>[K_alphaii] b_alphaii;  // population-level effects
  vector<lower=1e-05,upper=20>[K_alphaij0] b_alphaij0;  // population-level effects
  vector<lower=-2,upper=2>[K_alphaijT] b_alphaijT;  // population-level effects
  vector[Kc_shape] b_shape;  // population-level effects
  real Intercept_shape;  // temporary intercept for centered predictors
  real<lower=0,upper=1> zi;  // zero-inflation probability
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(b_B0 | 80, 15)
    - 1 * log_diff_exp(normal_lcdf(200 | 80, 15), normal_lcdf(1e-04 | 80, 15));
  lprior += normal_lpdf(b_Tmin | 20, 50)
    - 1 * log_diff_exp(normal_lcdf(22.5 | 20, 50), normal_lcdf(10 | 20, 50));
  lprior += normal_lpdf(b_Tmax | 30, 50)
    - 1 * log_diff_exp(normal_lcdf(40 | 30, 50), normal_lcdf(28.7 | 30, 50));
  lprior += normal_lpdf(b_alphaii | 0.00001, 1)
    - 1 * log_diff_exp(normal_lcdf(20 | 0.00001, 1), normal_lcdf(1e-05 | 0.00001, 1));
  lprior += normal_lpdf(b_alphaij0 | 0.00001, 1)
    - 1 * log_diff_exp(normal_lcdf(20 | 0.00001, 1), normal_lcdf(1e-05 | 0.00001, 1));
  lprior += normal_lpdf(b_alphaijT | 0, 0.05)
    - 1 * log_diff_exp(normal_lcdf(2 | 0, 0.05), normal_lcdf(-2 | 0, 0.05));
  lprior += student_t_lpdf(Intercept_shape | 3, 0, 2.5);
  lprior += beta_lpdf(zi | 1, 1);
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
    vector[N] nlp_alphaii = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_alphaij0 = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_alphaijT = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    // initialize linear predictor term
    vector[N] shape = rep_vector(0.0, N);
    nlp_B0 += X_B0 * b_B0;
    nlp_Tmin += X_Tmin * b_Tmin;
    nlp_Tmax += X_Tmax * b_Tmax;
    nlp_alphaii += X_alphaii * b_alphaii;
    nlp_alphaij0 += X_alphaij0 * b_alphaij0;
    nlp_alphaijT += X_alphaijT * b_alphaijT;
    shape += Intercept_shape + Xc_shape * b_shape;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (C_1[n] * (nlp_B0[n] * C_2[n] * (C_2[n] - nlp_Tmin[n]) * sqrt(nlp_Tmax[n] - C_2[n])) * (1 / ((1 + nlp_alphaii[n] * C_3[n] + (nlp_alphaij0[n] + nlp_alphaijT[n] * C_4[n]) * C_5[n]))));
    }
    shape = exp(shape);
    for (n in 1:N) {
      target += zero_inflated_neg_binomial_lpmf(Y[n] | mu[n], shape[n], zi);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_shape_Intercept = Intercept_shape - dot_product(means_X_shape, b_shape);
}


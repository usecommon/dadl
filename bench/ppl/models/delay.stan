// Gaussian peak on a flat background, Poisson counts per bin (times in ns).
data {
  int<lower=1> K;
  vector[K] centres;
  real<lower=0> width;
  array[K] int<lower=0> counts;
  vector[4] centre;
}
parameters {
  real log_n;
  real mu;
  real log_sigma;
  real log_b;
}
transformed parameters {
  real n = exp(log_n);
  real sigma = exp(log_sigma);
  real b = exp(log_b);
}
model {
  log_n ~ normal(centre[1], 10);
  mu ~ normal(centre[2], 10);
  log_sigma ~ normal(centre[3], 10);
  log_b ~ normal(centre[4], 10);
  counts ~ poisson(b + n * width * exp(-0.5 * square((centres - mu) / sigma)) / (sigma * sqrt(2 * pi())));
}

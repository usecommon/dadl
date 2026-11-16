// DADL dark inter-arrival likelihood (dadl.models.DarkAfterpulseModel), times in microseconds.
data {
  int<lower=1> N;
  vector<lower=0>[N] s;
  vector[3] centre;
}
transformed data {
  real s_sum = sum(s);
}
parameters {
  real log_r;
  real logit_p;
  real log_tau;
}
transformed parameters {
  real r_per_us = exp(log_r);
  real p = inv_logit(logit_p);
  real tau_us = exp(log_tau);
}
model {
  log_r ~ normal(centre[1], 10);
  logit_p ~ normal(centre[2], 10);
  log_tau ~ normal(centre[3], 10);
  target += sum(log(r_per_us * (1 - p) + p * (r_per_us + 1 / tau_us) * exp(-s / tau_us))) - r_per_us * s_sum;
}

data {
  int<lower=1> trials;
  int<lower=1> frames;
  matrix[trials,frames] signals;
  int<lower=0,upper=1> choices[trials];
}

parameters {
  real<lower=0,upper=1> prior_C;
  real<lower=-1,upper=1> gam;
  real<lower=0> bound;
  real<lower=0> noise;
  real<lower=0> scale;
  real<lower=0,upper=1> lapse[2];
  real<lower=0> temperature;

  // rather than integrate out the noise, we'll infer it
  matrix[trials,frames] noise_eps;
}

model {
  // Priors on parameters
  prior_C ~ beta(2, 2);
  gam ~ uniform(-1, 1);
  bound ~ gamma(2, 0.3333);
  noise ~ gamma(1, 4);
  scale ~ gamma(1, 0.05);
  lapse ~ beta(1, 10);
  temperature ~ gamma(1, 0.24);

  // Reparameterization of per-trial noise
  to_vector(noise_eps) ~ normal(0, 1);

  for(t in 1:trials){
    real log_post = log(prior_C) - log(1-prior_C);
    for(f in 1:frames){
      if(log_post >= bound) {
        log_post = bound;
      } else if(log_post <= -bound) {
        log_post = -bound;
      } else{
        // integration block
        log_post = log_post * (1-gam) + noise_eps[t,f]*noise + signals[t,f] / scale;
      }
    }
    real p_choice = 1.0/(1.0+exp(-log_post / temperature));
    p_choice = lapse[1] + (1-lapse[1]-lapse[2])*p_choice;
    choices[t] ~ bernoulli(p_choice);
  }
}
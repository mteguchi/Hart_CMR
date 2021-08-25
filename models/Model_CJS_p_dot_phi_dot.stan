// from https://github.com/stan-dev/example-models/tree/master/BPA
// This is p(dot) phi(dot) model.  I removed some deleted lines to
// make the code more legible.
//
// Added dt vector as a data input for unequal time between sampling

// 2019-09-23


functions {
// a function to find the first capture period of each individual
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }

// a function to find out the last capture period of each individual
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }

// a function to calculate the probability of non-captures
  matrix prob_uncaptured(int nind, int n_occasions,
                         matrix p, matrix phi) {
    matrix[nind, n_occasions] chi;

    for (i in 1:nind) {
      chi[i, n_occasions] = 1.0;   // last occasion is set to 1 (present)

      // this loop goes back in time
      for (t in 1:(n_occasions - 1)) {
        int t_curr = n_occasions - t;   // "current" occasion
        int t_next = t_curr + 1;        // "next" occasion

        // p[i, t_next-1]  can be p[i, t_curr]?
        chi[i, t_curr] = (1 - (phi[i, t_curr]))
                        + (phi[i, t_curr]) * (1 - p[i, t_next - 1]) * chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int<lower=0> nind;            // Number of individuals
  int<lower=2> n_occasions;     // Number of capture occasions
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
  real<lower = 0> dt[n_occasions-1];  // time between captures in years
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];
  int<lower=0,upper=n_occasions> last[nind];

  for (i in 1:nind)
    first[i] = first_capture(y[i]);
  for (i in 1:nind)
    last[i] = last_capture(y[i]);
}

parameters {
  real<lower=0,upper=1> mean_phi;    // Mean survival
  real<lower=0,upper=1> mean_p;      // Mean recapture
}

transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi;
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] p;
  matrix<lower=0,upper=1>[nind, n_occasions] chi;

  // Constraints
  for (i in 1:nind) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:n_occ_minus_1) {
      phi[i, t] = mean_phi ^ dt[t];
      p[i, t] = mean_p;
    }
  }

  chi = prob_uncaptured(nind, n_occasions, p, phi);
}

model {
  // Priors
  // Uniform priors are implicitly defined.
  //  mean_phi ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);

  // Likelihood
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (t in (first[i] + 1):last[i]) {
        1 ~ bernoulli(phi[i, t - 1]);
        y[i, t] ~ bernoulli(p[i, t - 1]);
      }
      1 ~ bernoulli(chi[i, last[i]]);
    }
  }
}

library(MASS)

source('R/matlab_connect.R')

RAND_MODE_NATIVE <- 1 # native
RAND_MODE_CACHED <- 2 # from saved file
RAND_MODE_MATLAB <- 3 # from MATLAB

RAND_MODE <- RAND_MODE_MATLAB

MATLAB_CONNECT <- FALSE

if (MATLAB_CONNECT) {
  matlab <- MaltlabConnection()
  matlab$start_server()
  matlab$connect()
}

## MATLAB compatibility functions

load_rand <- function(dist = "normal", rows = 1, cols = 1) {
  filename <- sprintf("./data/rand_%s_%d_%d.csv", dist, rows, cols)
  print(paste("loading ", filename))
  as.matrix(read.csv(filename, header = FALSE))
}

# matlab compatiliblity functions
ones <- function(nrow = 1, ncol = 1) { matrix(1, nrow, ncol) }
zeros <- function(nrow = 1, ncol = 1) { matrix(0, nrow, ncol) }
eye <- function(d) { diag(1, d) }
cdiv <- function(x, y) { t(t(x)/y) }

rand <- function(rows, cols) { 
  switch(RAND_MODE,
         RAND_MODE_NATIVE = matrix(runif(rows * cols), rows, cols),
         RAND_MODE_CACHED = load_rand("uniform", rows, cols),
         RAND_MODE_MATLAB = matlab$rand(rows, cols)
  )
}

randn <- function(rows, cols) { 
  switch(RAND_MODE,
         RAND_MODE_NATIVE = matrix(rnorm(rows * cols), rows, cols),
         RAND_MODE_CACHED = load_rand("normal", rows, cols),
         RAND_MODE_MATLAB = matlab$randn(rows, cols)
  )
}

pinv <- function(x) {
  #ginv(x)
  #pracma::pinv(x)
  corpcor::pseudoinverse(x)
}

checkpoint_id <- 0

matrix_checkpoint <- function(x) {
  return()
  
  # default tolerance = .Machine$double.eps ^ 0.5
  DEFAULT_Tolerance <- .Machine$double.eps ^ 0.5
  tolerance <- .Machine$double.eps ^ 0.1
  
  elementwise.all.equal <- Vectorize(function(x, y) { isTRUE(all.equal(x, y, tolerance)) })

  filename <- file.path("/tmp", sprintf("checkpoint_%03d", checkpoint_id))
  x_ref <- as.matrix(read.csv(filename, header = FALSE))
  if (isTRUE(all.equal(matrix(x), matrix(x_ref), tolerance = tolerance))) {
    #message(sprintf("%s: OK", filename))
  } else {
    if (nrow(x) != nrow(x_ref)) {
      browser()
    }
    m <- elementwise.all.equal(x, x_ref)
    #message(sprintf("%s: msimatch (%d/%d items)", filename, length(m[m==FALSE]), length(m)))
    #browser()
  }
  checkpoint_id <<- checkpoint_id + 1
}

base <- function(k, m, d) {
  mm <- m ^ c((d - 1):0)
  v <- rep(NA, d)
  
  for (i in 1:d) {
    v[i] <- trunc(k / mm[i])
    k <- k - mm[i] * v[i]
  }
  return(v + 1)
}

rand_matrix <- function(x, rows, cols) {
  matrix(sample(x, rows * cols, replace = TRUE), rows, cols)
}

FactorialHmm_OLD <- function()
{
  env <- environment()
  
  alpha <- NA

  this <- list(
    env = env,
    
    # get/set value
    get_value = function(key) {
      return(get(key, envir = env))
    },
    set_value = function(key, value) {
      return(assign(key, value, envir = env))
    },
    
    estimate = function(observations, num_chains = 2, num_states = num_chains, T = NULL, max_iterations = 100, tolerance = 0.0001) {
      
      matlab$set.seed(666)
      
      p <- ncol(observations)
      N <- nrow(observations)
      
      if (is.null(T))
        T <- N
      
      N <- N / T
      
      matrix_checkpoint(observations)      
      
      covariance_matrix <- cov(observations)
      if (length(covariance_matrix) > 1) {
        covariance_matrix <- diag(diag(covariance_matrix))
      }
      
      XX <- t(observations) %*% observations / (N * T)
      matrix_checkpoint(XX)      
      
      rand_matrix <- randn(num_chains * num_states, p)
      mu <- ((rand_matrix %*% sqrt(covariance_matrix)) / num_chains) + (ones(num_states * num_chains, 1) %*% (colMeans(observations) / num_chains));
      
      matrix_checkpoint(mu)      
      
      priors <- rand(num_states, num_chains)
      priors <- cdiv(priors, colSums(priors))
      
      matrix_checkpoint(priors)      
      
      transition_matrix <- rand(num_states * num_chains, num_states)
      transition_matrix <- transition_matrix / rowSums(transition_matrix)
      
      matrix_checkpoint(transition_matrix)      
      
      log_likelihood_trace <- c()
      log_likelihood <- 0
      log_likelihood_base <- 0
      
      dd <- zeros(num_states ^ num_chains, num_chains)
      for (i in 1:(num_states ^ num_chains)) {
        dd[i,] <- base(i - 1, num_states, num_chains)
      }
      
      matrix_checkpoint(dd)
      
      alpha <- zeros(T,num_states ^ num_chains)
      B <- zeros(T,num_states ^ num_chains)
      beta <- zeros(T,num_states ^ num_chains)
      gamma <- zeros(T,num_states ^ num_chains)
      B2 <- zeros(T * num_states ^ num_chains,num_states ^ num_chains)
      
      eta <- zeros(T * num_states * num_chains,num_states * num_chains)
      
      Mub <- zeros(num_states ^ num_chains,p)
      Pb <- ones(num_states ^ num_chains,num_states ^ num_chains)
      Pib <- ones(num_states ^ num_chains,1)
      collapse <- zeros(num_states ^ num_chains, num_chains * num_states)
      collapse2 <- zeros(num_states ^ num_chains, num_chains * num_states * num_chains * num_states)
      
      for (i in 1:(num_states ^ num_chains)) {
        dd[i,] <- base(i - 1, num_states, num_chains)
        for (j in 1:num_chains) {
          Mub[i,] <- as.matrix(Mub[i,] + mu[(j - 1) * num_states + dd[i,j],])
          Pib[i,] <- as.matrix(Pib[i,] * priors[dd[i,j],j])
          collapse[i, (j - 1) * num_states + dd[i,j]] <- 1
          for (l in 1:num_chains) {
            collapse2[i,((j - 1) * num_states + dd[i,j] - 1) * num_chains * num_states + (l - 1) * num_states + dd[i,l]] <- 1
          }
        }
        for (j in 1:(num_states ^ num_chains)) {
          for (l in 1:num_chains) {
            Pb[i,j] <- Pb[i,j] %*% transition_matrix[(l - 1) * num_states + dd[i,l], dd[j,l]]
          }
        }
      }
      
      matrix_checkpoint(Mub)
      matrix_checkpoint(Pb)
      matrix_checkpoint(Pib)
      
      k1 <- (2*pi) ^ (-p/2)
      
      for (iteration in 1:max_iterations) {
        # FORWARD-BACKWARD %%% EXACT E STEP
        
        Gamma <- c();
        GammaX <- zeros(num_states*num_chains, p);
        Eta <- zeros(num_states*num_chains, num_states*num_chains); 
        Scale <- zeros(T, 1);
        Xi <- zeros(num_states*num_chains, num_states);

        # expand
        Mub <- zeros(num_states^num_chains, p);
        Pb <- ones(num_states^num_chains, num_states^num_chains);
        Pib <- ones(num_states^num_chains, 1);

        for (i in 1:(num_states^num_chains)) {
          
          for (l in 1:num_chains) {
            Mub[i,] <- as.matrix(Mub[i,] + mu[(l-1) * num_states + dd[i,l],])
            Pib[i,] <- as.matrix(Pib[i,] %*% priors[dd[i,l],l])
          }
          
          for (j in 1:(num_states^num_chains)) {
            for (l in 1:num_chains) {
              Pb[i,j] <- Pb[i,j] %*% transition_matrix[(l-1) * num_states + dd[i,l], dd[j,l]]
            }
          }
        }
        
        matrix_checkpoint(Mub)
        matrix_checkpoint(Pb)
        matrix_checkpoint(Pib)
        
        for (n in 1:N) {
          # print(sprintf("max_iterations = %d, N = %d", iteration, n))
          B <- zeros(T, num_states^num_chains)
          
          ### TODO: check we need to do a solve(covariance_matrix) here instead
          iCov <- ginv(covariance_matrix)
          k2 <- k1/sqrt(det(covariance_matrix))

          for (l in 1:(num_states^num_chains)) {
            d <- ones(T,1) %*% Mub[l,] - observations[((n-1) *T +1):(n*T),]
            B[,l] <- k2 * exp(-0.5 * rowSums((d %*% iCov) * d))
          }

          scale <- zeros(T, 1)
          alpha[1,] <- t(Pib) * B[1,]
          
          # forward probs
          scale[1] <- sum(alpha[1,])
          alpha[1,] <- alpha[1,] / scale[1]
          for (i in 2:T) {
            alpha[i,] <- (alpha[i-1,] %*% Pb) * B[i,]
            scale[i] <- sum(alpha[i,])
            alpha[i,] <- alpha[i,] / scale[i]
          }
         
          # backward probs
          beta[T,] <- ones(1, num_states^num_chains) / scale[T]
          for (i in (T-1):1) {
            beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(Pb)/ scale[i]
          }

          # joint prob
          gamma <- alpha * beta
          gamma <- gamma / rowSums(gamma)
          
          gamma1 <- gamma %*% collapse
          for (i in 1:T) {
            for (j in 1:num_chains) {
              d1 <- ((j-1) * num_states+1):(j*num_states)
              gamma1[i,d1] <- gamma1[i,d1] / sum(gamma1[i,d1])
            }
          }
          
          xi <- zeros(num_chains*num_states, num_states);
          for (i in 1:(T-1)) {
            tt <- t(alpha[i,] %*% collapse) %*% ((beta[i+1,] * B[i+1,]) %*% collapse)
            for (j in 1:num_chains) {
              d1 <- ((j-1)*num_states+1):(j*num_states)
              t2 <- transition_matrix[d1,] * tt[d1,d1]
              xi[d1,] <- as.matrix(xi[d1,] + t2/sum(t2))
            }
          }
          
          tt <- gamma %*% collapse2
          
          for (i in 1:T) {
            d1 <- ((i-1)*num_states*num_chains+1) : (i*num_states*num_chains)
            eta[d1,] <- matrix(tt[i,], nrow = num_states*num_chains)
            for (j in 1:num_chains) {
              d2 <- ((i-1)*num_states*num_chains+(j-1)*num_states+1) : ((i-1)*num_states*num_chains + j*num_states)
              for (l in 1:num_chains) {
                if (j==l) {
                  eta[d2, ((j-1)*num_states+1) : (j*num_states)] <- diag(gamma1[i, ((j-1)*num_states+1) : (j*num_states)])
                }
                else {
                  d3 <- sum(sum(eta[d2, ((l-1)*num_states+1):(l*num_states)]))
                  eta[d2, ((l-1)*num_states+1): (l*num_states)] <- eta[d2, ((l-1)*num_states+1):(l*num_states)] / d3
                }
              }
            }
            Eta <- Eta + eta[d1,]
            GammaX <- GammaX + as.matrix(gamma1[i,]) %*% observations[(n-1) * T+i,]
          }
          
          matrix_checkpoint(GammaX)
          
          Scale <- Scale + log(scale)
          Gamma = cbind(Gamma, gamma1)
          Xi <- Xi+xi
        }
        Eta <- (Eta +t(Eta))/2
        
        matrix_checkpoint(Eta)
        
        # Calculate Likelihood and determine convergence
        
        log_likelihood_prev <- log_likelihood
        log_likelihood <- sum(Scale)
        log_likelihood_trace <- c(log_likelihood_trace, log_likelihood)
        print(sprintf('iteration %i log likelihood = %f (base: %f, old: %f) %d', iteration, log_likelihood, log_likelihood_base, log_likelihood_prev,
                      (log_likelihood - log_likelihood_base) < (1 + tolerance) * (log_likelihood_prev-log_likelihood_base) ))
        
        #print(sprintf("%f  <  %f", log_likelihood - log_likelihood_base, (1 + tolerance) * (log_likelihood_prev-log_likelihood_base)))
        
        if (iteration <= 2) {
          log_likelihood_base <- log_likelihood
        } else if (log_likelihood < log_likelihood_prev-2) {
          print('Large likelihood violation \n')
        } else if (log_likelihood < log_likelihood_prev) {
          print('v')
        } else if ((log_likelihood - log_likelihood_base) < (1 + tolerance) * (log_likelihood_prev-log_likelihood_base) || !is.finite(log_likelihood)) {
          print('DONE')
          break
        }
        
        # num_chains STEP 
        
        # outputs -- using SVD as generally ill-conditioned (mu=pinv(Eta)*GammaX);
        
        mu <- pinv(Eta) %*% GammaX;
        
        matrix_checkpoint(Eta)
        matrix_checkpoint(pinv(Eta))
        #matrix_checkpoint(mu)
        
        # covariance
        covariance_matrix <- XX - as.matrix(t(GammaX)) %*% as.matrix(mu/(N*T))
  
        covariance_matrix <- (covariance_matrix + t(covariance_matrix))/2;
        covariance_matrix_det <- det(covariance_matrix);
        while (covariance_matrix_det <= 0) {
          print('Covariance problem')
          covariance_matrix <- covariance_matrix + eye(p) %*% (1.5*(-covariance_matrix_det)^(1/p)+.Machine$double.eps)
          covariance_matrix_det <- det(covariance_matrix);
        }
        
        chain_id <- 1
        #  transition matrix 
        for (i in 1:(num_states * num_chains)) {
          d1 <- sum(Xi[i,])
          if(d1 == 0) {
            transition_matrix[i,] <- ones(1,num_states)/num_states
          } else {
            transition_matrix[i,] <- Xi[i,]/d1
            # 4 chains
            # 3 states
            # print(paste("CHAIN:", chain_id))
            chain_id <- chain_id + 1
          }
        }

        #  priors
        priors <- zeros(num_states, num_chains)
        for (i in 1:N) {
          priors <- priors + matrix(Gamma[(i-1)*T+1, ], nrow = num_states)
        }
        priors <- priors/N
          
        # print(sprintf('computing num_chains step: %i\n\n', flops)        
      }
      
      result <- list(mu = mu,
                     priors = priors,
                     covariance_matrix = covariance_matrix,
                     transition_matrix = transition_matrix,
                     log_likelihood_trace = log_likelihood_trace)
    },
    
    forward_probs = function(observations, state_matrix, transition_matrix, emission_matrix) {
      trellis <-
        matrix(, nrow = nrow(state_matrix), ncol = nrow(observations))
      for (s in 1:nrow(state_matrix)) {
        trellis[s, 1] <-
          state_matrix$p[s] * emission_matrix[s, observations[1]]
      }
      
      for (t in 2:length(observations)) {
        for (current_state in 1:nrow(state_matrix)) {
          transition_probs <- rep(NA, nrow(state_matrix))
          for (prev_state in 1:nrow(state_matrix)) {
            transition_probs[prev_state] <-
              trellis[prev_state, t - 1] * transition_matrix[prev_state, current_state]
          }
          trellis[current_state, t] <-
          sum(transition_probs) * emission_matrix[current_state, observations[t]]
        }
      }
      trellis
    },
    
    backward_probs = function(observations) {
      trellis <- matrix(, nrow = nrow(state_matrix), ncol = length(dataset))
      for (s in 1:nrow(state_matrix)) {
        trellis[s, length(dataset)] <- 1
      }
      
      for (t in (length(dataset)):1) {
        for (current_state in 1:nrow(state_matrix)) {
          print(sprintf("T %d: ----------------------------------------------------------", t))
          transition_probs <- rep(NA, nrow(state_matrix))
          for (next_state in 1:nrow(state_matrix)) {
            transition_probs[next_state] <- (transition_matrix[current_state, next_state] * emission_matrix[next_state, dataset[t]] * trellis[next_state, t]) 
            print(sprintf("      transition_matrix(S%d | S%d) = %g, transition_matrix(%s | S%d) = %g, p = %g = %g", 
                          next_state, 
                          current_state, 
                          transition_matrix[current_state, next_state],
                          dataset[t], 
                          next_state, 
                          emission_matrix[next_state, dataset[t]],
                          trellis[next_state, t],
                          transition_probs[next_state]))
          }
          trellis[current_state, t-1] <- sum(transition_probs) #* emission_matrix[next_state, dataset[t+1]]
        }
      }
      trellis
    }
  )
  
  assign('this', this, envir = env)
  
  class(this) <- append(class(this), "FactorialHmm_OLD")
  
  return(this)
}

FactorialHmm <- function()
{
  env <- environment()
  
  alpha <- NA
  
  this <- list(
    env = env,
    
    # get/set value
    get_value = function(key) {
      return(get(key, envir = env))
    },
    set_value = function(key, value) {
      return(assign(key, value, envir = env))
    },
    
    fit = function(observations, num_chains = 2, num_states = num_chains, T = NULL, max_iterations = 100, tolerance = 0.0001) {
      message("FIT")
    },
    
    sample = function(seq_length) {
      message("SAMPLE")
    }
  )
    
  assign('this', this, envir = env)
  
  class(this) <- append(class(this), "FactorialHmm")
  
  return(this)
}

  
dataset <- as.matrix(read.csv("./data/m200x3.csv", header = FALSE))
#dataset <- as.matrix(read.csv("./data/X", header = FALSE))

#dataset <- rand_matrix(1:100, 200, 3)

fhmm <- FactorialHmm()

result <- fhmm$fit(dataset, num_chains = 1, num_states = 4, max_iterations = 121)
print(result)

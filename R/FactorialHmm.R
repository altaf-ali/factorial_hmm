source("R/debug.R")
source("R/matlab.R")

FactorialHmm <- setRefClass("FactorialHmm",
  contains = "DebugTrace",
  fields = list(
    observations = "data.frame",
    num_chains = "numeric",
    num_states = "numeric",
    max_iterations = "numeric",
    tolerance = "numeric",
    sequence_length = "numeric",
    rng = "ANY" # random number generator
  ),
  
  methods = list(
    initialize = function(num_chains = 1, num_states = 2, max_iterations = 100) {
      dbg_trace(num_chains, num_states)
      
      num_chains <<- num_chains
      num_states <<- num_states
      max_iterations <<- max_iterations
      tolerance <<- 0.0001
      
      sequence_length <<- 0
      
      rng <<- MatlabRandomNumberGenerator$new()
    },
    
    finalize = function() {
      rng$close()
      rng <<- NULL  
    },

    fit = function(observations) {

      dbg_trace(dim(observations))
      rng$set.seed(666)
      
      p <- ncol(observations)
      sequence_length <<- nrow(observations)
      
      N <- sequence_length
      N <- N / sequence_length

      covariance_matrix <- cov(observations)
      if (length(covariance_matrix) > 1) {
        covariance_matrix <- diag(diag(covariance_matrix))
      }
      
      XX <- t(observations) %*% observations / (N * sequence_length)
      
      rand_matrix <- rng$randn(num_chains * num_states, p)
      mu <- ((rand_matrix %*% sqrt(covariance_matrix)) / num_chains) + (ones(num_states * num_chains, 1) %*% (colMeans(observations) / num_chains));
      
      priors <- rng$runif(num_states, num_chains)
      priors <- cdiv(priors, colSums(priors))
      
      transition_matrix <- rng$runif(num_states * num_chains, num_states)
      transition_matrix <- transition_matrix / rowSums(transition_matrix)
      
      log_likelihood_trace <- c()
      log_likelihood <- 0
      log_likelihood_base <- 0
      
      dd <- zeros(num_states ^ num_chains, num_chains)
      for (i in 1:(num_states ^ num_chains)) {
        dd[i,] <- base(i - 1, num_states, num_chains)
      }
      
      eta <- zeros(sequence_length * num_states * num_chains,num_states * num_chains)
      
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
      
      k1 <- (2*pi) ^ (-p/2)
      
      for (iteration in 1:max_iterations) {
        # FORWARD-BACKWARD %%% EXACT E STEP
        
        Gamma <- c();
        GammaX <- zeros(num_states*num_chains, p);
        Eta <- zeros(num_states*num_chains, num_states*num_chains); 
        Scale <- zeros(sequence_length, 1);
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
        
        for (n in 1:N) {
          B <- zeros(sequence_length, num_states^num_chains)
          
          ### TODO: check we need to do a solve(covariance_matrix) here instead
          iCov <- MASS::ginv(covariance_matrix)
          k2 <- k1/sqrt(det(covariance_matrix))
          
          for (l in 1:(num_states^num_chains)) {
            d <- ones(sequence_length,1) %*% Mub[l,] - observations[((n-1) *sequence_length +1):(n*sequence_length),]
            B[,l] <- k2 * exp(-0.5 * rowSums((d %*% iCov) * d))
          }
          
          # forward probs
          probs <- forward_probs(B, Pb, Pib)
          scale <- probs$scale
          alpha <- probs$alpha
          
          # backward probs
          beta <- backward_probs(B, Pb, scale)

          # joint prob
          gamma <- alpha * beta
          gamma <- gamma / rowSums(gamma)
          
          gamma1 <- gamma %*% collapse
          for (i in 1:sequence_length) {
            for (j in 1:num_chains) {
              d1 <- ((j-1) * num_states+1):(j*num_states)
              gamma1[i,d1] <- gamma1[i,d1] / sum(gamma1[i,d1])
            }
          }
          
          xi <- zeros(num_chains*num_states, num_states);
          for (i in 1:(sequence_length-1)) {
            tt <- t(alpha[i,] %*% collapse) %*% ((beta[i+1,] * B[i+1,]) %*% collapse)
            for (j in 1:num_chains) {
              d1 <- ((j-1)*num_states+1):(j*num_states)
              t2 <- transition_matrix[d1,] * tt[d1,d1]
              xi[d1,] <- as.matrix(xi[d1,] + t2/sum(t2))
            }
          }
          
          tt <- gamma %*% collapse2
          
          for (i in 1:sequence_length) {
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
            GammaX <- GammaX + as.matrix(gamma1[i,]) %*% observations[(n-1) * sequence_length+i,]
          }
          
          Scale <- Scale + log(scale)
          Gamma = cbind(Gamma, gamma1)
          Xi <- Xi+xi
        }
        Eta <- (Eta +t(Eta))/2
        
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
        
        # M STEP 
        
        # outputs -- using SVD as generally ill-conditioned (mu=pinv(Eta)*GammaX);
        
        mu <- pinv(Eta) %*% GammaX;
        
        # covariance
        covariance_matrix <- XX - as.matrix(t(GammaX)) %*% as.matrix(mu/(N*sequence_length))
        
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
          priors <- priors + matrix(Gamma[(i-1)*sequence_length+1, ], nrow = num_states)
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
   
    forward_probs = function(B, Pb, Pib)  {
      alpha <- zeros(sequence_length, num_states ^ num_chains)
      scale <- zeros(sequence_length, 1)
      alpha[1,] <- t(Pib) * B[1,]
      
      scale[1] <- sum(alpha[1,])
      alpha[1,] <- alpha[1,] / scale[1]
      for (i in 2:sequence_length) {
        alpha[i,] <- (alpha[i-1,] %*% Pb) * B[i,]
        scale[i] <- sum(alpha[i,])
        alpha[i,] <- alpha[i,] / scale[i]
      }
      
      return(list(scale = scale, alpha = alpha))
    },

    backward_probs = function(B, Pb, scale) {
      beta <- zeros(sequence_length, num_states ^ num_chains)
      beta[sequence_length,] <- ones(1, num_states ^ num_chains) / scale[sequence_length]
      for (i in (sequence_length-1):1) {
        beta[i,] <- (beta[i+1,] * B[i+1,]) %*% t(Pb)/ scale[i]
      }
      
      return(beta)
    },
    
    sample = function() {
      dbg_trace()
    }
  )
)

dataset <- as.matrix(read.csv("./data/m200x3.csv", header = FALSE))
#dataset <- as.matrix(read.csv("./data/X", header = FALSE))

model <- FactorialHmm$new(num_chains = 2, num_states = 4)
model$fit(dataset)

model$rng$shutdown()



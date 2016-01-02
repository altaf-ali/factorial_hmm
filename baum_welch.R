# 
# p_s <- .85
# p_t <- .15
# 
# p_As <- .4
# p_Bs <- .6
# 
# p_At <- .5
# p_Bt <- .5
# 
# alpha_ABBA_1_s <- p_s * p_As 
# alpha_ABBA_1_t <- p_t * p_At
# 
# alpha_ABBA_2_s <- alpha_ABBA_1_s p_s * p_As 
# alpha_ABBA_2_t <- p_t * p_At
# 
# 
# transition_matrix <- data.frame(S1 = c(1, .5, 0), S2 = c(0, .5, .5), SF = c(0, 0, .5),
#                                 row.names = c("S0", "S1", "S2"))
# emission_matrix <- data.frame(R = c(0.33, 0.33), W = c(0.33, 0.33), B = c(0.33, 0.33),
#                               row.names = c("S1", "S2"))
# 
# # A0,1 * Pr(R|S1) * A1,1 * Pr(W|S1) * A1,1 * Pr(B|S1) * A1,2 * Pr(B|S2) * A2,F
# 
# 1 * 0.33 * 0.5 * 0.33 * 0.5 * 0.33 * 0.5 * 0.33 * 0.5

# P1 <- transition_matrix["S0", "S1"] * emission_matrix["S1", "R"] *
#   transition_matrix["S1", "S1"] * emission_matrix["S1", "W"] *
#   transition_matrix["S1", "S1"] * emission_matrix["S1", "B"] *
#   transition_matrix["S1", "S2"] * emission_matrix["S1", "B"] *
#   transition_matrix["S2", "SF"] 
#   


# > transition_matrix
#        S1   S2
# S1    0.6   0.4
# S2    0.3   0.7

# > emission_matrix
#           R     W     B
# S1      0.3   0.4   0.3
# S2      0.4   0.3   0.3

# > state_matrix
#       p
# S1    0.8
# S2    0.2

log.sum.exp<- function(x) {
  # Computes log(sum(exp(x))
  # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}

dataset = c("R", "W", "B", "B")
#dataset = c("START", dataset, "STOP")
#dataset = c("R", "W")

state_matrix <- data.frame(p = c(.8, .2), row.names = c("S1", "S2"))
transition_matrix <- data.frame(S1 = c(.6, .3), S2 = c(.4, .7), row.names = c("S1", "S2"))
emission_matrix <- data.frame(R = c(.3, .4), W = c(.4, .3), B = c(.3, .3), row.names = c("S1", "S2"))

forward_probs <- function(dataset, state_matrix, transition_matrix, emission_matrix) {
  trellis <- matrix(, nrow = nrow(state_matrix), ncol = length(dataset))
  for (s in 1:nrow(state_matrix)) {
    trellis[s, 1] <- state_matrix$p[s] * emission_matrix[s, dataset[1]]
  }
  
  for (t in 2:length(dataset)) {
    for (current_state in 1:nrow(state_matrix)) {
      transition_probs <- rep(NA, nrow(state_matrix))
      for (prev_state in 1:nrow(state_matrix)) {
        transition_probs[prev_state] <- trellis[prev_state, t-1] * transition_matrix[prev_state, current_state]
      }
      trellis[current_state, t] <- sum(transition_probs) * emission_matrix[current_state, dataset[t]]
    }
  }
  trellis
}
forward_probs(dataset, state_matrix, transition_matrix, emission_matrix)

backward_probs <- function(dataset, state_matrix, transition_matrix, emission_matrix) {
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
        print(sprintf("      P(S%d | S%d) = %g, P(%s | S%d) = %g, p = %g = %g", 
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
backward_probs(dataset, state_matrix, transition_matrix, emission_matrix)

X = read.csv("./data/X", header = FALSE)
#cov(X)
#Cov=diag(diag(cov(X))); 
#XX=X'*X/(N*T);

max_iterations <- 2
for (i in 1:max_iterations) {
  print(sprintf("iteration %d", i))
  
  #Gamma=[];
  #GammaX=zeros(M*K,p);
  #Eta=zeros(K*M,K*M); 
  #Xi=zeros(M*K,K);
}


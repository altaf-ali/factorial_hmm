# trans = [0.95,0.05;
#          0.10,0.90];
# emis = [1/6 1/6 1/6 1/6 1/6 1/6;
#         1/10 1/10 1/10 1/10 1/10 1/2];
# 
# [seq,states] = hmmgenerate(100,trans,emis);
# pStates = hmmdecode(seq,trans,emis);
# [seq,states] = hmmgenerate(100,trans,emis,...
#                            'Symbols',{'one','two','three','four','five','six'})
# pStates = hmmdecode(seq,trans,emis,...
#                     'Symbols',{'one','two','three','four','five','six'});
# 
# x <- c(1,2,1,1,3,4,4,1,2,4,1,4,3,4,4,4,3,1,3,2,3,3,3,4,2,2,3)
# p <- matrix(nrow = 4, ncol = 4, 0)
# for (t in 1:(length(x) - 1)) p[x[t], x[t + 1]] <- p[x[t], x[t + 1]] + 1
# for (i in 1:4) p[i, ] <- p[i, ] / sum(p[i, ])

nSim          = 10000
States        = c("One","Two")
Symbols       = 1:6
transProbs    = matrix(c(.99,.01,.02,.98), c(length(States),length(States)), byrow=TRUE)

transProbs    = matrix(c(.75,.25,.4,.6), c(length(States),length(States)), byrow=TRUE)

emissionProbs = matrix(c(rep(1/6,6),c(rep(.1,5),.5)), c(length(States),length(Symbols)), byrow=TRUE)
hmm = initHMM(States, Symbols, transProbs=transProbs, emissionProbs=emissionProbs)
sim = simHMM(hmm,nSim)

baumWelch(hmm, sim[[2]])


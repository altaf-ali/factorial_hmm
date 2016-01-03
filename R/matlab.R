RandomNumberGenerator <- setRefClass("RandomNumberGenerator",
   methods = list(
     set.seed = function(seed, kind = NULL) { base::set.seed(seed, kind) },
     runif = function(rows, cols, min = 0, max = 1) { matrix(base::runif(rows * cols, min, max), rows, cols) },
     randn = function(rows, cols, mean = 0, sd = 1) { matrix(base::rnorm(rows * cols, mean, sd), rows, cols) },
     shutdown = function() { }
   )
)

MatlabRandomNumberGenerator <- setRefClass("MatlabRandomNumberGenerator", 
  contains = "RandomNumberGenerator",
  fields = list(
    connection = "ANY"
  ),
 
  methods = list(
    initialize = function() { 
      connection <<- NULL
      R.matlab::Matlab$startServer(matlab = '/Applications/MATLAB_R2015b.app/bin/matlab')

      matlab <- R.matlab::Matlab()
      status <- open(matlab)
      if (status)
        connection <<- matlab
    },
    
    finalize = function() { shutdown() },
     
    shutdown = function() {
      if (is.null(connection))
        return()
      
      R.matlab::close.Matlab(connection)
      connection <<- NULL
    },
    
    set.seed = function(seed) { R.matlab::evaluate(connection, sprintf("rng(%s)", seed)) },
     
    runif = function(rows, cols) { .rand("rand", rows, cols) },
    randn = function(rows, cols) { .rand("randn", rows, cols) },
     
    .rand = function(func, rows, cols) {
      R.matlab::evaluate(connection, sprintf("x = %s(%s, %s)", func, rows, cols)) 
      R.matlab::getVariable(connection, "x")$x
    }
  )
)

# matlab compatiliblity functions
ones <- function(nrow = 1, ncol = 1) { matrix(1, nrow, ncol) }
zeros <- function(nrow = 1, ncol = 1) { matrix(0, nrow, ncol) }
eye <- function(d) { diag(1, d) }
cdiv <- function(x, y) { t(t(x)/y) }

base <- function(k, m, d) {
  mm <- m ^ c((d - 1):0)
  v <- rep(NA, d)
  
  for (i in 1:d) {
    v[i] <- trunc(k / mm[i])
    k <- k - mm[i] * v[i]
  }
  return(v + 1)
}

pinv <- function(x) {
  #ginv(x)
  #pracma::pinv(x)
  corpcor::pseudoinverse(x)
}



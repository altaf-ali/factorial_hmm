# MATLAB server

#library(R.matlab)

MaltlabConnection <- function()
{
  env <- environment()
  
  this <- list(
    env = env,
    
    # get/set value
    get_value = function(key) { 
      if (exists(key, envir = env))
        return(get(key, envir = env))
      else 
        return()
    },
    set_value = function(key, value) { return(assign(key, value, envir = env)) },
    
    # server start
    start_server = function() {
      R.matlab::Matlab$startServer(matlab = '/Applications/MATLAB_R2015b.app/bin/matlab')
    },
    
    # connection management
    connect = function() { 
      connection <- R.matlab::Matlab()
      status <- open(connection)
      if (status)
        this$set_value("connection", connection)
    },
    
    get_connection = function() { this$get_value("connection") },

    # random numbers
    set.seed = function(seed) { R.matlab::evaluate(this$get_connection(), sprintf("rng(%s)", seed)) },
    
    # set format
    format = function(fmt) { R.matlab::evaluate(this$get_connection(), sprintf("format %s", fmt)) },
    
    get_rand = function(func, rows, cols) {
      R.matlab::evaluate(this$get_connection(), sprintf("x = %s(%s, %s)", func, rows, cols)) 
      R.matlab::getVariable(this$get_connection(), "x")$x
    },
    
    rand = function(rows, cols) { this$get_rand("rand", rows, cols) },
    randn = function(rows, cols) { this$get_rand("randn", rows, cols) }
  )
    
  assign('this', this, envir = env)
  
  class(this) <- append(class(this), "MaltlabConnection")
  
  return(this)
}

if (FALSE) {
  matlab <- MaltlabConnection()
  matlab$start_server()
  matlab$connect()
  matlab$set.seed(666)
  x = matlab$rand(5,3)
  print(x)
}

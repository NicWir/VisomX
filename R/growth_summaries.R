summary.drFitSpline <- function (object,...)
  {
    # object of class drFitSpline
    data.frame(object$parameters)
}

summary.gcFitSpline <- function(object,...)
  {
    
    # object of class gcFitSpline
    
    contents.fitted.spline  <- c("mu.spline", "lambda.spline", "A.spline", "integral.spline")
    
    if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){
      table<-rep(NA,length(contents.fitted.spline))
    }
    else{
      table <- c(object$parameters$mu, object$parameters$lambda,  object$parameters$A, object$parameters$integral)
    }
    
    table               <- data.frame(t(table))
    colnames(table)     <- contents.fitted.spline
    summary.gcFitSpline <- table
}

summary.gcFitModel <- function(object, ...)
  {
    # object of class gcFitModel
    
    contents.fitted.param     = c("mu.model", "lambda.model", "A.model", "integral.model",
                                  "stdmu.model", "stdlambda.model", "stdA.model",
                                  "ci90.mu.model.lo", "ci90.mu.model.up",
                                  "ci90.lambda.model.lo", "ci90.lambda.model.up",
                                  "ci90.A.model.lo", "ci90.A.model.up",
                                  "ci95.mu.model.lo", "ci95.mu.model.up",
                                  "ci95.lambda.model.lo", "ci95.lambda.model.up",
                                  "ci95.A.model.lo", "ci95.A.model.up")
    
    
    if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){
      
      
      table<-rep(NA,length(contents.fitted.param))
    }
    else{
      table <- c(object$parameter$mu[1], object$parameter$lambda[1],  object$parameter$A[1], object$parameter$integral,
                 object$parameter$mu[2], object$parameter$lambda[2],  object$parameter$A[2], 
                 object$parameter$mu[1]-1.645*object$parameter$mu[2], object$parameter$mu[1]+1.645*object$parameter$mu[2],
                 object$parameter$lambda[1]-1.645*object$parameter$lambda[2], object$parameter$lambda[1]+1.645*object$parameter$lambda[2],
                 object$parameter$A[1]-1.645*object$parameter$A[2], object$parameter$A[1]+1.645*object$parameter$A[2],
                 object$parameter$mu[1]-1.96*object$parameter$mu[2], object$parameter$mu[1]+1.96*object$parameter$mu[2],
                 object$parameter$lambda[1]-1.96*object$parameter$lambda[2], object$parameter$lambda[1]+1.96*object$parameter$lambda[2],
                 object$parameter$A[1]-1.96*object$parameter$A[2], object$parameter$A[1]+1.96*object$parameter$A[2])
      
    }
    table <- data.frame(t(table))
    colnames(table) <- contents.fitted.param
    summary.gcFitModel <- table
}

summary.drFit <- function(object, ...)
  {
    # object of class drFit
    data.frame(object$drTable)
}

summary.gcBootSpline <- function(object, ...)
  {
    # object of class gcBootSpline
    contents.bootstrap        <- c("mu.bt", "lambda.bt", "A.bt", "integral.bt", "stdmu.bt", "stdlambda.bt", "stdA.bt", "stdintegral.bt",
                                   "ci90.mu.bt.lo", "ci90.mu.bt.up", "ci90.lambda.bt.lo", "ci90.lambda.bt.up",
                                   "ci90.A.bt.lo", "ci90.A.bt.up", "ci90.integral.bt.lo", "ci90.integral.bt.up",
                                   "ci95.mu.bt.lo", "ci95.mu.bt.up", "ci95.lambda.bt.lo", "ci95.lambda.bt.up",
                                   "ci95.A.bt.lo", "ci95.A.bt.up", "ci95.integral.bt.lo", "ci95.integral.bt.up")
    
    
    if (object$bootFlag==FALSE){
      table<-rep(NA,length(contents.bootstrap))
    }
    else{
      mu          <- mean(object$mu, na.rm=TRUE)
      lambda      <- mean(object$lambda, na.rm=TRUE)
      A           <- mean(object$A, na.rm=TRUE)
      integral    <- mean(object$integral, na.rm=TRUE)
      
      mu.sd       <- sd(object$mu, na.rm=TRUE)
      lambda.sd   <- sd(object$lambda, na.rm=TRUE)
      A.sd        <- sd(object$A, na.rm=TRUE)
      integral.sd <- sd(object$integral, na.rm=TRUE)
      
      table <- c(mu, lambda, A, integral, mu.sd, lambda.sd, A.sd, integral.sd,
                 mu-1.645*mu.sd, mu+1.645*mu.sd, lambda-1.645*lambda.sd,     lambda+1.645*lambda.sd,
                 A-1.645*A.sd,   A+1.645*A.sd,   integral-1.645*integral.sd, integral+1.645*integral.sd,
                 mu-1.96*mu.sd,  mu+1.96*mu.sd,  lambda-1.96*lambda.sd,      lambda+1.96*lambda.sd,
                 A-1.96*A.sd,    A+1.96*A.sd,    integral-1.96*integral.sd,  integral+1.96*integral.sd)
    }
    
    table               <- data.frame(t(table))
    colnames(table)     <- contents.bootstrap
    summary.gcBootSpline <- table
}

summary.drBootSpline <- function(object, ...)
  {
    # object of class drBootSpline
    contents.bootstrap <- c("meanEC50", "sdEC50", "ci90EC50.lo", "ci90EC50.up", "ci95EC50.lo", "ci95EC50.up", 
                            "meanEC50.orig", "ci90EC50.orig.lo", "ci90EC50.orig.up", "ci95EC50.orig.lo", "ci95EC50.orig.up")
    if (object$bootFlag==FALSE){
      table                <- rep(NA,length(contents.bootstrap))
    }
    else{
      m.test <- mean(object$ec50.boot, na.rm=TRUE)
      s.test <- sd(object$ec50.boot, na.rm=TRUE)	
      EC50   <- c(m.test, s.test, m.test-1.645*s.test, m.test+1.645*s.test, m.test-1.96*s.test, m.test+1.96*s.test)
      if (object$control$log.x.dr==TRUE){
        EC50.orig <- c(exp(m.test)-1, exp(m.test-1.645*s.test)-1, exp(m.test+1.645*s.test)-1, exp(m.test-1.96*s.test)-1, exp(m.test+1.96*s.test)-1)
      }
      else
      {
        EC50.orig <- c(m.test, m.test-1.645*s.test, m.test+1.645*s.test, m.test-1.96*s.test, m.test+1.96*s.test)
      }
      
      table <- c(EC50, EC50.orig)
    }
    
    table                <- data.frame(t(table))
    colnames(table)      <- contents.bootstrap
    summary.drBootSpline <- table
}

summary.gcFit <- function(object,...)
  {
    # object of class gcFit
    data.frame(object$gcTable)
}

summary.gcFitLinear <- function(object,...)
  {
  # object of class gcFitLinear
  
  contents.fitted.param     = c("mu.linfit", "lambda.linfit",
                                "stdmu.linfit", "tmu.start.linfit", 
                                "tmu.end.linfit", "r2mu.linfit")
  
  
  if ((is.na(object$fitFlag)==TRUE)|(object$fitFlag==FALSE)){
    table<-rep(NA,length(contents.fitted.param))
  }
  else{
    table <- c(object$par[3], object$par[5],  
               object$par[4], object$par[6],
               object$par[7], object$rsquared)
    
  }
  table <- data.frame(t(table))
  colnames(table) <- contents.fitted.param
  summary.gcFitLinear <- table
}
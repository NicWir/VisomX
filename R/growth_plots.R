#' @rdname plot
#'
setMethod("plot", c("gcFitLinear", "missing"),
          function(x, y, log="", which=c("fit", "diagnostics"), ...) 
            {
            which <- match.arg(which)
            
            switch(which,
                   fit = {
                     
                     plot(x$"raw.data" ~ x$"raw.time", xlab="Time", ylab="Ln(density)",
                          log=log, las=1, main = "Exponential fit", ...)
                     points(x$raw.data[x$ndx] ~ x$raw.time[x$ndx], pch=21, col="black", bg="red")
                     
                     ## lag phase
                     lag <- x$par["lag"]
                     
                     time <- seq(lag, max(x$"raw.time"), length=200)
                     coef_ <- x$par
                     lines(time, x$FUN(time, c(y0=unname(coef_["y0_lm"]), mumax=unname(coef_["mumax"])))[,"y"], lty=2, lwd=2, col=ggplot2::alpha("firebrick3", 0.7), ...)
                     lines(c(min(x$"raw.time"[1]), lag), rep(x$"raw.data"[1], 2), lty=2, lwd=2, col=ggplot2::alpha("firebrick3", 0.7))
                   },
                   diagnostics = {
                     opar <- par(no.readonly = TRUE)
                     on.exit(par(opar))
                     par(mfrow=c(1,2))
                     
                     ## residuals vs. fitted
                     obs <- obs(x)
                     sim <- x$FUN(x$"raw.time", x$par)
                     plot(fit.linear[["fit"]][["residuals"]] ~ fitted(fit.linear[["fit"]]), xlab="fitted", ylab="residuals")
                     abline(h=0, col="grey")
                     ## normal q-q-plot
                     qqnorm(fit.linear[["fit"]][["residuals"]])
                     qqline(fit.linear[["fit"]][["residuals"]])
                   }
            )
          }
)

plot.gcFitModel <- function(x, add=FALSE, raw=TRUE, slope=TRUE, pch=1, colData=1, colModel=1, colLag = 1, cex=1, ...)
  {
    # x an object of class gcFitModel
    
    # /// check input parameters
    if (is.logical(add)==FALSE)   stop("Need logical value for: add")
    if (is.logical(raw)==FALSE)   stop("Need logical value for: raw")
    if (is.logical(slope)==FALSE) stop("Need logical value for: slope")
    if (is.numeric(pch)==FALSE)   stop("Need numeric value for: pch")
    if (is.numeric(cex)==FALSE)   stop("Need numeric value for: cex")
    

    # /// check if a data fit is available
    if ((is.na(x$fitFlag)==TRUE)|(x$fitFlag==FALSE)){
      warning("plot.gcFitModel: no data fit available!")
    }
    else{
      if (raw==TRUE){
        if (add==TRUE){
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==FALSE)){
            try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch,cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colModel, type="l") )
          }
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==TRUE)){
            try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colModel, type="l") )
          }
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==FALSE)){
            try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colModel, type="l" ) )
          }
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==TRUE)){
            try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colModel, type="l") )
          }
        }
        else{ # of if (add==TRUE){
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==FALSE)){
            try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==TRUE)){
            try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="Growth [Ln(y(t)/y0)]", col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==FALSE)){
            try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l" ) )
          }
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==TRUE)){
            try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="Growth [Ln(y(t)/y0)]", col=colData, pch=pch, cex=cex) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
        }
      }
      else{ # of if (raw==TRUE)
        if (add==TRUE){
          # /// try to plot data fit
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==FALSE)){
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==TRUE)){
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==FALSE)){
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==TRUE)){
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
        }
        else{ # of if (add==TRUE)
          # /// try to plot data fit
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==FALSE)){
            try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", col=0) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colModel, type="l") )
          }
          
          if ((x$control$log.x.gc==FALSE) && (x$control$log.y.model==TRUE)){
            try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="Growth [Ln(y(t)/y0)]", col=0) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
          
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==FALSE)){
            try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", col=0 ) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l" ) )
          }
          
          if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.model==TRUE)){
            try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="Growth [Ln(y(t)/y0)]", col=0) )
            try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colModel, type="l") )
          }
        }
      }
      
      # /// add tangent at maximum slope
      if (slope==TRUE){
        mu     <- as.numeric(x$parameter$mu[1])
        lambda <- as.numeric(x$parameter$lambda[1])
        bla    <- (x$fit.time)*mu
        bla    <- bla+(-x$parameter$mu[1]*x$parameter$lambda[1])
        try(lines(x$fit.time, bla, lw=2, lty=2, col=colModel))
      }
    }
  }

plot.drBootSpline <- function (x,
                               pch = 1,
                               colData = 1,
                               colSpline = 1,
                               cex = 0.5,
                               ...)
{
  # x an object of class drBootSpline
  
  # /// initialize "Empty Plot" function
  empty.plot  <- function(text = "Empty plot", main = "") {
    plot(
      c(0, 1, 0, 1, 0),
      c(0, 1, 1, 0, 0),
      type = "l",
      axes = FALSE,
      xlab = "",
      ylab = "",
      lwd = 1,
      col = "gray",
      main = main
    )
    lines(c(0, 0),
          c(0, 1),
          type = "l",
          lwd = 1,
          col = "gray")
    lines(c(1, 1),
          c(1, 0),
          type = "l",
          lwd = 1,
          col = "gray")
    text(0.5, 0.1, text, col = "gray")
  }
  
  # /// check input parameters
  if (FALSE %in% (colData %in% c(colors(), 0:8)))
    stop("colData needs to be numeric from 0:8 or a string from colors()")
  if (FALSE %in% (colSpline %in% c(colors(), 0:8)))
    stop("colSpline needs to be numeric from 0:8 or a string from colors()")
  if (is.numeric(pch) == FALSE)
    stop("Need numeric value for: pch")
  if (is.numeric(cex) == FALSE)
    stop("Need numeric value for: cex")
  
  if (x$bootFlag == FALSE) {
    empty.plot()
  }
  else{
    colSpline   <-
      rep(colSpline, (x$control$nboot.dr %/% length(colSpline)) + 1)
    conc.log    <- log(x$raw.conc + 1)
    test.log    <- log(x$raw.test + 1)
    conc        <- x$raw.conc
    test        <- x$raw.test
    
    global.minx <- min(min(x$boot.conc))
    global.maxx <- max(max(x$boot.conc))
    global.miny <- min(min(x$boot.test))
    global.maxy <- max(max(x$boot.test))
    
    dev.new()
    # initialize plot
    if ((x$control$log.x.dr == TRUE) &&
        (x$control$log.y.dr == FALSE)) {
      plot(
        c(global.minx, global.maxx),
        c(global.miny, global.maxy),
        type = "n",
        xlab = "ln(1+concentration)",
        ylab = "response"
      )
    }
    else{
      if ((x$control$log.x.dr == FALSE) &&
          (x$control$log.y.dr == FALSE)) {
        plot(
          c(global.minx, global.maxx),
          c(global.miny, global.maxy),
          type = "n",
          xlab = "concentration",
          ylab = "response"
        )
      }
      else{
        if ((x$control$log.x.dr == TRUE) && (x$control$log.y.dr == TRUE)) {
          plot(
            c(global.minx, global.maxx),
            c(global.miny, global.maxy),
            type = "n",
            xlab = "ln(1+concentration)",
            ylab = "ln(1+response)"
          )
        }
        else{
          if ((x$control$log.x.dr == FALSE) && (x$control$log.y.dr == TRUE)) {
            plot(
              c(global.minx, global.maxx),
              c(global.miny, global.maxy),
              type = "n",
              xlab = "concentration",
              ylab = "ln(1+response)"
            )
          }
        }
      }
    }
    
    # /// plot raw data
    points(
      x$raw.conc,
      x$raw.test,
      col = colData,
      pch = pch,
      cex = cex
    )
    
    # /// loop over all fitted splines and plot drFitSpline objects
    for (i in 1:x$control$nboot.dr) {
      plot(
        x$boot.drSpline[[i]],
        add = TRUE,
        ec50line = FALSE,
        pch = 0,
        colSpline = colSpline[i],
        colData = 0,
        cex = cex
      )
    }
    
    dev.new()
    if (sum(!is.na(x$ec50.boot)) == length(x$ec50.boot)) {
      hist(
        x$ec50.boot,
        col = "gray",
        main = as.character(x$gcID),
        xlab = "EC50"
      )
    }
    else{
      empty.plot()
    }
    
  } # /// of if (x$bootFlag==FALSE){
  
}

plot.drFit <- function(x, ...)
{
  # x an object of class drFit
  
  n <- length(x$drFittedSplines)
  
  # /// plot all drFitSpline objects
  for (i in 1:n) {
    #x11()
    try(plot(x$drFittedSplines[[i]]))
    title(as.character(x$drFittedSplines[[i]]$drID))
  }
  
}

plot.drFitSpline <-
  function (x,
            add = FALSE,
            ec50line = TRUE,
            pch = 1,
            colSpline = 1,
            colData = 1,
            cex = 1,
            ...)
  {
    # x an object of class drFitSpline
    
    # /// check input parameters
    if (FALSE %in% (colData %in% c(colors(), 0:8)))
      stop("colData needs to be numeric from 0:8 or a string from colors()")
    if (FALSE %in% (colSpline %in% c(colors(), 0:8)))
      stop("colSpline needs to be numeric from 0:8 or a string from colors()")
    
    if (is.logical(add) == FALSE)
      stop("Need logical value for: add")
    if (is.logical(ec50line) == FALSE)
      stop("Need logical value for: ec50line")
    if (is.numeric(pch) == FALSE)
      stop("Need numeric value for: pch")
    if (is.numeric(cex) == FALSE)
      stop("Need numeric value for: cex")
    
    if (add == FALSE) {
      if ((x$control$log.x.dr == TRUE) && (x$control$log.y.dr == TRUE)) {
        plot(
          log(x$raw.conc + 1),
          log(x$raw.test + 1),
          pch = pch,
          cex = cex,
          col = colData,
          xlab = "ln(1+concentration)",
          ylab = "ln(1+response)"
        )
      }
      else
      {
        if ((x$control$log.x.dr == FALSE) && (x$control$log.y.dr == TRUE)) {
          plot(
            x$raw.conc,
            log(x$raw.test + 1),
            pch = pch,
            cex = cex,
            col = colData,
            xlab = "concentration",
            ylab = "ln(1+response)"
          )
        }
        else
        {
          if ((x$control$log.x.dr == TRUE) && (x$control$log.y.dr == FALSE)) {
            plot(
              log(x$raw.conc + 1),
              x$raw.test,
              pch = pch,
              cex = cex,
              col = colData,
              xlab = "ln(1+concentration)",
              ylab = "response"
            )
          }
          else
          {
            if ((x$control$log.x.dr == FALSE) && (x$control$log.y.dr == FALSE)) {
              plot(
                x$raw.conc,
                x$raw.test,
                pch = pch,
                cex = cex,
                col = colData,
                xlab = "concentration",
                ylab = "response"
              )
            }
          }
        }
      }
    }
    else{
      if ((x$control$log.x.dr == TRUE) && (x$control$log.y.dr == TRUE)) {
        points(
          log(x$raw.conc + 1),
          log(x$raw.test + 1),
          pch = pch,
          cex = cex,
          col = colData
        )
      }
      else
      {
        if ((x$control$log.x.dr == FALSE) && (x$control$log.y.dr == TRUE)) {
          points(
            x$raw.conc,
            log(x$raw.test + 1),
            pch = pch,
            cex = cex,
            col = colData
          )
        }
        else
        {
          if ((x$control$log.x.dr == TRUE) && (x$control$log.y.dr == FALSE)) {
            points(
              log(x$raw.conc + 1),
              x$raw.test,
              pch = pch,
              cex = cex,
              col = colData
            )
          }
          else
          {
            if ((x$control$log.x.dr == FALSE) && (x$control$log.y.dr == FALSE)) {
              points(
                x$raw.conc,
                x$raw.test,
                pch = pch,
                cex = cex,
                col = colData
              )
            }
          }
        }
      }
    }
    
    try(lines(
      x$fit.conc,
      x$fit.test,
      type = "l",
      lwd = 2,
      col = colSpline
    ))
    
    if (ec50line == TRUE) {
      #vertical lines
      totmin = min(min(x$fit.conc), min(x$fit.test))
      lines(c(x$parameters$EC50, x$parameters$EC50),
            c(totmin - 1, x$parameters$yEC50),
            lty = 2)
      #horizontal
      lines(c(-1, x$parameters$EC50),
            c(x$parameters$yEC50, x$parameters$yEC50),
            lty = 2)
    }
  }

plot.gcBootSpline <- function(x, pch=1, colData=1, colSpline=ggplot2::alpha("black", 0.1), cex=1, ...)
{
  # x an object of class gcBootSpline
  
  # /// initialize "Empty Plot" function
  empty.plot <- function(text="Empty plot",main=""){
    plot(c(0,1,0,1,0),c(0,1,1,0,0), type="l", axes=FALSE, xlab="", ylab="", lwd=1, col="gray",main=main)
    lines(c(0,0),c(0,1), type="l", lwd=1, col="gray")
    lines(c(1,1),c(1,0), type="l", lwd=1, col="gray")
    text(0.5,0.1,text, col="gray")
  }
  
  # /// check input parameters
  if (is.numeric(pch)==FALSE)   stop("Need numeric value for: pch")
  if (is.numeric(cex)==FALSE)   stop("Need numeric value for: cex")
  

  if (x$bootFlag==FALSE){
    empty.plot()
  }
  else{
    colSpline <- rep(colSpline, (x$control$nboot.gc%/%length(colSpline))+1)
    
    lambda    <- x$lambda
    mu        <- x$mu  
    A         <- x$A
    integral  <- x$integral
    
    log.x     <- x$control$log.x.gc
    log.y     <- x$control$log.y.gc
    
    global.minx <- min(min(x$boot.time,na.rm=TRUE),na.rm=TRUE)
    global.maxx <- max(max(x$boot.time,na.rm=TRUE),na.rm=TRUE)
    global.miny <- min(min(x$boot.data,na.rm=TRUE),na.rm=TRUE)
    global.maxy <- max(max(x$boot.data,na.rm=TRUE),na.rm=TRUE)
    
    # initialize plot
    if ((log.x==TRUE)&&(log.y==FALSE)){
      plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="ln(1+time)",ylab="growth y(t)")
    }
    else{
      if ((log.x==FALSE)&&(log.y==FALSE)){
        plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="time",ylab="growth y(t)")
      }
      else{
        if ((log.x==TRUE)&&(log.y==TRUE)){
          plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="ln(1+time)",ylab="Growth [Ln(y(t)/y0)]")
        }
        else{
          if ((log.x==FALSE)&&(log.y==TRUE)){
            plot(c(global.minx, global.maxx), c(global.miny, global.maxy), pch="",xlab="time",ylab="Growth [Ln(y(t)/y0)]")
          }
        }
      }
    }
    
    # /// plot data
    points(x$raw.time, x$raw.data, col=colData, pch=pch, cex=cex)
    
    # /// plot all gcFitSpline objects
    for(i in 1:x$control$nboot.gc){
      plot(x$boot.gcSpline[[i]],add=TRUE, raw=FALSE, slope=FALSE, pch=0, colSpline=colSpline[i], cex=cex)
    }
    
    # /// plot histograms of growth parameters
    dev.new()
    par(mfrow=c(2,2))
    
    if (sum(!is.na(lambda))>1){
      try(hist(lambda, col="gray",xlab="lambda", main=expression(lambda)))
    }
    else{
      empty.plot("Empty plot!")
    }
    
    if (sum(!is.na(mu))>1){ try(hist(mu , col="gray", xlab="mu", main=expression(mu))) } else { empty.plot("Empty plot!", main=expression(mu)) }
    if (sum(!is.na(A))>1){ try(hist(A, col="gray", xlab="A", main=expression(A))) } else { empty.plot("Empty plot!", main=expression(A)) }
    if (sum(!is.na(integral))>1){ try(hist(integral, col="gray", xlab="integral", main=expression(Integral))) } else { empty.plot("Empty plot!", main=expression(Integral))}
  }
  par(mfrow=c(1,1))
}

plot.gcFit <- function(x, opt="m",raw=TRUE, slope=FALSE, pch=1, colModel=1, colSpline=2, colData=1, cex=1, ...)
{
  # x an object of class gcFit
  
  # /// check input parameters
  if ( FALSE%in%(colData%in%c(colors(),0:8))) stop("colData needs to be numeric from 0:8 or a string from colors()")
  if ( FALSE%in%(colSpline%in%c(colors(),0:8))) stop("colSpline needs to be numeric from 0:8 or a string from colors()")
  if (is.logical(raw)==FALSE)       stop("Need logical value for: raw")
  if (is.logical(slope)==FALSE)     stop("Need logical value for: slope")
  if (is.numeric(pch)==FALSE)       stop("Need numeric value for: pch")
  if (is.numeric(cex)==FALSE)       stop("Need numeric value for: cex")
  if (!(opt%in%c("m","s")))     stop("Need 'm' or 's' for: opt")
  
  # /// recycle plot options
  n         <- dim(x$gcTable)[1]
  pch       <- rep(pch,       (n%/%length(pch))       +1)
  colModel  <- rep(colModel,  (n%/%length(colModel))  +1)
  colSpline <- rep(colSpline, (n%/%length(colSpline)) +1)
  colData   <- rep(colData,   (n%/%length(colData))   +1)
  
  # /// determine number of different tests
  distinct <- summary(x$gcTable[,1]) 
  k        <- length(distinct)
  
  # /// loop over all different tests
  for (j in 1:k){
    
    # /// initialize plot
    if(j>1){dev.new()}
    
    ind <- which(x$raw.data[,1]==names(distinct)[j])
    tspan <- x$raw.data
    
    data <- as.matrix((x$raw.data[ind,])[-1:-3])
    time <- as.matrix((x$raw.time[ind,]))
    
    if (x$control$log.x.gc==TRUE) time <- log(time + 1)
    if (x$control$log.y.gc==TRUE) data <- log(data + 1)
    
    tspan <- c(min(time, na.rm=TRUE), max(time, na.rm=TRUE))
    yspan <- c(min(data, na.rm=TRUE), max(data, na.rm=TRUE))
    
    scale  <- 1.025
    ts     <- ((scale-1)* diff(tspan))/2
    ys     <- ((scale-1)* diff(yspan))/2
    
    tspan0 <- tspan+c(-ts,ts)
    yspan0 <- yspan+c(-ys,ys)
    
    if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
      plot(tspan0, yspan0, xlab="time", ylab="growth y(t)", type="n")
    }
    if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
      plot(tspan0, yspan0, xlab="time", ylab="Growth [Ln(y(t)/y0)]", type="n")
    }
    if ((x$control$log.x.gc==TRUE) && (x$control$log.y.gc==FALSE)){
      plot(tspan0, yspan0, xlab="log(1+time)", ylab="growth y(t)", type="n")
    }
    if ((x$control$log.x.gc==TRUE) && (x$control$log.y.gc==TRUE)){
      plot(tspan0, yspan0, xlab="log(1+time)", ylab="Growth [Ln(y(t)/y0)]", type="n")
    }
    
    counter <- 0
    id      <- 0
    leg     <- rep("",distinct[j])
    
    # plot parametric fit
    if (opt=="m"){
      for (i in 1:n){
        if ((x$gcFittedModels[[i]]$gcID[1]==names(distinct)[j])&&( (x$gcFittedModels[[i]]$reliable==TRUE)||(is.null(x$gcFittedModels[[i]]$reliable)==TRUE) ) ){
          counter     <- counter + 1
          id[counter] <- i
          for (m in 1:length(x$gcFittedModels[[i]]$gcID)){
            leg[counter]=paste(leg[counter], as.character((x$gcFittedModels[[i]]$gcID[m])))
          }
          plot(x$gcFittedModels[[i]], add=TRUE, raw=raw, slope=slope, pch=pch[i], colData=colData[i], colModel=colModel[i], cex=cex)
        }
      }
      legend(x="topleft", pch=pch[id], leg[1:counter], col=colData[id], cex=cex, bty="n")
    }
    
    
    # plot spline fit
    if (opt=="s"){
      for (i in 1:n){
        if ((x$gcFittedSplines[[i]]$gcID[1]==names(distinct)[j])&&( (x$gcFittedSplines[[i]]$reliable==TRUE)||(is.null(x$gcFittedSplines[[i]]$reliable)==TRUE) ) ){
          counter     <- counter + 1
          id[counter] <- i
          for (m in 1:length(x$gcFittedSplines[[i]]$gcID)){
            leg[counter]=paste(leg[counter], as.character((x$gcFittedSplines[[i]]$gcID[m])))
          }
          plot(x$gcFittedSplines[[i]], add=TRUE, raw=raw, slope=slope, pch=pch[i], colData=colData[i], colSpline=colSpline[i], cex=cex)
        }
      }
      legend(x="topleft", pch=pch[id], leg[1:counter], col=colData[id], cex=cex, bty="n")
    }
    
    title(sub=names(distinct)[j])
  } # of for (j in 1:k){
} 

plot.gcFitSpline <- function(x, add=FALSE, raw=TRUE, slope=TRUE, pch=1, colData=1, colSpline="firebrick3", cex=1, lwd = 2, ...)
 {
  
  # x an object of class gcFitSpline
  
  # /// check input parameters
  if (is.logical(add)==FALSE)   stop("Need logical value for: add")
  if (is.logical(raw)==FALSE)   stop("Need logical value for: raw")
  if (is.logical(slope)==FALSE) stop("Need logical value for: slope")
  if (is.numeric(pch)==FALSE)   stop("Need numeric value for: pch")
  if (is.numeric(cex)==FALSE)   stop("Need numeric value for: cex")
  
  # /// check if a data fit is available
  if ((is.na(x$fitFlag)==TRUE)|(x$fitFlag==FALSE)){
    warning("plot.gcFitModel: no data fit available!")
  }
  else{
    if (raw==TRUE){
      if (add==TRUE){
        # /// try to plot raw data and data fit
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
          try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd) )
        }
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
          try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd) )
        }
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
          try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd ) )
        }
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
          try( points(x$raw.time, x$raw.data, sub=x$name.fit, col=colData, pch=pch, cex=cex)  )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd) )
          
        }
      }
      else{ # of if (add==TRUE)
        # /// try to plot raw data and data fit
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
          try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
          try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="time", ylab="Growth [Ln(y(t)/y0)]", col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
          try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd ) )
        }
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
          try( plot(x$raw.time, x$raw.data, sub=x$name.fit, xlab="log(1+time)", ylab="Growth [Ln(y(t)/y0)]", col=colData, pch=pch, cex=cex) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
      }
    }
    else{ # of (raw==TRUE)
      if (add==TRUE){
        # /// try to plot data fit
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
        
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
        
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd ) )
        }
        
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
      }
      else{ # of if (add==TRUE)
        # /// try to plot data fit
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==FALSE)){
          try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="growth y(t)", type="n") )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd) )
        }
        
        if ((x$control$log.x.gc==FALSE) && (x$control$log.y.gc==TRUE)){
          try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="time", ylab="Growth [Ln(y(t)/y0)]", type="n") )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd) )
        }
        
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==FALSE)){
          try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="growth y(t)", type="n" ) )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit, col=colSpline, type="l", lwd=lwd ) )
        }
        
        if ((x$control$log.x.gc==TRUE)  && (x$control$log.y.gc==TRUE)){
          try( plot(x$fit.time, x$fit.data, sub=x$name.fit, xlab="log(1+time)", ylab="Growth [Ln(y(t)/y0)]", type="n") )
          try( lines(x$fit.time, x$fit.data, sub=x$name.fit,  col=colSpline, type="l", lwd=lwd) )
        }
      }
    }
    
    # /// add tangent at maximum slope
    if (slope==TRUE){
      mu     <- as.numeric(x$parameters$mu)
      lambda <- as.numeric(x$parameters$lambda)
      
      time <- seq(lambda, max(x$"fit.time"), length=200)
      y_tangent <- x$parameters["b.tangent"][[1]]+time*mu
      try(lines(time, y_tangent, lty=2, lwd=2, col=ggplot2::alpha(colSpline, 0.7), ...))
      try(lines(c(min(x$"raw.time"[1]), lambda), rep(x$"raw.data"[1], 2), lty=2, lwd=2, col=ggplot2::alpha(colSpline, 0.7)))

    }
  }
}
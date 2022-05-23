growth.read_data <- function(data, data.format = "col", csvsep = ";")
{
  if (!is.character(data)) {
    dat <- data
  } else {
    # Read table file
    if (file.exists(data)) {
      if (stringr::str_replace_all(data, ".{1,}\\.", "") == "csv") {
        dat <-
          utils::read.csv(
            data,
            sep = csvsep,
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (stringr::str_replace_all(data, ".{1,}\\.", "") == "xls" |
                 stringr::str_replace(data, ".{1,}\\.", "") == "xlsx") {
        dat <- data.frame(read_excel(data, col_names = F))
      } else if (stringr::str_replace_all(data, ".{1,}\\.", "") == "tsv") {
        dat <-
          utils::read.csv(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else if (stringr::str_replace_all(data, ".{1,}\\.", "") == "txt") {
        dat <-
          utils::read.table(
            data,
            sep = "\t",
            header = T,
            stringsAsFactors = F,
            fill = T,
            na.strings = "",
            quote = "",
            comment.char = "",
            check.names = F
          )
      } else {
        stop(
          "No compatible file format provided.
             Supported formats are: \\.txt (tab delimited), \\.csv (delimiters can be specified with the argument \"csvsep = \", \\.tsv, \\.xls, and \\.xlsx"
        )
      }
    } else {
      stop(paste0("File \"", data, "\" does not exist."), call. = F)
    }
  }
  if(data.format == "col"){
    dat <- t(dat)
  }
  if(data.format == "col"){
    message("Sample data are stored in columns.")
  } else {
    message("Sample data are stored in rows")
  }
  if(!(any(grepl("time", unlist(dat[,1]))))){
    if(data.format == "col"){
      stop("Could not find 'time' in column 1 of dataset.")
    } else {
      stop("Could not find 'time' in row 1 of dataset.")
    }
  }
  # Create time matrix
  time.ndx <- grep("time", unlist(dat[,1]))
  if(length(time.ndx)==1){
    time <- as.numeric(unlist(dat[time.ndx[1],4:ncol(dat)]))
    t.mat <- data.matrix(data.frame(matrix(
      data = rep(time, nrow(dat)-1),
      nrow = nrow(dat)-1,
      byrow = T
    )))
  } else {
    time <- list()
    time[[1]] <- as.numeric(unlist(dat[time.ndx[1],4:ncol(dat)]))
    t.mat <- data.matrix(data.frame(matrix(
      data = rep(time[[1]], time.ndx[2]-time.ndx[1]-1),
      nrow = time.ndx[2]-time.ndx[1]-1,
      byrow = T
    )))
    for (i in 2:(length(time.ndx))){
      time[[i]] <- as.numeric(unlist(dat[time.ndx[i],4:ncol(dat)]))
      t.mat <- rbind(t.mat,
                     data.matrix(data.frame(
                       matrix(
                         data = rep(time[[i]], times = if (is.na(time.ndx[i + 1])) {
                           nrow(dat) - time.ndx[i]
                         } else {
                           time.ndx[i + 1] - time.ndx[i] - 1
                         }),
                         nrow = if (is.na(time.ndx[i + 1])) {
                           nrow(dat) - time.ndx[i]
                         } else {
                           time.ndx[i + 1] - time.ndx[i] - 1
                         },
                         byrow = T)
                     )
                     )
      )
    } # end of for (i in 2:(length(time.ndx)))
  } # end of else {}
  
  # Create data matrix
  if(length(time.ndx)==1){
    dat.mat <- data.frame(dat[(time.ndx[1]+1):nrow(dat),])
  } else {
    dat.mat <- data.frame(dat[(time.ndx[1]+1) : (time.ndx[2]-1), ])
    for (i in 2:(length(time.ndx))){
      dat.mat <- rbind(dat.mat, 
                       data.frame(dat[ if (is.na(time.ndx[i + 1])) {
                         (time.ndx[i]+1) : nrow(dat)
                       } else {
                         (time.ndx[i]+1) : (time.ndx[i+1] - 1)
                       } , ])
      )
    } # end of for (i in 2:(length(time.ndx)))
  } # end of else {}
  colnames(dat.mat)[1:3] <- c("Sample", "replicate", "concentration")
  dat.mat <- as.data.frame(unclass(dat.mat), stringsAsFactors = TRUE)

  dataset <- list("time" = t.mat,
                  "data" = dat.mat)
  invisible(dataset)
  }

growth.control <-
  function (neg.nan.act = FALSE,
            clean.bootstrap = TRUE,
            suppress.messages = FALSE,
            fit.opt = "a",
            min.density = NA,
            log.x.gc = FALSE,
            log.y.gc = TRUE,
            log.y.model = FALSE,
            interactive = TRUE,
            nboot.gc = 0,
            smooth.gc = 0.6,
            model.type = c("logistic",
                           "richards", "gompertz", "gompertz.exp"),
            have.atleast = 6,
            parameter = 34, # parameter used for creating dose response curve. # 34 is µ determined with spline fit
            smooth.dr = NULL,
            log.x.dr = FALSE,
            log.y.dr = FALSE,
            nboot.dr = 0,
            lag.method = "intersect")
{
  if ((is.character(fit.opt) == FALSE)) 
    stop("value of fit.opt must be character and of one element")
  if (is.character(model.type) == FALSE) 
    stop("value of model.type must be character")
  if ((is.logical(neg.nan.act) == FALSE) | (length(neg.nan.act) != 
                                            1)) 
    stop("value of neg.nan.act must be logical and of one element")
  if ((is.logical(clean.bootstrap) == FALSE) | (length(clean.bootstrap) != 
                                                1)) 
    stop("value of clean.bootstrap must be logical and of one element")
  if ((is.logical(suppress.messages) == FALSE) | (length(suppress.messages) != 
                                                  1)) 
    stop("value of suppress.messages must be logical and of one element")
  if ((is.logical(log.x.gc) == FALSE) | (length(log.x.gc) != 
                                         1)) 
    stop("value of log.x.gc must be logical and of one element")
  if ((is.logical(log.y.gc) == FALSE) | (length(log.y.gc) != 
                                         1)) 
    stop("value of log.y.gc must be logical and of one element")
  if ((is.logical(interactive) == FALSE) | (length(interactive) != 
                                            1)) 
    stop("value of interactive must be logical and of one element")
  if ((is.logical(log.x.dr) == FALSE) | (length(log.x.dr) != 
                                         1)) 
    stop("value of log.x.dr must be logical and of one element")
  if ((is.logical(log.y.dr) == FALSE) | (length(log.y.dr) != 
                                         1)) 
    stop("value of log.y.dr must be logical and of one element")
  if ((is.numeric(nboot.gc) == FALSE) | (length(nboot.gc) != 
                                         1) | (nboot.gc < 0)) 
    stop("value of nboot.gc must be numeric (>=0) and of one element")
  if ((is.numeric(have.atleast) == FALSE) | (length(have.atleast) != 
                                             1) | (have.atleast < 6)) 
    stop("value of have.atleast must be numeric (>=6) and of one element")
  if ((is.numeric(parameter) == FALSE) | (length(parameter) != 
                                          1)) 
    stop("value of parameter must be numeric and of one element")
  if ((is.numeric(nboot.dr) == FALSE) | (length(nboot.dr) != 
                                         1) | (nboot.dr < 0)) 
    stop("value of nboot.dr must be numeric (>=0) and of one element")
  if (((is.numeric(smooth.gc) == FALSE) && (is.null(smooth.gc) == 
                                            FALSE))) 
    stop("value of smooth.gc must be numeric or NULL")
  if (((is.numeric(smooth.dr) == FALSE) && (is.null(smooth.dr) == 
                                            FALSE))) 
    stop("value of smooth.dr must be numeric or NULL")
  grofit.control <- list(neg.nan.act = neg.nan.act, clean.bootstrap = clean.bootstrap, 
                         suppress.messages = suppress.messages, fit.opt = fit.opt, min.density = min.density,
                         log.x.gc = log.x.gc, log.y.gc = log.y.gc, log.y.model=log.y.model, interactive = interactive, 
                         nboot.gc = round(nboot.gc), smooth.gc = smooth.gc, smooth.dr = smooth.dr, 
                         have.atleast = round(have.atleast), parameter = round(parameter), 
                         log.x.dr = log.x.dr, log.y.dr = log.y.dr, nboot.dr = round(nboot.dr), 
                         model.type = model.type, lag.method = lag.method)
  class(grofit.control) <- "grofit.control"
  grofit.control
}

growth.fit <- function (time, data, t0 = 0, ec50 = FALSE, 
                        neg.nan.act = FALSE, clean.bootstrap = TRUE,
                        suppress.messages = FALSE, fit.opt = "b", min.density = NA,
                        log.x.gc = FALSE, log.y.gc = TRUE, log.y.model = FALSE,
                        interactive = TRUE, nboot.gc = 100,
                        smooth.gc= 0.6, model.type=c("logistic",
                                                      "richards","gompertz", "gompertz.exp"),
                        have.atleast = 6, parameter = 34, smooth.dr = NULL,
                        log.x.dr = FALSE, log.y.dr = FALSE, nboot.dr = 0, lag.method = "intersect") 
{
  if(!(class(test_data)=="list")){
    if (is.numeric(as.matrix(time)) == FALSE) 
      stop("Need a numeric matrix for 'time'")
    if (is.numeric(as.matrix(data[-1:-3])) == FALSE) 
      stop("Need a numeric matrix for 'data'")
    if (is.logical(ec50) == FALSE) 
      stop("Need a logical value for 'ec50'")
  } else {
    time <- data$time
    data <- data$data
  }
  if(is.numeric(min.density) && min.density != 0){
    if(!all(as.vector(data[,4]) < min.density)){
      stop(paste0("Start density values need to be below 'min.density'.\nThe minimum start value in your dataset is: ", 
                  min(as.vector(data[,4]))), call. = FALSE)
    }
  }
  control <- growth.control(neg.nan.act = neg.nan.act, clean.bootstrap = clean.bootstrap, 
                            suppress.messages = suppress.messages, fit.opt = fit.opt, min.density = min.density,
                            log.x.gc = log.x.gc, log.y.gc = log.y.gc, log.y.model = log.y.model, interactive = interactive, 
                            nboot.gc = round(nboot.gc), smooth.gc = smooth.gc, smooth.dr = smooth.dr, 
                            have.atleast = round(have.atleast), parameter = round(parameter), 
                            log.x.dr = log.x.dr, log.y.dr = log.y.dr, nboot.dr = round(nboot.dr), 
                            model.type = model.type, lag.method = lag.method)
  nboot.gc <- control$nboot.gc
  nboot.dr <- control$nboot.dr
  out.gcFit <- NA
  out.drFit <- NA
  
  # /// fit of growth curves -----------------------------------
  out.gcFit <- growth.gcFit(time, data, control, t0)
  
  # /// Estimate EC50 values -----------------------------------	
  if (ec50 == TRUE) {
    out.drFit <- growth.drFit(summary(out.gcFit), control)
    EC50.table <- out.drFit$drTable
    boot.ec <- out.drFit$boot.ec
  }
  else {
    ec50.fit <- NULL
  }
  grofit <- list(time = time, data = data, gcFit = out.gcFit, 
                 drFit = out.drFit, control = control)
  class(grofit) <- "grofit"
  grofit
}

growth.gcFit <- function(time, data, control=grofit.control(), t0)
{
  # /// check if start density values are above min.density in all samples
  if(is.numeric(control$min.density) && control$min.density != 0){
    if(!all(as.vector(data[,4]) < control$min.density)){
      stop(paste0("Start density values need to be below 'min.density'.\nThe minimum start value in your dataset is: ", 
                  min(as.vector(data[,4]))), call. = FALSE)
    }
  }
  # /// check input parameters
  if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
  
  # /// check number of datasets
  if ( (dim(time)[1])!=(dim(data)[1]) ) stop("gcFit: Different number of datasets in data and time")
  
  # /// check fitting options
  if (!all(control$fit.opt %in% c("s","m","b","a", "l"))){
    options(warn=1)
    warning("fit.opt must be set to 's', 'm', 'l', 'b', or 'a'. Changed to 'a'!")
    fit.opt="a"
    options(warn=0)
  }
  
  # /// Initialize some parameters
  out.table       <- NULL
  used.model      <- NULL
  fitpara.all     <- list()
  fitnonpara.all  <- list()
  fitlinear.all <- list()
  boot.all        <- list()
  fitted.param    <- NULL
  fitted.nonparam <- NULL
  bootstrap.param <- NULL
  reliability_tag_linear <- NA
  reliability_tag_param <- NA
  reliability_tag_nonpara <- NA
  
  
  # /// loop over all wells
  for (i in 1:dim(data)[1]){
    # /// conversion, to handle even data.frame inputs
    acttime    <- as.numeric(as.matrix(time[i,]))
    actwell <- as.numeric(as.matrix((data[i,-1:-3])))
    gcID    <- as.matrix(data[i,1:3])
    wellname <- paste(as.character(data[i,1]), as.character(data[i,2]),as.character(data[i,3]), sep="-")
    if ((control$suppress.messages==FALSE)){
      cat("\n\n")
      cat(paste("= ", as.character(i), ". growth curve =================================\n", sep=""))
      cat("----------------------------------------------------\n")
    }
    # /// Linear regression on log-transformed data
    if (("l" %in% control$fit.opt) || ("a"  %in% control$fit.opt)){
      fitlinear          <- growth.gcFitLinear(acttime, actwell, gcID = gcID, control = control, t0 = t0)
      fitlinear.all[[i]] <- fitlinear
    }
    else{
      # /// generate empty object
      fitlinear          <- list(raw.x = acttime,
                                 raw.y = actwell,
                                 gcID = gcID,
                                 FUN = NA, fit = NA,
                                 par = c(y0 = NA, y0_lm = NA, mumax = NA, lag = NA, tmax_start = NA, tmax_end = NA),
                                 ndx = NA,
                                 rsquared = NA,
                                 control = control,
                                 fitFlag = FALSE)
      class(fitlinear)   <- "gcFitLinear"
      fitlinear.all[[i]] <- fitlinear
    }
    # /// plot linear fit
    if ((control$interactive == TRUE)) {
      if (("l" %in% control$fit.opt) || ("a"  %in% control$fit.opt)) {
        if (fitlinear$fitFlag == TRUE) {
          answer_satisfied <- "n"
          reliability_tag_linear <- NA
          while ("n" %in% answer_satisfied) {
            try(plot(fitlinear, log = "y"))
            mtext(side = 3,
                  outer = F,
                  cex = 1,
                  wellname)
            answer_satisfied <- readline("Are you satisfied with the linear fit (y/n)?")
            if ("n" %in% answer_satisfied) {
              test_answer <- readline("Skip (enter 'n'), or adjust t0, h, quota, and min.density (see ?growth.gcFitLinear). >\n Enter: t0, h, quota, min.density")
              if ("n" %in% test_answer) {
                cat("\n Tagged the linear fit of this sample as unreliable !\n\n")
                reliability_tag_linear              <- FALSE
                fitlinear.all[[i]]$reliable    <- FALSE
                answer_satisfied <- "y"
              } # end if ("n" %in% test_answer)
              else {
                new_params <- unlist(strsplit(test_answer, split = ","))
                t0_new <- as.numeric(new_params[1])
                h_new <- as.numeric(new_params[2])
                quota_new <- as.numeric(new_params[3])
                min.density_new <- as.numeric(new_params[4])
                control.new <- control
                if(is.numeric(min.density_new)){
                  if(!all(as.vector(data[,4]) < control$min.density)){
                    warning(paste0("Start density values need to be below 'min.density'.\nThe minimum start value in your dataset is: ", 
                                min(as.vector(data[,4])),". 'min.density' was not adjusted."), call. = FALSE)
                  } else {
                    control.new$min.density <- min.density_new
                  }
                }
                fitlinear <- growth.gcFitLinear(acttime, actwell,
                                                         gcID = gcID, control = control.new,
                                                         t0 = t0_new, h = h_new, quota = quota_new)
                fitlinear.all[[i]] <- fitlinear
              } #end else
            } # end if ("n" %in% test_answer)
            else{
              reliability_tag_linear <- TRUE
              fitlinear.all[[i]]$reliable <- TRUE
              cat("Sample was (more ore less) o.k.\n")
            } # end else
          } # end while ("n" %in% answer_satisfied)
        } # end if (fitlinear$fitFlag == TRUE)
      } # end if (("l" %in% control$fit.opt) || ("a"  %in% control$fit.opt))
    } # end if ((control$interactive == TRUE))
    else {
      reliability_tag_linear <- TRUE
    }
    # /// Parametric fit
    if (("m" %in% control$fit.opt) || ("b" %in% control$fit.opt) || ("a"  %in% control$fit.opt)){
      fitpara          <- growth.gcFitModel(acttime, actwell, gcID, control)
      fitpara.all[[i]] <- fitpara
    }
    else{
      # /// generate empty object
      fitpara          <- list(raw.x = acttime, raw.y = actwell, gcID = gcID, fit.x = NA, fit.y = NA, parameters = list(A=NA, mu=NA, lambda=NA, integral=NA),
                               model = NA, nls = NA, reliable=NULL, fitFlag=FALSE, control = control)
      class(fitpara)   <- "gcFitModel"
      fitpara.all[[i]] <- fitpara
    }
    
    # /// Non parametric fit
    if (("s" %in% control$fit.opt) || ("b" %in% control$fit.opt) || ("a"  %in% control$fit.opt)){
      nonpara             <- growth.gcFitSpline(acttime, actwell, gcID, control, t0 = t0)
      fitnonpara.all[[i]] <- nonpara
    }
    else{
      # /// generate empty object
      nonpara             <- list(raw.x = acttime, raw.y = actwell, gcID = gcID, fit.x = NA, fit.y = NA, parameters = list(A= NA, mu=NA, lambda=NA, integral=NA),
                                  parametersLowess=list(A= NA, mu=NA, lambda=NA), spline = NA, reliable=NULL, fitFlag=FALSE, control = control)
      class(nonpara)      <- "gcFitSpline"
      fitnonpara.all[[i]] <- nonpara
    }
    
    
    # /// plotting parametric fit 
    if ((control$interactive == TRUE)) {
      if (("m" %in% control$fit.opt) || ("a"  %in% control$fit.opt) || ("b"  %in% control$fit.opt)) {
        if (fitpara$fitFlag == TRUE) {
          plot.gcFitModel(fitpara, colData=1, colModel=2, colLag = 3, cex=1.0, raw=T)
          legend(x="bottomright", legend=fitpara$model, col="red", lty=1)
          title("Parametric fit")
          mtext(side=3, outer = F, cex=1, wellname)        
          }
        # /// here a manual reliability tag is set in the interactive mode
        reliability_tag_paarm <- NA
        answer <- readline("Are you satisfied with the model fit (y/n)?")
        if ("n" %in% answer) {
          cat("\n Tagged the parametric fit of this sample as unreliable !\n\n")
          reliability_tag_param              <- FALSE
          fitpara.all[[i]]$reliable    <- FALSE
        }
        else{
          reliability_tag_param              <- TRUE
          fitpara.all[[i]]$reliable    <- TRUE
          cat("Sample was (more ore less) o.k.\n")
        }
      }
      else{
        reliability_tag_param <- FALSE
      }
      # /// plotting nonparametric fit par(mfrow=c(1,2))
      
      
      if (("s" %in% control$fit.opt) || ("a"  %in% control$fit.opt) || ("b"  %in% control$fit.opt)) {
        if (nonpara$fitFlag == TRUE) {
          answer_satisfied <- "n"
          reliability_tag_nonpara <- NA
          while ("n" %in% answer_satisfied) {
            plot.gcFitSpline(nonpara, add=FALSE, raw=TRUE, colData=1, cex=1)
            legend(x="bottomright", legend="spline fit", col=ggplot2::alpha(2, 0.5), lty=1)  
            title("Nonparametric fit")
            mtext(side=3, outer = F, cex=1, wellname)  
            answer_satisfied <- readline("Are you satisfied with the spline fit (y/n)?")
            if ("n" %in% answer_satisfied) {
                  test_answer <- readline("Skip (enter 'n'), or smooth.gc, t0, and min.density (see ?growth.control). >\n Enter: smooth.gc, t0, min.density >> ")
                  if ("n" %in% test_answer) {
                    cat("\n Tagged the linear fit of this sample as unreliable !\n\n")
                    reliability_tag_nonpara              <- FALSE
                    fitnonpara.all[[i]]$reliable    <- FALSE
                    answer_satisfied <- "y"
                  } # end if ("n" %in% test_answer)
                  else{
                  new_params <- unlist(strsplit(test_answer, split = ","))
                  if(!is.na(as.numeric(new_params[2])) && as.numeric(new_params[2]) != ""){
                    t0_new <- as.numeric(new_params[2])
                  } else {
                    t0_new <- t0
                  }
                  smooth.gc_new <- as.numeric(new_params[1])
                  control.new <- control
                  if(!is.na(smooth.gc_new) && smooth.gc_new != ""){
                    control.new$smooth.gc <- smooth.gc_new
                  }
                  min.density_new <- as.numeric(new_params[3])
                  if(is.numeric(min.density_new)){
                    if(!all(as.vector(data[,4]) < control$min.density)){
                      warning(paste0("Start density values need to be below 'min.density'.\nThe minimum start value in your dataset is: ", 
                                     min(as.vector(data[,4])),". 'min.density' was not adjusted."), call. = FALSE)
                    } else {
                      control.new$min.density <- min.density_new
                    }
                  }
                  nonpara <- growth.gcFitSpline(acttime, actwell, gcID, control.new, t0=t0_new)
                  fitnonpara.all[[i]] <- nonpara
                } #end else
              } # end if ("n" %in% answer_satisfied)
              else{
                reliability_tag_nonpara <- TRUE
                fitnonpara.all[[i]]$reliable <- TRUE
                cat("Sample was (more ore less) o.k.\n")
              } # end else
            } # end while ("n" %in% answer_satisfied)
          } # end if (nonpara$fitFlag == TRUE)
        } # end if (("s" %in% control$fit.opt) || ("a"  %in% control$fit.opt) || ("b"  %in% control$fit.opt))
    } # end of if((control$interactive == TRUE))
    else{
      reliability_tag_param <- TRUE
      reliability_tag_nonpara <- TRUE
    }
    # /// Beginn Bootstrap
    if ((("s" %in% control$fit.opt) || ("a"  %in% control$fit.opt) || ("b"  %in% control$fit.opt)) &&
        (control$nboot.gc > 0) && (reliability_tag_nonpara ==TRUE)){
      bt            <- growth.gcBootSpline(acttime, actwell, gcID, control)
      boot.all[[i]] <- bt	
    } # /// end of if (control$nboot.gc ...)
    else{
      # /// create empty gcBootSpline  object
      bt            <- list(raw.x=acttime, raw.y=actwell, gcID =gcID, boot.x=NA, boot.y=NA, boot.gcSpline=NA,
                            lambda=NA, mu=NA, A=NA, integral=NA, bootFlag=FALSE, control=control)
      class(bt)     <- "gcBootSpline"
      boot.all[[i]] <- bt	
    }
    reliability_tag <- any(reliability_tag_linear, reliability_tag_nonpara, reliability_tag_param)
    # create output table
    description     <- data.frame(TestId=data[i,1], AddId=data[i,2],concentration=data[i,3], 
                                  reliability=reliability_tag, used.model=fitpara$model, 
                                  log.x=control$log.x.gc, log.y=control$log.y.gc, nboot.gc=control$nboot.gc)
    
    fitted          <- cbind(description, summary.gcFitLinear(fitlinear), summary(fitpara), summary(nonpara), summary(bt))
    
    out.table       <- rbind(out.table, fitted)
    
  } # /// end of for (i in 1:dim(data)[1])
  
  gcFit           <- list(raw.time = time, raw.data = data, gcTable = out.table, gcFittedLinear = fitlinear.all, gcFittedModels = fitpara.all, gcFittedSplines = fitnonpara.all, gcBootSplines = boot.all, control=control)
  
  class(gcFit)    <- "gcFit"
  gcFit
  
  
  
}

growth.gcFitModel <- function(time, data, gcID ="undefined", control=grofit.control())
{
  # /// check input parameters
  if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
  if (!any(c("m","b","a") %in% control$fit.opt)) stop("Fit option is not set for a model fit. See grofit.control()")
  
  # /// conversion to handle even data.frame inputs
  time <- as.vector(as.numeric(as.matrix(time)))
  data    <- as.vector(as.numeric(as.matrix(data)))
  
  # /// check length of input data
  if (length(time)!=length(data)) stop("gcFitModel: length of time and data input vectors differ!")
  
  # /// check if there are enough data points
  if (length(data)<5){
    warning("gcFitModel: There is not enough valid data. Must have at least 5 unique values!")
    gcFitModel   <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = NA, fit.data = NA, parameters = list(A=NA, mu=NA, lambda=NA, integral=NA), model = NA, nls = NA, reliable=NULL, fitFlag=FALSE, control = control)
    class(gcFitModel) <- "gcFitModel"
    return(gcFitModel)
  }
  else{
      gcFitModel <- grofit.param(time, data, gcID, control)
  }
  return(gcFitModel)
}
  
# internal #
grofit.param <- function(time, data, gcID = "undefined", control)
{
    time.in <- time
    data.in <- data
    # Perform log transformation of data
    if (control$log.y.model == TRUE) {
      data <- log(data/data[1])
    }
    
    # /// determine which values are not valid
    bad.values <- (is.na(time))|(time<0)|(is.na(data))|(data<0)
    
    # /// remove bad values or stop program
    if (TRUE%in%bad.values){
      if (control$neg.nan.act==FALSE){
        time    <- time[!bad.values]
        data    <- data[!bad.values]
      }
      else{
        stop("Bad values in gcFitModel")
      }
    }
    
    # fitting parametric growth curves
    y.model     <- NULL
    bestAIC     <- NULL
    best        <- NULL
    used        <- NULL
    
    # starting values for parametric fitting from spline fit
    nonpara     <- growth.gcFitSpline(time.in, data.in, gcID, control)
    mu.low      <- nonpara$parametersLowess$mu
    lambda.low  <- nonpara$parametersLowess$lambda
    A.low       <- nonpara$parametersLowess$A
    
    # /// determine length of model names
    l               <- 10
    for (i in 1:length(control$model.type)) {
      l[i] <- nchar((control$model.type)[i])
    }
    lmax <- max(l)
    
    # /// loop over all parametric models requested
    for (i in 1:length(control$model.type)) {
      if (control$suppress.messages == FALSE) {
        cat(paste("--> Try to fit model", (control$model.type)[i]))
      }
      initmodel    <- paste("init", (control$model.type)[i], sep = "")
      formulamodel <-
        as.formula(paste(
          "data ~ ",
          (control$model.type)[i],
          "(time, A, mu, lambda, addpar)",
          sep = ""
        ))
      if ((exists((control$model.type)[i])) && (exists(initmodel))) {
        init.model  <-
          do.call(initmodel,
                  list(
                    y = data,
                    time = time,
                    A = A.low,
                    mu = mu.low,
                    lambda = lambda.low
                  ))
        try(y.model <-
              nls(formulamodel, start = init.model), silent = TRUE)
        if (!(TRUE %in% is.null(y.model))) {
          AIC       <- AIC(y.model)
        }
        
        if (control$suppress.messages == FALSE) {
          if (class(y.model) == "nls") {
            if (y.model$convInfo$isConv == TRUE) {
              message(paste(rep(".", lmax + 3 - l[i])), " OK")
            }
            else{
              warning(paste(
                rep(".", lmax + 3 - l[i]),
                " nls() failed to converge with stopCode ",
                as.character(y.model$convInfo$stopCode)
              ))
            }
          }
          else{
            message(paste(rep(".", lmax + 3 - l[i])),
                    " ERROR in nls(). For further information see help(gcFitModel)")
          }
        }
        if (FALSE %in% is.null(AIC)) {
          if (is.null(best) || AIC < bestAIC) {
            bestAIC <- AIC
            best    <- y.model
            used    <- (control$model.type)[i]
          }
        }
      } # of if ( (exists((control$model.type)[i])) ...
      else{
        cat((control$model.type)[i])
        cat("\n")
        cat(initmodel)
        cat("\n")
        stop("The model definition above does not exist! Spelling error?")
      }
      y.model <- NULL
    }
    
    if (control$suppress.messages == FALSE)
      cat("\n")
      cat(paste0("Best fitting model: ", sub("\\(.+", "", sub("data", "", paste(formula(best), collapse = "")))))
    # /// extract parameters from data fit
    if (is.null(best) == FALSE) {
      Abest      <- summary(best)$parameters["A", 1:2]
      mubest     <- summary(best)$parameters["mu", 1:2]
      if("addpar" %in% rownames(summary(best)$parameters)){
        fitparbest <- summary(best)$parameters["addpar", 1:2]
      }
      fitFlag    <- TRUE
      lambdabest <- summary(best)$parameters["lambda", 1:2]
      
      best.spline <- stats::smooth.spline(time, fitted.values(best), spar = 0)
      best.deriv1 <-  stats::predict(best.spline, deriv=1)
      mumax.index <- which.max(best.deriv1$y)
      
      y.max <- fitted.values(best)[mumax.index]
      t.max <- time[mumax.index]
      b.tangent <- y.max - max(best.deriv1$y) * t.max
        
      if (length(time) == length(as.numeric(fitted.values(best)))) {
        Integralbest <-
          low.integrate(time, as.numeric(fitted.values(best)))
      }
      else{
        Integralbest <- NA
      }
    }
    else{
      warning("gcFitModel: Unable to fit this curve parametrically!")
      Abest        <- c(NA, NA)
      mubest       <- c(NA, NA)
      lambdabest   <- c(NA, NA)
      Integralbest <- NA
      fitFlag      <- FALSE
      b.tangent <- NA
    }
      
    gcFitModel <-
      list(
        raw.time = time,
        raw.data = data,
        gcID = gcID,
        fit.time = time,
        fit.data = as.numeric(fitted.values(best)),
        parameters = list(
          A = if (control$log.y.model == TRUE) {
            # Correct ln(N/N0) transformation for max density value
            data[1] * exp(Abest)
          } else {
            Abest
          } ,
          mu = mubest,
          lambda = lambdabest,
          b.tangent = b.tangent,

          fitpar = if(exists("fitparbest")){
            fitparbest
          } else {
            NULL
          },
          integral = Integralbest
        ),
        model = used,
        nls = best,
        reliable = NULL,
        fitFlag = fitFlag,
        control = control
      )
    
    class(gcFitModel) <- "gcFitModel"
    
    invisible(gcFitModel)
  }

growth.gcFitSpline <- function (time, data, gcID = "undefined", control = grofit.control(), t0 = 0) 
{
  if(!is.null(t0) && !is.na(t0) && t0 != ""){
    t0 <- as.numeric(t0)
  } else {
    t0 <- 0
  }
  
  if (is(control) != "grofit.control") 
    stop("control must be of class grofit.control!")
  if (!any(control$fit.opt %in% c("s", "b", "a"))) 
    stop("Fit option is not set for a spline fit. See grofit.control()")
  time <- as.vector(as.numeric(as.matrix(time)))
  data <- as.vector(as.numeric(as.matrix(data)))

  if (length(time) != length(data)) 
    stop("gcFitSpline: length of input vectors differ!")
  bad.values <- (is.na(time)) | (is.na(data)) | 
    (!is.numeric(time)) | (!is.numeric(data))
  if (TRUE %in% bad.values) {
    if (control$neg.nan.act == FALSE) {
      time <- time[!bad.values]
      data <- data[!bad.values]
    }
    else {
      stop("Bad values in gcFitSpline")
    }
  }
  if (length(data) < 5) {
    cat("gcFitSpline: There is not enough valid data. Must have at least 5!")
    gcFitSpline <- list(raw.time = time, raw.data = data, 
                        gcID = gcID, fit.time = NA, fit.data = NA, parameters = list(A = NA, 
                                                                                     mu = NA, lambda = NA, integral = NA), parametersLowess = list(A = NA, 
                                                                                                                                                   mu = NA, lambda = NA), spline = NA, reliable = NULL, 
                        fitFlag = FALSE, control = control)
    class(gcFitSpline) <- "gcFitSpline"
    return(gcFitSpline)
  }
  else {
    if (control$log.x.gc == TRUE) {
      bad.values <- (time < 0)
      if (TRUE %in% bad.values) {
        if (control$neg.nan.act == FALSE) {
          time <- time[!bad.values]
          data <- data[!bad.values]
        }
        else {
          stop("Bad values in gcFitSpline")
        }
      }
      time <- log(1 + time)
    }
    if (control$log.y.gc == TRUE) {
      data.log <- log(data/data[1])
      bad.values <- (data.log < 0)
      if (TRUE %in% bad.values) {
        if (control$neg.nan.act == FALSE) {
          time <- time[!bad.values]
          data.log <- data.log[!bad.values]
        }
        else {
          stop("Bad values in gcFitSpline")
        }
      }
    }
    time.raw <- time
    data.raw <- if (control$log.y.gc == TRUE) {
      data.log
    } else {
      data
    }
    if(!is.na(control$min.density) && control$min.density != 0){
      if (control$log.y.gc == TRUE) {
        min.density <- log(control$min.density / data[1])
        time <- time[max(which.min(abs(time-t0)), which.min(abs(data.log-min.density))):length(time)]
        data.log <- data.log[max(which.min(abs(time.raw-t0)), which.min(abs(data.log-min.density))) : length(data.log)]
      } else {
        min.density <- control$min.density
        time <- time[max(which.min(abs(time-t0)), which.min(abs(data-min.density))):length(data)]
        data <- data[max(which.min(abs(time.raw-t0)), which.min(abs(data-min.density))) : length(data)]
      }
    } else {
      if (control$log.y.gc == TRUE) {
        data.log <- data.log[which.min(abs(time-t0)):length(data.log)]
      } 
      time <- time[which.min(abs(time-t0)):length(time)]
      data <- data[which.min(abs(time.raw-t0)):length(data)]
    }
    halftime <- (min(time) + max(time))/2
    try(y.spl <- smooth.spline(time, y = if(control$log.y.gc == TRUE){
      data.log
    } else {
      data
    }, spar = control$smooth.gc))
    if (is.null(y.spl) == TRUE) {
      warning("Spline could not be fitted to data!")
      if (is.null(control$smooth.gc) == TRUE) {
        cat("This might be caused by usage of smoothing parameter NULL\n")
        fit.nonpara <- list(raw.x = time, raw.y = data, 
                            fit.x = NA, fit.y = NA, parameters = list(A = NA, 
                                                                      mu = NA, lambda = NA, integral = NA), spline = NA, 
                            parametersLowess = list(A = NA, mu = NA, lambda = NA), 
                            spline = NA, reliable = NULL, fitFlag = FALSE, 
                            control = control)
        class(gcFitSpline) <- "gcFitSpline"
        return(gcFitSpline)
      }
    }
    dydt.spl <- predict(y.spl, time, deriv = 1)
    mumax.index <- which.max(dydt.spl$y) # index of data point with maximum growth rate
    t.max <- dydt.spl$x[mumax.index] # time of maximum growth rate
    dydt.max <- max(dydt.spl$y) # maximum value of first derivative of spline fit (i.e., greatest slope in growth curve spline fit)
    y.max <- y.spl$y[mumax.index] # cell density at time of max growth rate
    mu.spl <- dydt.max # maximum growth rate
    b.spl <- y.max - dydt.max * t.max # the y-intercept of the tangent at µmax
    lambda.spl <- -b.spl/mu.spl # lag time
    integral <- low.integrate(y.spl$x, y.spl$y)
    low <- lowess(time, y = if(control$log.y.gc == TRUE){
      data.log
    } else {
      data
    }, f = 0.25)
    y.low <- low$y
    x.low <- low$x
    dydt.low <- diff(y.low)/diff(time)
    mu.low <- max(dydt.low)
    mumax.index.low <- which.max(dydt.low)
    t.max.low <- x.low[mumax.index.low]
    y.max.low <- y.low[mumax.index.low]
    b.low <- y.max.low - mu.low * t.max.low
    lambda.low <- (-1) * b.low/mu.low
  }
  gcFitSpline <-
    list(
      raw.time = time.raw,
      raw.data = data.raw,
      gcID = gcID,
      fit.time = y.spl$x,
      fit.data = y.spl$y,
      parameters = list(
        A = if (control$log.y.gc == TRUE) {
          # Correct ln(N/N0) transformation for max density value
          data[1] * exp(max(y.spl$y))
        } else {
          max(y.spl$y)
        },
        mu = mu.spl,
        lambda = lambda.spl,
        b.tangent = b.spl,
        integral = integral),
      parametersLowess = list(
        A = max(y.low),
        mu = mu.low,
        lambda = lambda.low
      ),
      spline = y.spl,
      spline.deriv1 = dydt.spl,
      reliable = NULL,
      fitFlag = TRUE,
      control = control
    )
  class(gcFitSpline) <- "gcFitSpline"
  gcFitSpline
}

#' Fit Exponential Growth Model with a Heuristic Linear Method
#'
#' Determine maximum growth rates from the log-linear part of a growth curve using
#' a heuristic approach similar to the ``growth rates made easy''-method of
#' Hall et al. (2013).
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Fit linear regressions to all subsets of \code{h} consecutive data
#'     points (sliding window). If for example \eqn{h=5}, fit a linear regression to points
#'     1 \dots 5, 2 \dots 6, 3 \dots 7 and so on. The method seeks the highest
#'     rate of exponential growth, so the dependent variable is of course
#'     log-transformed.
#'   \item Find the subset with the highest slope \eqn{b_{max}} and
#'     include also the data points of adjacent subsets that have a slope of
#'     at least \eqn{quota \cdot b_{max}},
#'     e.g. all data sets that have at least 95\% of the maximum slope.
#'   \item Fit a new linear model to the extended data window identified in step 2.
#' }
#'
#' @param time vector of independent variable.
#' @param y vector of dependent variable (concentration of organisms).
#' @param h width of the window (number of data).
#' @param quota part of window fits considered for the overall linear fit
#'   (relative to max. growth rate)
#'
#' @return object with parameters of the fit. The lag time is currently estimated
#' as the intersection between the fit and the horizontal line with \eqn{y=y_0},
#' where \code{y0} is the first value of the dependent variable. The intersection
#' of the fit with the abscissa is indicated as \code{y0_lm} (lm for linear model).
#' These identifieres and their assumptions may change in future versions.
#'
#' @references Hall, BG., Acar, H, Nandipati, A and Barlow, M (2014) Growth Rates Made Easy.
#' Mol. Biol. Evol. 31: 232-38, \doi{10.1093/molbev/mst187}
#'
#' @family fitting functions
#'
#' @examples
#' data(bactgrowth)
#'
#' splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#' dat <- splitted.data[[1]]
#'
#' plot(value ~ time, data=dat)
#' fit <- fit_easylinear(dat$time, dat$value)
#'
#' plot(fit)
#' plot(fit, log="y")
#' plot(fit, which="diagnostics")
#'
#' fitx <- fit_easylinear(dat$time, dat$value, h=8, quota=0.95)
#'
#' plot(fit, log="y")
#' lines(fitx, pch="+", col="blue")
#'
#' plot(fit)
#' lines(fitx, pch="+", col="blue")
#' 
#' @export 
#'
growth.gcFitLinear <- function(time, data, gcID = "undefined", t0 = 0, h = NULL, quota = 0.95, control) 
  {
  bad.values <- ((is.na(time))|(is.na(data)))
  data.in <- data[!bad.values]
  time.in <- time[!bad.values]
  if(!is.null(t0) && !is.na(t0) && t0 != ""){
    t0 <- as.numeric(t0)
  } else {
    t0 <- 0
  }

  if(!is.null(h) && !is.na(h) && h != ""){
    h <- as.numeric(h)
  } else {
    # extract period of growth (from defined t0)
    t.growth <- time[which.min(abs(time-t0)):which.max(data)]
    # determine number of data points in period until maximum density
    n.spl <- length(t.growth)
    # Calculate h via log-transformation of the number of data points
    h <- round((log(n.spl+5, base=2.4))/0.92)
  }
  
  if(!is.null(quota) && !is.na(quota) && quota != ""){
    if(quota > 1){
      quota <- as.numeric(quota)/100
    } else {
      quota <- as.numeric(quota)
    }
  } else {
    quota <- 0.95
  }
  data.log <- log(data/data[1])
  if(!is.na(control$min.density) && control$min.density != 0){
    min.density <- log(control$min.density / data[1])
  }
  bad.values <- ((is.na(data.log))|(is.infinite(data.log))|data.log<0|(is.na(time))|(is.na(data)))
  # /// remove bad values or stop program
  if (TRUE%in%bad.values){
    if (control$neg.nan.act==FALSE){
      time    <- time[!bad.values]
      data.log    <- data.log[!bad.values]
      data <- data[!bad.values]
    }
    else{
      stop("Bad values in gcFitModel")
    }
  }
  if (any(duplicated(time))) stop("time variable must not contain duplicated values")
  
  obs <- data.frame(time, data)
  obs$ylog <- data.log
  
  ## number of values
  N <- nrow(obs)
  
  ## repeat for all windows and save results in 'ret'
  ret <- matrix(0, nrow = N - h, ncol = 6)
  for(i in 1:(N - h)) {
    ret[i, ] <- c(i, with(obs, suppressWarnings(lm_parms(lm_window(time, ylog, i0 = i, h = h)))))
  }
  # add time values as column in ret
  ret <- data.frame(ret, time = time[ret[,1]], data = obs$ylog[ret[,1]])
  # duplicate ret for further tuning of fit
  if(exists("min.density")){
    ret.check <- ret[max(which.min(abs(time-t0)), which.min(abs(obs$ylog-min.density))) : nrow(ret),] # consider only slopes from defined t0
  } else{
    ret.check <- ret[which.min(abs(time-t0)):nrow(ret),] # consider only slopes from defined t0
  }
  
  ## Determine indices of windows with high growth rate
  success <- FALSE
  while (!success){
    index.max <- which.max(ret.check[, 3]) 
    if(ret.check[index.max,5] >= 0.95 & ret.check[index.max,6] <= 0.1){ # prerequisites for suitable µmax candidate: R2 >= 0.95 and RSD <= 0.05
      slope.max <- ret.check[index.max,3]
      success <- TRUE
    } else {
      ret.check <- ret.check[-index.max,]
    }
  }
  index.max.ret <- ret.check[which(ret.check[,3]==slope.max),1] # index of maximum slope in original fit table
  slope.quota <- quota * slope.max
  if(exists("min.density")){
    candidates <- which(ret[, 3] >= slope.quota & # indices of slopes greater than slope.quota
                          ret[, 5] >= 0.90 & # R2 criterion for candidates
                          ret[, 6] <= 0.1 & # RSD criterion for candidates
                          ret[, 7] >= t0 & # consider only slopes after defined t0
                          ret[, 8] >= min.density # cosider only slopes at densities higher than "min.density"
                        ) 
  } else{
    candidates <- which(ret[, 3] >= slope.quota & # indices of slopes greater than slope.quota
                          ret[, 5] >= 0.90 & # R2 criterion for candidates
                          ret[, 6] <= 0.1 & # RSD criterion for candidates
                          ret[, 7] >= t0) # consider only slopes after defined t0
    }
  
  
  if(length(candidates) > 0) {
    tp <- seq(min(candidates), max(candidates) + h-1)
    m <- lm_window(obs$time, obs$ylog, min(tp), length(tp)) # linear model
    p  <- c(lm_parms(m), n=length(tp)) # # slope equation parameters (linear model)
  } else {
    p <- c(a=0, b=0, se=0, r2=0, cv=0, n=0)
    m = NULL
  }

  if(length(candidates) > 0) {
    ## get time window of exponential fit
    tmax_start <- obs$time[tp[1]]
    tmax_end <- obs$time[tp[length(tp)]]
    
    y0_lm    <- unname(coef(m)[1]) # y-intercept of tangent
    y0_data  <- obs$ylog[1] # y0 in dataset
    mumax <- unname(coef(m)[2])
    
    ## estimate lag phase
    lambda <- (y0_data - y0_lm) / mumax
    
    # correct y0 values for Ln(y(t)/y0) 
    y0_lm <- obs$data[1] * exp(y0_lm)
    y0_data <- obs$data[1]
    
    # get indices of time points used in linear fit
    ndx <- seq(min(match(ret[candidates, "time"], time.in)), 
               max(match(ret[candidates, "time"], time.in)) + h-1)
  } else {
    tmax_start <- NULL
    tmax_end <- NULL
    
    y0_lm    <- NULL
    y0_data  <- NULL
    mumax <- NULL
    
    ## estimate lag phase
    lambda <- NULL
    
    # correct y0 values for Ln(y(t)/y0) 
    y0_lm <- NULL
    y0_data <- NULL
    
    # get indices of time points used in linear fit
    ndx <- NULL
  }
  mu.se <- p[3] # standard error of slope
  
  

  gcFitLinear <- list(
    raw.time = time.in,
    raw.data = data.in,
    filt.time = obs$time,
    filt.data = obs$data,
    log.data = obs$ylog,
    gcID = gcID,
    FUN = grow_exponential,
    fit = m,
    par = c(
      y0 = y0_data,
      y0_lm = y0_lm,
      mumax = mumax,
      mu.se = mu.se,
      lag = lambda,
      tmax_start = tmax_start,
      tmax_end = tmax_end
    ),
    ndx = ndx,
    rsquared = p["r2"],
    control = control,
    fitFlag = TRUE
  )
  
  class(gcFitLinear) <- "gcFitLinear"
  
  invisible(gcFitLinear)
}

growth.gcBootSpline <- function (time, data, gcID = "undefined", control = grofit.control()) 
{
  if (is(control) != "grofit.control") 
    stop("control must be of class grofit.control!")
  if (control$nboot.gc == 0) 
    stop("Number of bootstrap samples is zero! See grofit.control()")
  time <- as.vector(as.numeric(as.matrix(time)))
  data <- as.vector(as.numeric(as.matrix(data)))
  if (length(time) != length(data)) 
    stop("gcBootSpline: length of input vectors differ!")
  bad.values <- (is.na(time)) | (is.na(data)) | is.infinite(data) |
    (!is.numeric(time)) | (!is.numeric(data))
  if (TRUE %in% bad.values) {
    if (control$neg.nan.act == FALSE) {
      time <- time[!bad.values]
      data <- data[!bad.values]
    }
    else stop("Bad values in gcBootSpline")
  }
  if (length(data) < 6) {
    warning("gcBootSpline: There is not enough valid data. Must have at least 6 unique values!")
    gcBootSpline <- list(raw.time = time, raw.data = data, 
                         gcID = gcID, boot.time = NA, boot.data = NA, boot.gcSpline = NA, 
                         lambda = NA, mu = NA, A = NA, integral = NA, bootFlag = FALSE, 
                         control = control)
    class(gcBootSpline) <- "gcBootSpline"
    return(gcBootSpline)
  }
  A <- NA
  mu <- NA
  lambda <- NA
  integral <- NA
  boot.y <- array(NA, c(control$nboot.gc, length(time)))
  boot.x <- array(NA, c(control$nboot.gc, length(time)))
  nonpara <- list()
  control.change <- control
  control.change$fit.opt <- "s"
  if (control$nboot.gc > 0) {
    for (j in 1:control$nboot.gc) {
      choose <- sort(sample(1:length(time), length(time), replace = TRUE))
      while (length(unique(choose)) < 5) {
        choose <- sort(sample(1:length(time), length(time), 
                         replace = TRUE))
      }
      time.cur <- time[choose]
      data.cur <- data[choose]
      if(IQR(time.cur) > 0){
        nonpara[[j]] <- growth.gcFitSpline(time.cur, data.cur, gcID, 
                                    control.change)
        boot.y[j, 1:length(nonpara[[j]]$fit.data)] <- nonpara[[j]]$fit.data
        boot.x[j, 1:length(nonpara[[j]]$fit.time)] <- nonpara[[j]]$fit.time
        lambda[j] <- nonpara[[j]]$parameters$lambda
        mu[j] <- nonpara[[j]]$parameters$mu
        A[j] <- nonpara[[j]]$parameters$A
        integral[j] <- nonpara[[j]]$parameters$integral
      }
    }
    lambda[which(!is.finite(lambda))] <- NA
    mu[which(!is.finite(lambda))] <- NA
    A[which(!is.finite(lambda))] <- NA
    integral[which(!is.finite(lambda))] <- NA
    if (control$clean.bootstrap == TRUE) {
      lambda[which(lambda < 0)] <- NA
      mu[which(mu < 0)] <- NA
      A[which(A < 0)] <- NA
      integral[which(integral < 0)] <- NA
    }
  }
  if (control$log.x.gc == TRUE) {
    bad.values <- (time < 0)
    if (TRUE %in% bad.values) {
        time <- time[!bad.values]
        data <- data[!bad.values]
    }
    time.log <- log(1 + time)
  }
  if (control$log.y.gc == TRUE) {
    data.log <- log(data/data[1])
    bad.values <- (data.log < 0)
    if (TRUE %in% bad.values) {
        time <- time[!bad.values]
        data.log <- data.log[!bad.values]
    }
  }
  gcBootSpline <- list(
    raw.time = if (control$log.x.gc == TRUE) {
      time.log
    } else {
      time
    },
    time,
    raw.data = if (control$log.y.gc == TRUE) {
      data.log
    } else {
      data
    },
    gcID = gcID,
    boot.time = boot.x,
    boot.data = boot.y,
    boot.gcSpline = nonpara,
    lambda = lambda,
    mu = mu,
    A = A,
    integral = integral,
    bootFlag = TRUE,
    control = control
  )
  class(gcBootSpline) <- "gcBootSpline"
  gcBootSpline
}

growth.drFit <- function (gcFitData, control = grofit.control()) 
{
  if (is(control) != "grofit.control") 
    stop("control must be of class grofit.control!")
  EC50.table <- NULL
  all.EC50 <- NA
  table.tests <- table((gcFitData[, 1])[which((gcFitData[, 
                                                         4] == TRUE) & (is.na(gcFitData[, control$parameter]) == 
                                                                          FALSE))])
  distinct <- names(table.tests)
  EC50 <- list()
  EC50.boot <- list()
  validdata <- cbind(as.character(distinct), table.tests)
  colnames(validdata) <- c("TestID", "Number")
  rownames(validdata) <- rep("     ", times = dim(validdata)[1])
  if (control$suppress.messages == FALSE) {
    cat("\n")
    cat("=== EC 50 Estimation ==============================\n")
    cat("---------------------------------------------------\n")
    cat("--> Checking data ...\n")
    cat(paste("--> Number of distinct tests found:", as.character(length(distinct))), 
        "\n")
    cat("--> Valid datasets per test: \n")
    print(validdata, quote = FALSE)
  }
  if (TRUE %in% (table.tests < control$have.atleast)) {
    cat(paste("Warning: following tests have not enough ( <", 
              as.character(control$have.atleast - 1), ") datasets:\n"))
    cat(distinct[(table.tests < control$have.atleast)])
    cat("These tests will not be regarded\n")
    distinct <- distinct[table.tests >= control$have.atleast]
  }
  if ((length(distinct)) == 0) {
    cat(paste("There are no tests having enough ( >", as.character(control$have.atleast - 
                                                                     1), ") datasets!\n"))
  }
  else {
    for (i in 1:length(distinct)) {
      conc <- (gcFitData[, 3])[which(gcFitData[, 1] == 
                                       distinct[i])]
      test <- (gcFitData[, control$parameter])[which(gcFitData[, 
                                                               1] == distinct[i])]
      drID <- distinct[i]
      EC50[[i]] <- growth.drFitSpline(conc, test, drID, control)
      if (control$nboot.dr > 0) {
        EC50.boot[[i]] <- drBootSpline(conc, test, drID, 
                                       control)
      }
      else {
        EC50.boot[[i]] <- list(raw.x = conc, raw.y = test, 
                               drID = drID, boot.x = NA, boot.y = NA, boot.drSpline = NA, 
                               ec50.boot = NA, bootFlag = FALSE, control = control)
        class(EC50.boot[[i]]) <- "drBootSpline"
      }
      description <- data.frame(Test = distinct[i], log.x = control$log.x.dr, 
                                log.y = control$log.y.dr, Samples = control$nboot.dr)
      out.row <- cbind(description, summary(EC50[[i]]), 
                       summary(EC50.boot[[i]]))
      EC50.table <- rbind(EC50.table, out.row)
    }
  }
  drFit <- list(raw.data = gcFitData, drTable = EC50.table, 
                drBootSplines = EC50.boot, drFittedSplines = EC50, control = control)
  class(drFit) <- "drFit"
  drFit
}

growth.drFitSpline <- function (conc, test, drID = "undefined", control = grofit.control()) 
{
  if (is(control) != "grofit.control") 
    stop("control must be of class grofit.control!")
  test <- as.vector(as.numeric(as.matrix(test)))
  conc <- as.vector(as.numeric(as.matrix(conc)))
  if (is.vector(conc) == FALSE || is.vector(test) == FALSE) 
    stop("drFitSpline: dose or response data must be a vector !")
  if (control$neg.nan.act == FALSE) {
    missings <- is.na(conc) | is.na(test) | !is.numeric(conc) | 
      !is.numeric(test)
    conc <- conc[!missings]
    test <- test[!missings]
    negs <- (conc < 0) | (test < 0)
    conc <- conc[!negs]
    test <- test[!negs]
  }
  else {
    if (sum(is.na(conc) | is.na(test))) 
      stop("drFitSpline: NA values encountered. Program terminated")
    if ((sum((conc < 0)) > 0) | (sum((test < 0)) > 0)) 
      stop("drFitSpline: Negative values encountered. Program terminated")
    if ((FALSE %in% is.numeric(conc)) || (FALSE %in% is.numeric(test))) 
      stop("drFitSpline: Non numeric values encountered. Program terminated")
  }
  if (length(test) < 6) {
    warning("drFitSpline: There is not enough valid data. Must have at least 6 unique values!")
    drFitSpline <- list(raw.conc = conc, raw.test = test, 
                        drID = drID, fit.conc = NA, fit.test = NA, spline = NA, 
                        parameters = list(EC50 = NA, yEC50 = NA, EC50.orig = NA, 
                                          yEC50.orig = NA), fitFlag = FALSE, reliable = NULL, 
                        control = control)
    class(drFitSpline) <- "drFitSpline"
    return(drFitSpline)
  }
  if (length(test) < control$have.atleast) {
    warning("drFitSpline: number of valid data points is below the number specified in 'have.atleast'. See grofit.control().")
    drFitSpline <- list(raw.conc = conc, raw.test = test, 
                        drID = drID, fit.conc = NA, fit.test = NA, spline = NA, 
                        parameters = list(EC50 = NA, yEC50 = NA, EC50.orig = NA, 
                                          yEC50.orig = NA), fitFlag = FALSE, reliable = NULL, 
                        control = control)
    class(drFitSpline) <- "drFitSpline"
    return(drFitSpline)
  }
  if (control$log.x.dr == TRUE) 
    conc.log <- log(conc + 1)
  if (control$log.y.dr == TRUE) 
    test.log <- log(test + 1)
  if (control$log.x.dr == TRUE) {
    conc.fit <- log(conc + 1)
  }
  else {
    conc.fit <- conc
  }
  if (control$log.y.dr == TRUE) {
    test.fit <- log(test + 1)
  }
  else {
    test.fit <- test
  }
  spltest <- NULL
  fitFlag <- TRUE
  try(spltest <- smooth.spline(conc.fit, test.fit, spar = control$smooth.dr))
  if (is.null(spltest) == TRUE) {
    cat("Spline could not be fitted to data!\n")
    fitFlag <- FALSE
    if (is.null(control$smooth.dr) == TRUE) {
      cat("This might be caused by usage of smoothing parameter NULL\n")
    }
    stop("Error in drFitSpline")
  }
  conc.min <- min(conc.fit)
  conc.max <- max(conc.fit)
  c.pred <- seq(conc.min, conc.max, length.out = 1000)
  ytest <- predict(spltest, c.pred)
  yEC.test <- (max(ytest$y) - min(ytest$y))/2 + min(ytest$y)
  last.test <- max(ytest$y)
  kec.test <- 1
  for (k in 1:(length(c.pred) - 1)) {
    d1 <- (ytest$y[k] - yEC.test)
    d2 <- (ytest$y[k + 1] - yEC.test)
    if (((d1 <= 0) && (d2 >= 0)) | ((d1 >= 0) && (d2 <= 0))) {
      kec.test <- k
      break
    }
  }
  EC.test <- c.pred[kec.test]
  if (control$suppress.messages == FALSE) {
    cat("\n\n=== Dose response curve estimation ================\n")
    cat("--- EC 50 -----------------------------------------\n")
    cat(paste("-->", as.character(drID)))
    cat("\n")
    cat(paste(c("xEC50", "yEC50"), c(EC.test, yEC.test)))
  }
  if ((control$log.x.dr == TRUE) && (control$log.y.dr == FALSE)) {
    if (control$suppress.messages == FALSE) {
      cat("\n--> Original scale \n")
      cat(paste(c("xEC50", "yEC50"), c(exp(EC.test) - 1, 
                                       yEC.test)))
    }
    EC.orig <- c(exp(EC.test) - 1, yEC.test)
  }
  else {
    if ((control$log.x.dr == FALSE) && (control$log.y.dr == 
                                        TRUE)) {
      if (control$suppress.messages == FALSE) {
        cat("\n--> Original scale \n")
        cat(paste(c("xEC50", "yEC50"), c(EC.test, exp(yEC.test) - 
                                           1)))
      }
      EC.orig <- c(EC.test, exp(yEC.test) - 1)
    }
    else {
      if ((control$log.x.dr == TRUE) && (control$log.y.dr == 
                                         TRUE)) {
        if (control$suppress.messages == FALSE) {
          cat("\n--> Original scale \n")
          cat(paste(c("xEC50", "yEC50"), c(exp(EC.test) - 
                                             1, exp(yEC.test) - 1)))
        }
        EC.orig <- c(exp(EC.test) - 1, exp(yEC.test) - 
                       1)
      }
      else {
        if ((control$log.x.dr == FALSE) && (control$log.y.dr == 
                                            FALSE)) {
          EC.orig <- c(EC.test, yEC.test)
        }
      }
    }
  }
  if (control$suppress.messages == FALSE) {
    cat("\n\n\n")
  }
  drFitSpline <- list(raw.conc = conc, raw.test = test, drID = drID, 
                      fit.conc = ytest$x, fit.test = ytest$y, spline = spltest, 
                      parameters = list(EC50 = EC.test[1], yEC50 = yEC.test, 
                                        EC50.orig = EC.orig[1], yEC50.orig = EC.orig[2]), 
                      fitFlag = fitFlag, reliable = NULL, control = control)
  class(drFitSpline) <- "drFitSpline"
  drFitSpline
}

lm_parms <- function (m) 
{
  sm <- summary(m)
  a <- sm$coefficients[1, 1]
  b <- sm$coefficients[2, 1]
  b.se <- sm$coefficients[2, 2]
  r2 <- sm$r.squared
  c(a = a, b = b, b.se = b.se, r2 = r2, b.rsd = b.se/b)
}

lm_window <- function (x, y, i0, h = 5) 
{
  x <- x[i0 - 1 + (1:h)]
  y <- y[i0 - 1 + (1:h)]
  m <- lm(y ~ x)
  return(m)
}

grow_exponential <- function (time, parms) 
{
  if (is.null(names(parms))) {
    y0 <- parms[1]
    mumax <- parms[2]
  }
  else {
    y0 <- parms["y0"]
    mumax <- parms["mumax"]
  }
  y <- y0 * exp(mumax * time)
  return(as.matrix(data.frame(time = time, y = y)))
}

low.integrate <- function (x, y) 
{
  if (is.vector(x) == FALSE || is.vector(y) == FALSE) 
    stop("low.integrate: two vectors x and y are needed !")
  if (length(x) != length(y)) 
    stop("low.integrate: x and y have to be of same length !")
  y.spl <- NULL
  try(y.spl <- smooth.spline(x, y))
  if (is.null(y.spl) == TRUE) {
    warning("Spline could not be fitted to data!")
    stop("Error in low.integrate")
  }
  f <- function(t) {
    p <- predict(y.spl, t)
    f <- p$y
  }
  low.integrate <- integrate(f, min(x), max(x))$value
}

logistic <- function (time, A, mu, lambda, addpar = NULL) 
{
  A <- A[1]
  mu <- mu[1]
  lambda <- lambda[1]
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  y <- A/(1 + exp(4 * mu * (lambda - time)/A + 2))
  logistic <- y
}

initlogistic <- function (time, y, A, mu, lambda) 
{
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(y) == FALSE) 
    stop("Need numeric vector for: y")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  A <- max(y)
  mu <- mu[1]
  lambda <- lambda[1]
  initlogistic <- list(A = A, mu = mu, lambda = lambda, addpar = NULL)
}

initrichards <- function (time, y, A, mu, lambda) 
{
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(y) == FALSE) 
    stop("Need numeric vector for: y")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  nu <- 0.1
  A <- max(y)
  mu <- mu[1]
  lambda <- lambda[1]
  initrichards <- list(A = A, mu = mu, lambda = lambda, addpar = nu)
}

richards <- function (time, A, mu, lambda, addpar) 
{
  A <- A[1]
  mu <- mu[1]
  lambda <- lambda[1]
  nu <- addpar[1]
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  if (is.numeric(nu) == FALSE) 
    stop("Need numeric vector for: addpar[1]")
  y <- A * (1 + nu * exp(1 + nu) * exp(mu * (1 + nu)^(1 + 1/nu) * 
                                         (lambda - time)/A))^(-1/nu)
  richards <- y
}

gompertz <- function (time, A, mu, lambda, addpar = NULL) 
{
  A <- A[1]
  mu <- mu[1]
  lambda <- lambda[1]
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  e <- exp(1)
  y <- A * exp(-exp(mu * e * (lambda - time)/A + 1))
  gompertz <- y
}

initgompertz <- function (time, y, A, mu, lambda) 
{
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(y) == FALSE) 
    stop("Need numeric vector for: y")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  A <- max(y)
  mu <- mu[1]
  lambda <- lambda[1]
  initgompertz <- list(A = A, mu = mu, lambda = lambda, addpar = NULL)
}

initgompertz.exp <- function (time, y, A, mu, lambda) 
{
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(y) == FALSE) 
    stop("Need numeric vector for: y")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  alfa <- 0.1
  tshift <- max(time)/10
  A <- max(y)
  mu <- mu[1]
  lambda <- lambda[1]
  initgompertz.exp <- list(A = A, mu = mu, lambda = lambda, 
                           addpar = c(alfa, tshift))
}

gompertz.exp <- function (time, A, mu, lambda, addpar) 
{
  A <- A[1]
  mu <- mu[1]
  lambda <- lambda[1]
  alpha <- addpar[1]
  tshift <- addpar[2]
  if (is.numeric(time) == FALSE) 
    stop("Need numeric vector for: time")
  if (is.numeric(mu) == FALSE) 
    stop("Need numeric vector for: mu")
  if (is.numeric(lambda) == FALSE) 
    stop("Need numeric vector for: lambda")
  if (is.numeric(A) == FALSE) 
    stop("Need numeric vector for: A")
  if (is.numeric(alpha) == FALSE) 
    stop("Need numeric vector for: addpar[1]")
  if (is.numeric(tshift) == FALSE) 
    stop("Need numeric vector for: addpar[2]")
  e <- exp(1)
  y <- A * exp(-exp(mu * e * (lambda - time)/A + 1)) + A * 
    exp(alpha * (time - tshift))
  gompertz.exp <- y
}
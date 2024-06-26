#######################################################
# Print - S3 Method
#######################################################


#' Basic print of objects of class 'sdeTMB'
#' @returns A huge amount of information
#' @export
print.sdeTMB = function(object,...) {
  
  obj = object$print()
  #
  return(invisible(obj))
  
}

#' Basic print of objects of class 'sdeTMB'
#' @returns A huge amount of information
#' @export
print.sdeTMB.fit = function(fit) {
  
  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  
  return(invisible(mat))
  
}


#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'sdeTMB'
#' @returns A huge amount of information
#' @export
summary.sdeTMB = function(object,correlation=FALSE) {
  
  obj = object$summary(correlation)
  #
  return(invisible(obj))
  
}

#' @returns summary of fit object from \code{obj$estimate()}
#' @export
summary.sdeTMB.fit = function(fit,correlation=FALSE) {
  
  if (!is.logical(correlation)) {
    stop("correlation must be logical")
  }
  
  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)
  if (correlation){
    cat("\nCorrelation Matrix\n")
    S = diag(1/private$fit$sd.fixed)
    cor = S %*% private$fit$cov.fixed %*% S
    rownames(cor) = names(fit$par.fixed)
    colnames(cor) = names(fit$par.fixed)
    cor <- format(round(cor, 2), nsmall = 2, digits = 6)
    cor[!lower.tri(cor,diag=T)] <- ""
    print(cor, drop = FALSE, quote = FALSE)
    # cor[!lower.tri(cor)] <- ""
    # print(cor[-1, -(dim(cor)[1]), drop = FALSE], quote = FALSE)
  }
  
  return(invisible(list(parameters=mat)))
}

#######################################################
# GGPLOT2 FUNCTIONS FOR USE IN PLOTS
#######################################################

getggplot2theme = function() {
  mytheme =
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text("Avenir Next Condensed",size=12),
      legend.text = ggplot2::element_text(size=12),
      axis.text = ggplot2::element_text(size=12),
      strip.text = ggplot2::element_text(face="bold",size=12),
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      legend.box = "vertical",
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust=0.5)
    )
  return(mytheme)
}

getggplot2colors = function(n) {
  hues = seq(15, 375, length = n + 1)
  ggcolors = hcl(h = hues, l = 65, c = 100)[1:n]
  return(ggcolors)
}

#######################################################
# Plot - S3 Method
#######################################################

#' Basic summary of objects of class 'sdeTMB.pred' from predict
#' @returns A huge amount of information
#' @export
plot.sdeTMB.pred = function(pred.data,
                          n.ahead=0,
                          state.name=NULL,
                          ...) {
  
  # check n.ahead
  if(!any(n.ahead %in% unique(pred$n.ahead))){
    stop("n.ahead not found in the prediction data frame")
  }
  
  # set state name to plot
  if(is.null(state.name)){
    state.name = colnames(pred)[6]
  }
  
  # filter and plot
  bool = pred$n.ahead==n.ahead
  x = pred[bool,"t_{k+i}"]
  y = pred[bool,state.name]
  plot(x=x, y=y, type="l",...)
  
  
  return(invisible(NULL))
}


#' Basic summary of objects of class 'sdeTMB'
#' @param plot.obs a vector to indicate which observations should be plotted for. If multiple
#' are chosen a list of plots for each observation is returned.
#' @param pacf logical to indicate whether or not the partial autocorrelations should be returned.
#' The default is FALSE in which case a histogram is returned instead.
#' @param extended logical. if TRUE additional information is printed
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @returns A huge amount of information
#' @export
plot.sdeTMB = function(object,
                     plot.obs=1,
                     pacf=FALSE,
                     extended=FALSE,
                     ggtheme=getggplot2theme()) {
  
  object$plot(plot.obs=plot.obs, pacf=pacf, extended=extended, ggtheme=ggtheme)
  
  return(invisible(NULL))
}

#' Basic summary of objects of class 'sdeTMB'
#' @param plot.obs a vector to indicate which observations should be plotted for. If multiple
#' are chosen a list of plots for each observation is returned.
#' @param pacf logical to indicate whether or not the partial autocorrelations should be returned.
#' The default is FALSE in which case a histogram is returned instead.
#' @param extended logical. if TRUE additional information is printed
#' @param ggtheme ggplot2 theme to use for creating the ggplot.
#' @returns A list of plots
#' @export
#' @importFrom stats frequency fft spec.taper
plot.sdeTMB.fit = function(fit,
                         plot.obs=1,
                         pacf=FALSE,
                         extended=FALSE,
                         ggtheme=getggplot2theme()) {
  
  # get private fields from object
  private = fit$.__object__$.__enclos_env__$private
  
  if (!is.logical(pacf)) {
    stop("pacf must be logical")
  }
  if (!is.logical(extended)) {
    stop("extended must be logical")
  }
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("ggtheme must be a ggplot2 theme")
  }
  if(!is.numeric(plot.obs)){
    stop("plot.obs must be a numeric (or integer) value")
  }
  if(plot.obs > private$number.of.states){
    message("Can't plot state ", plot.obs, " because there's only ", private$number.of.states, " state(s). Setting plot.obs = ",private$number.of.states)
    plot.obs = private$number.of.states
  }
  
  if(private$method == "laplace"){
    message("'plot' is not available for method 'laplace' yet.")
    return(NULL)
  }
  
  # retrieve user default parameter settings
  
  # use ggplot to plot
  mycolor = getggplot2colors(2)[2]
  # if (use.ggplot) {
  plots = list()
  t = fit$residuals$normalized[["t"]]
  for (i in seq_along(private$obs.eqs.trans)) {
    e = fit$residuals$normalized[[private$obs.names[i]]]
    id = !is.na(e)
    e = e[id]
    t = t[id]
    nam = private$obs.names[i]
    
    # time vs residuals
    plot.res =
      ggplot2::ggplot(data=data.frame(t,e)) +
      ggplot2::geom_point(ggplot2::aes(x=t,y=e),color=mycolor) +
      ggtheme +
      ggplot2::labs(
        title = paste("Time Series of Residuals: "),
        y = "",
        x = "Time"
      )
    # quantile plots
    plot.qq =
      ggplot2::ggplot(data=data.frame(e)) +
      ggplot2::stat_qq(ggplot2::aes(sample=e),color=mycolor) +
      ggplot2::stat_qq_line(ggplot2::aes(sample=e),lty="dashed") +
      ggtheme +
      ggplot2::labs(
        title = paste("Normal Q-Q Plot: "),
        y = "Sample Quantiles",
        x = "Theoretical Quantiles"
      )
    # histogram
    plot.hist =
      ggplot2::ggplot(data=data.frame(e)) +
      ggplot2::geom_histogram(ggplot2::aes(x=e,y=after_stat(density)),bins=100,color="black",fill=mycolor) +
      ggtheme +
      ggplot2::labs(
        title = paste("Histogram: "),
        y = "",
        x = ""
      )
    # acf
    myacf = stats::acf(e,na.action=na.pass,plot=FALSE)
    plot.acf =
      ggplot2::ggplot(data=data.frame(lag=myacf$lag[-1],acf=myacf$acf[-1])) +
      ggplot2::geom_errorbar(ggplot2::aes(x=lag,ymax=acf,ymin=0),width=0,color=mycolor) +
      ggplot2::geom_hline(yintercept=0) +
      ggplot2::geom_hline(yintercept=c(-2/sqrt(myacf$n.used),2/sqrt(myacf$n.used)),color="blue",lty="dashed") +
      ggtheme +
      coord_cartesian(xlim=c(0,ceiling(max(myacf$lag)/10)*10)) +
      ggplot2::labs(
        title = paste("Auto-Correlation: "),
        y = "",
        x = "Lag"
      )
    # pacf
    mypacf = stats::pacf(e,na.action=na.pass,plot=FALSE)
    plot.pacf = 
      ggplot2::ggplot(data=data.frame(lag=mypacf$lag[-1], pacf=mypacf$acf[-1])) +
      ggplot2::geom_errorbar(ggplot2::aes(x=lag,ymax=pacf,ymin=0),width=0,color=mycolor) +
      ggplot2::geom_hline(yintercept=0) +
      ggplot2::geom_hline(yintercept=c(-2/sqrt(mypacf$n.used),2/sqrt(mypacf$n.used)),color="blue",lty="dashed") +
      ggtheme +
      coord_cartesian(xlim=c(0,ceiling(max(mypacf$lag)/10)*10)) +
      ggplot2::labs(
        title = paste("Partial Auto-Correlation: "),
        y = "",
        x = "Lag"
      )
    # cpgram
    plot.cpgram =
      ggfortify::ggcpgram(e,colour=mycolor) +
      ggtheme +
      ggplot2::labs(
        title = paste("Cumulative Periodogram: "),
        y = "",
        x = "Lag"
      )
    #
    templist = list(
      plot.res,
      plot.hist,
      plot.qq,
      plot.acf,
      plot.cpgram
    )
    # store plot
    plots[[i]] = (plot.qq + plot.hist) / (plot.acf + plot.cpgram) +
      patchwork::plot_annotation(title=paste("Residual Analysis for ", nam),theme=ggtheme + ggplot2::theme(text=ggplot2::element_text("Avenir Next Condensed",size=18,face="bold")),)
    if(pacf){
      plots[[i]] = (plot.qq + plot.pacf) / (plot.acf + plot.cpgram) +
        patchwork::plot_annotation(title=paste("Residuals Analysis for ", nam),theme=ggtheme + ggplot2::theme(text=ggplot2::element_text("Avenir Next Condensed",size=18,face="bold")),)
    }
    # }
    # return plot list, and print one of the plots
    print(plots[[plot.obs]])
    return(invisible(plots))
  }
  
  # return
  return(invisible(NULL))
}

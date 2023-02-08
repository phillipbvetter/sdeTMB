#######################################################
# Print - S3 Method
#######################################################


#' Basic print of objects of class 'sdemTMB'
#' @returns A huge amount of information
#' @export
print.sdemTMB = function(object,...) {

  obj = object$print()
  #
  return(invisible(obj))

}

#' Basic print of objects of class 'sdemTMB'
#' @returns A huge amount of information
#' @export
print.sdemTMB.fit = function(fit) {

  mat = cbind(fit$par.fixed,fit$sd.fixed,fit$tvalue,fit$Pr.tvalue)
  colnames(mat) = c("Estimate","Std. Error","t value","Pr(>|t|)")
  cat("Coefficent Matrix \n")
  stats::printCoefmat(mat)

  return(invisible(mat))

}


#######################################################
# Summary - S3 Method
#######################################################

#' Basic summary of objects of class 'sdemTMB'
#' @returns A huge amount of information
#' @export
summary.sdemTMB = function(object,correlation=FALSE) {

  obj = object$summary(correlation)
  #
  return(invisible(obj))

}

#' @returns summary of fit object from \code{obj$estimate()}
#' @export
summary.sdemTMB.fit = function(fit,correlation=FALSE) {

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

#' Basic summary of objects of class 'sdemTMB'
#' @param extended logical. if TRUE additional information is printed
#' @returns A huge amount of information
#' @export
plot.sdemTMB = function(object,
                     plot.obs=1,
                     extended=FALSE,
                     use.ggplot=FALSE,
                     ggtheme=getggplot2theme()
) {

  object$plot(plot.obs=plot.obs, use.ggplot=use.ggplot, extended=extended, ggtheme=ggtheme)
  return(invisible(NULL))
}

#' Basic summary of objects of class 'sdemTMB'
#' @param extended logical. if TRUE additional information is printed
#' @returns A huge amount of information
#' @export
#' @importFrom stats frequency fft spec.taper
plot.sdemTMB.fit = function(fit,
                         plot.obs=1,
                         use.ggplot=FALSE,
                         extended=FALSE,
                         ggtheme=getggplot2theme()
){

  if (!is.logical(use.ggplot)) {
    stop("use.ggplot must be logical")
  }
  if (!is.logical(extended)) {
    stop("extended must be logical")
  }
  if (!(inherits(ggtheme,"theme") & inherits(ggtheme,"gg"))) {
    stop("ggtheme must be a ggplot2 theme")
  }

  private = fit$.__object__$.__enclos_env__$private

  if(private$method == "tmb"){
    message("'plot' is not available for method 'tmb' yet.")
    return(NULL)
  }

  mycolor = getggplot2colors(2)[2]
  if (use.ggplot) {
    plots = list()
    t = fit$residuals$normalized[["t"]]
    for (i in 1:private$m) {
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
        ggplot2::geom_histogram(ggplot2::aes(x=e,y=..density..),bins=100,color="black",fill=mycolor) +
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
        ggplot2::labs(
          title = paste("Auto-Correlation: "),
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
        patchwork::plot_annotation(title=paste("Residuals for ", nam),theme=ggtheme + ggplot2::theme(text=ggplot2::element_text("Avenir Next Condensed",size=18,face="bold")),)

    }
    # return plot list, and print one of the plots
    print(plots[[plot.obs]])
    return(invisible(plots))
  }

  # if we aren't using ggplot we can't return the plots, we just print one plot
  graphics::par(mfrow=c(2,2))
  e = fit$residuals$normalized[[plot.obs+1]]
  nam = private$obs.names[plot.obs]
  # qqnorm and line
  stats::qqnorm(e, main=paste("Normal Q-Q Plot:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1,-3:3)
  graphics::axis(2)
  graphics::grid()
  stats::qqline(e)
  # histogram
  graphics::hist(e,main=paste("Histogram:"),breaks=50)
  graphics::grid()
  graphics::axis(1)
  graphics::axis(2)
  # acf
  max.lag = round(10*log10(length(stats::na.omit(e))))
  plot(acf(e,plot=F,na.action=na.pass)[1:max.lag],main=paste("Auto-Correlation:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1,seq(0,max.lag,round(max.lag/10)))
  graphics::axis(2)
  # cpgram
  x <- e[!is.na(e)]
  x <- stats::spec.taper(scale(x, TRUE, FALSE), p = 0.1)
  y <- Mod(stats::fft(x))^2/length(x)
  y[1L] <- 0
  n <- length(x)
  x <- (0:(n/2)) * stats::frequency(e)/n
  if (length(x)%%2 == 0) {
    n <- length(x) - 1
    y <- y[1L:n]
    x <- x[1L:n]
  } else {
    y <- y[seq_along(x)]
  }
  xm <- stats::frequency(e)/2
  mp <- length(x) - 1
  crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
  plot(x, cumsum(y)/sum(y), type = "s", xlim = c(0, xm), ylim = c(0,1),
       xaxs = "i", yaxs = "i", xlab = "frequency", ylab = "",
       main=paste("Cumulative Periodogram:"), axes=FALSE)
  graphics::lines(c(0, xm * (1 - crit)), c(crit, 1), col = "blue", lty = 2)
  graphics::lines(c(xm * crit, xm), c(0, 1 - crit), col = "blue", lty = 2)
  # stats::cpgram(e,main=paste("Cumulative Periodogram:"),axes=FALSE)
  graphics::grid()
  graphics::axis(1)
  graphics::axis(2)
  # overall title
  graphics::title(main=paste("Residuals for ", nam),outer=TRUE,line=-1.25,cex.main=1.5)

  # return
  return(invisible(NULL))
}

#######################################################
# Predict - S3 Method
#######################################################

#' @title Predict state estimates on \code{sdemTMB.fit} object.
#' @param fit An \code{sdemTMB.fit} object from calling \code{estimate} method on \code{sdemTMB} object.
#' @param data data.frame used to predict and update the model on which \code{fit} is based.
#' @param n.step Integer. Number of prediction steps before updating. Default: 1
#' @param use.simulation Boolean. If TRUE simulates the stochastic model and returns quantiles. if FALSE uses
#' and returns the solution values for the first and second order moment differential equations. Default: FALSE
#' @param n.sim Integer. Number of simulations performed if \code{use.simulation} is TRUE. Default: 100
#' @param quantiles The returned state distribution quantiles if \code{use.simulation} is TRUE.
#' @returns Prediction from model
#' @export
predict.sdemTMB.fit = function(fit,
                            data=NULL,
                            x0=NULL,
                            p0=NULL,
                            n.step.ahead=1,
                            ode.timestep=NULL,
                            return.state.dispersion=FALSE,
                            covariance=FALSE,
                            give.only.n.step.ahead=FALSE,
                            use.simulation=FALSE,
                            sim.timestep=NULL,
                            returned.sim.dt=NULL,
                            n.sim=100,
                            quantiles=c(0.05,0.5,0.95),
                            return.full.simulations=FALSE
){



  ############################################################
  ############################################################
  ##################### INITIALIZATION #######################
  ############################################################
  ############################################################

  # extract private fields from object
  private = fit$.__object__$.__enclos_env__$private

  # if data is null use fit data
  if(is.null(data)){
    data = private$data
  }
  # if x0 is NULL use fit data
  if(is.null(x0)){
    x0 = private$initial.state$mean
  }
  # if p0 is NULL use fit data
  if(is.null(p0)){
    p0 = private$initial.state$cov
  }

  if(n.step.ahead < 1){
    stop("You must choose a positive integer for n.step.ahead")
  }
  if(n.sim < 1){
    stop("You must choose a positive integer for n.sim")
  }
  if(!is.logical(use.simulation)){
    stop("The provided value for use.simulation is not logical")
  }
  if(!is.logical(return.state.dispersion)){
    stop("The provided value for return.state.dispersion is not logical")
  }
  if(!is.logical(covariance)){
    stop("The provided value for covariance is not logical")
  }
  if(!is.logical(give.only.n.step.ahead)){
    stop("The provided value for give.only.n.step.ahead is not logical")
  }
  if(!is.logical(return.full.simulations)){
    stop("The provided value for get.full.simulations is not logical")
  }
  if(any(quantiles>1 | quantiles<0)){
    stop("Please choose quantiles betwen 0 and 1")
  }
  input.bool = private$input.names %in% names(data)
  if(!all(input.bool)){
    stop("The data frame does not contain any entry for the input(s): ",private$input.names[!input.bool])
  }

  # ODE time-step
  if (!is.null(ode.timestep)) {
    # must be numeric
    if (!is.numeric(dt)) {
      stop("ode.timestep should be a numeric value.")
    }
    # must have length 1
    if (length(dt)!=1) {
      stop("ode.timestep must have length 1.")
    }
  }
  if (is.null(ode.timestep)) {
    ode.timestep = min(diff(data$t))
  }
  if (ode.timestep - min(diff(data$t)) > 1e-12) {
    ode.timestep = min(diff(data$t))
    message(sprintf("'ode.timestep' is too large - reducing to min(diff(data$t)) = %s...", format(ode.timestep,digits=10,scientific=T)))
  }
  # correct time-step to give integer steps - if the step-size is N.01 or larger we use N+1 steps instead of N
  ode.N = floor(diff(data$t)/ode.timestep)
  bool = diff(data$t)/ode.timestep - ode.N > 1e-2
  ode.N[bool] = ode.N[bool] + 1
  ode.dt = diff(data$t)/ode.N

  # Simulation time-step
  if (is.null(sim.timestep)) {
    sim.timestep = min(diff(data$t))/10
  }
  if (sim.timestep - min(diff(data$t)) > 1e-12) {
    sim.timestep = min(diff(data$t))/10
    message(sprintf("'sim.timestep' is too large - reducing to min(diff(data$t))/10 = %s...", format(sim.timestep,digits=10,scientific=T)))
  }
  # correct time-step to give integer steps - if the step-size is N.01 or larger we use N+1 steps instead of N
  sim.N = floor(diff(data$t)/sim.timestep)
  bool = diff(data$t)/sim.timestep - sim.N > 1e-2
  sim.N[bool] = sim.N[bool] + 1
  sim.dt = diff(data$t)/sim.N

  # last index for predictions
  last.possible.index = nrow(data) - n.step.ahead
  if(last.possible.index < 1){
    stop("There are only ", nrow(data), " time points in the provided data, so it is not possible to choose n.step.ahead = ",n.step.ahead)
  }


  ############################################################
  ##################### KALMAN FILTERS #######################
  ############################################################

  if (private$method == "ekf" | private$method == "ukf") {

    ##################### helper functions for construction ####################
    fun.wrapper = function(name,states,inputs,return.name.type,body,return.name){
      paste(c(name,"=function(states,inputs){",states,inputs,return.name.type,";",body,"return(",return.name,")}"),collapse="")
    }
    fun.states = paste(private$state.names,"=","states[",1:private$n,"]",";",sep="",collapse="")
    fun.inputs = paste(private$input.names,"=","inputs[",1:length(private$input.names),"]",";",sep="",collapse="")
    # vectorized helper
    fun.states.vec = paste(private$state.names,"=","states[,",1:private$n,"]",";",sep="",collapse="")


    ##################### construct drift function f ####################
    f = c()
    for(i in 1:private$n){
      f[i] = sprintf("f[%s]=%s",i,paste(deparse(private$diff.terms[[i]]$dt),collapse=""))
    }
    f. = paste(c(paste(f,collapse=";"),";"),collapse="")
    txt = fun.wrapper("f",fun.states,fun.inputs,"f=c()",f.,"f")
    eval(parse(text=txt))

    # vectorized verrsion of f to evaluate many simulations with one call
    txt = fun.wrapper("fvector",fun.states.vec,fun.inputs,"f=c()",f.,"f")
    eval(parse(text=txt))
    #
    f_ = c()
    for(i in 1:private$n){
      f_[i] = sprintf("f=c(f,rep_len(%s,length.out=length(%s)))",paste(deparse(private$diff.terms[[i]]$dt),collapse=""),private$state.names[1])
    }
    f. = paste(c(paste(f_,collapse=";"),";"),collapse="")
    txt = fun.wrapper("fvector",fun.states.vec,fun.inputs,"f=c()",f.,"f")
    eval(parse(text=txt))

    ############### construct jacobian of drift function dfdx ############
    dfdx = matrix(0,nrow=private$n,ncol=private$n)
    for(i in 1:private$n){
      terms = lapply(private$state.names,function(x) Deriv::Deriv(private$diff.terms[[i]]$dt,x=x,cache.exp=FALSE))
      for(j in 1:private$n){
        dfdx[i,j] = sprintf("dfdx[%s,%s]=%s",i,j,paste(deparse(terms[[j]]),collapse=""))
      }
    }
    dfdx. = paste(c(paste(dfdx,collapse=";"),";"),collapse="")
    txt = fun.wrapper("dfdx",fun.states,fun.inputs,sprintf("dfdx=diag(%s)",private$n),dfdx.,"dfdx")
    eval(parse(text=txt))

    ##################### construct diffusion function g ####################
    g = matrix(0,nrow=private$n,ncol=private$ng)
    for(i in 1:private$n){
      for(j in 1:private$ng){
        g[i,j] = sprintf("g[%s,%s]=%s",i,j,paste(deparse(private$diff.terms[[i]][[j+1]]),collapse=""))
      }
    }
    g. = paste(c(paste(g,collapse=";"),";"),collapse="")
    txt = fun.wrapper("g",fun.states,fun.inputs,sprintf("g=matrix(0,nrow=%s,ncol=%s)",private$n,private$ng),g.,"g")
    eval(parse(text=txt))
    # vectorized version of g to evaluate many simulations with one call
    g_ = c()
    iii = 1
    for(i in 1:(private$n)){
      for(j in 1:private$ng){
        g_[iii] = sprintf("g=c(g,rep_len(%s,length.out=length(%s)))",paste(deparse(private$diff.terms[[i]][[j+1]]),collapse=""),private$state.names[1])
        iii = iii + 1
      }
    }
    g. = paste(c(paste(g_,collapse=";"),";"),collapse="")
    txt = fun.wrapper("gvector",fun.states.vec,fun.inputs,"g=c()",g.,"g")
    eval(parse(text=txt))

    # function to easily compute sum G %* dw
    g_dot_dw = function(gvec) {
      temp = c()
      start_id = 1 + (1:private$n-1)*private$ng*n.sim
      end_id = 1:private$n * private$ng * n.sim
      for(i in 1:private$n) {
        temp = cbind(temp, rowSums(matrix(gvec[start_id[i]:end_id[i]],ncol=private$ng)))
      }
      return(temp)
    }

    ##################### construct observation function h ####################
    h = c()
    for(i in 1:private$m){
      h[i] = sprintf("h[%s]=%s",i,paste(deparse(private$obs.eqs.trans[[i]]$rhs),collapse=""))
    }
    h. = paste(c(paste(h,collapse=";"),";"),collapse="")
    txt = fun.wrapper("h",fun.states,fun.inputs,"h=c()",h.,"h")
    eval(parse(text=txt))

    ##################### construct jacobian of observation function dhdx ####################
    dhdx = matrix(NA,nrow=private$m,ncol=private$n)
    for(i in 1:private$m){
      terms = lapply(private$state.names, function(x) Deriv::Deriv(private$obs.eqs.trans[[i]]$rhs,x=x,cache.exp=FALSE))
      for(j in 1:private$n){
        dhdx[i,j] = sprintf("dhdx[%s,%s]=%s",i,j,paste(deparse(terms[[j]]),collapse=""))
      }
    }
    dhdx. = paste(c(paste(dhdx,collapse=";"),";"),collapse="")
    txt = fun.wrapper("dhdx",fun.states,fun.inputs,sprintf("dhdx=matrix(0,nrow=%s,ncol=%s)",private$m,private$n),dhdx.,"dhdx")
    eval(parse(text=txt))

    ##################### observation variance function ####################
    obs_f = c()
    for(i in 1:private$m){
      obs_f[i] = sprintf("obs_f[%s]=%s",i,paste(deparse(private$obs.var.trans[[i]]$rhs),collapse=""))
    }
    obs_f. = paste(c(paste(obs_f,collapse=";"),";"),collapse="")
    txt = fun.wrapper("obs_f",fun.states,fun.inputs,"obs_f=c()",obs_f.,"obs_f")
    eval(parse(text=txt))

    ##################### assign parameter values #######################
    for(parname in names(private$fixed.pars)){
      # print(parname)
      assign(parname,private$parameters[[parname]]$init)
    }
    for(parname in names(fit$par.fixed)){
      # print(parname)
      assign(parname,fit$par.fixed[parname])
    }

    ##################### MAIN FOR-LOOP ########################
    ##################### MAIN FOR-LOOP ########################
    ##################### MAIN FOR-LOOP ########################


    ##########################################
    ####### Branch 1: Simulations ############
    ##########################################
    if(use.simulation){

      ##### CREATE DATA FRAME TO HOLD RESULTS #####

      # create data.frame to hold state predictions and quantiles
      df.out = data.frame(matrix(nrow=last.possible.index*(n.step.ahead+1), ncol=5+private$n+private$n*length(quantiles)))
      quantile.names = paste(rep(private$state.names,each=length(quantiles)),paste0("quantile(",quantiles,")"),sep="-")
      names(df.out) = c("k","k+i","n.step.ahead","t_{k}","t_{k+i}",private$state.names,quantile.names)
      ran = 0:(last.possible.index-1)
      df.out["k"] = rep(ran,each=n.step.ahead+1)
      df.out["k+i"] = df.out["k"] + rep(0:n.step.ahead,last.possible.index)
      df.out["n.step.ahead"] = rep(0:n.step.ahead,last.possible.index)
      df.out["t_{k}"] = rep(data$t[ran+1],each=n.step.ahead+1)
      df.out["t_{k+i}"] = data$t[df.out[,"k"]+1+rep(0:n.step.ahead,last.possible.index)]

      # create list of data.frames to hold simulations
      df.sim = stats::setNames(
        lapply(1:last.possible.index,function(x) stats::setNames(replicate(private$n,data.frame(matrix(nrow=n.sim),fix.empty.names=FALSE),simplify=FALSE),private$state.names)),
        paste0("t=",data$t[1:last.possible.index]))

      ##### KALMAN FILTER INITIALIZATION #####

      datmat = as.matrix(data)
      I.n = diag(private$n)
      V.obs = diag(private$n)
      dw_number = private$ng * n.sim
      len.g = private$n * private$ng * n.sim

      ##### LOOP OVER DATA TIME POINTS #####

      for (i in 1:last.possible.index){
        # store the posterior estimate (initial or previous loop)
        index_i_i = (i-1)*(n.step.ahead+1) + 1
        df.out[index_i_i, private$state.names] = x0

        # compute the n.step.ahead ahead priors (predictions): x{i+k|i}, k = 1..n.step.ahead
        for(k in 1:n.step.ahead){
          inputs = datmat[i+k-1,private$input.names]

          # solve moment odes using ode.N time-steps with ode.dt step size.
          for(j in 1:ode.N[i+k-1]){
            A. = dfdx(x0,inputs)
            G. = g(x0,inputs)
            x0 = x0 + f(x0,inputs) * ode.dt[i+k-1]
            p0 = p0 + (A. %*% p0 + p0 %*% t(A.) + G. %*% t(G.)) * ode.dt[i+k-1]
          }

          # store n.step.ahead ahead state x_{i+k|i} and covariance p_{i+k|i}
          index_i_k = (i-1) * (n.step.ahead+1) + k + 1
          df.out[index_i_k,private$state.names] = x0
        }

        # recover state x{i+1|i} and covariance p{i+1|i} one-step prior
        index_i_1 = (i-1) * (n.step.ahead+1) + 2
        x0 = unlist(df.out[index_i_1, private$state.names])

        # store x0 in first columns
        .x0 = matrix(x0,ncol=private$n,nrow=n.sim,byrow=T)
        for(ii in seq_len(private$n)){
          df.sim[[i]][[ii]][,1] = .x0[,ii]
        }
        counter = 1
        # simulate between t_{i} and t_{i+n.step.ahead}
        for(k in 1:n.step.ahead){
          inputs = datmat[i+k-1,private$input.names]
          sqrt_dt = sqrt(sim.dt[i+k-1])

          # simulate between t_{i+k-1} and t_{i+k}
          for(j in 1:sim.N[i+k-1]){
            dw = rep(stats::rnorm(dw_number),times=private$n)
            .f = matrix(fvector(.x0,inputs),ncol=private$n) * sim.dt[i+k-1]
            .g_times_dw = g_dot_dw(gvector(.x0,inputs) * dw * sqrt_dt)
            .x0 = .x0 + .f + .g_times_dw
            # store results
            counter = counter + 1
            for(ii in seq_len(private$n)){
              df.sim[[i]][[ii]][,counter] = .x0[,ii]
            }
          }
        }

        # recover state x{i+1|i} and covariance p{i+1|i} one-step prior
        index_i_1 = (i-1) * (n.step.ahead+1) + 2
        x0 = unlist(df.out[index_i_1, private$state.names])

        # data:update to update from prior to posterior
        inputs = datmat[i,private$input.names]
        diag(V.obs) = obs_f(x0,inputs)
        y = datmat[i+1,private$obs.names]
        nas = !is.na(y)
        s = sum(nas)
        if (s > 0){
          ynew = y[nas]
          E. = diag(private$n)[nas,,drop=FALSE]
          e. = ynew - E. %*% h(x0,inputs)
          C. = E. %*% dhdx(x0,inputs)
          V. = E. %*% V.obs %*% t(E.)
          R. = C. %*% p0 %*% t(C.) + V.
          Ri. = solve(R.)
          K. = p0 %*% t(C.) %*% Ri.
          # update states
          x0  = x0 + K. %*% e.
          p0  = (I.n - K. %*% C.) %*% p0 %*% t(I.n - K. %*% C.) + K. %*% V. %*% t(K.)
        }
      }

      # return
      return(list(pred=df.out,simulations=df.sim))
    }

    ##########################################
    ### Branch 2: kalman filter predictions ##
    ##########################################
    if(!use.simulation){

      ##### CREATE OUTPUT DATA FRAME #####
      df.out = data.frame(matrix(nrow=last.possible.index*(n.step.ahead+1), ncol=5+private$n+private$n^2))
      disp_names = sprintf(rep("cor[%s,%s]",private$n^2),rep(private$state.names,each=private$n),rep(private$state.names,private$n))
      disp_names[seq.int(1,private$n^2,by=private$n+1)] = sprintf(rep("var[%s]",private$n),private$state.names)
      var_bool = !stringr::str_detect(disp_names,"cor")
      if(covariance){
        disp_names = sprintf(rep("cov[%s,%s]",private$n^2),rep(private$state.names,each=private$n),rep(private$state.names,private$n))
        disp_names[seq.int(1,private$n^2,by=private$n+1)] = sprintf(rep("var[%s]",private$n),private$state.names)
      }
      names(df.out) = c("k","k+i","n.step.ahead","t_{k}","t_{k+i}",private$state.names,disp_names)
      ran = 0:(last.possible.index-1)
      df.out["k"] = rep(ran,each=n.step.ahead+1)
      df.out["k+i"] = df.out["k"] + rep(0:n.step.ahead,last.possible.index)
      df.out["n.step.ahead"] = rep(0:n.step.ahead,last.possible.index)
      df.out["t_{k}"] = rep(data$t[ran+1],each=n.step.ahead+1)
      df.out["t_{k+i}"] = data$t[df.out[,"k"]+1+rep(0:n.step.ahead,last.possible.index)]

      ##### KALMAN FILTER INITIALIZATION #####

      datmat = as.matrix(data)
      I.n = diag(private$n)
      V.obs = diag(private$n)

      ##### LOOP OVER DATA TIME POINTS #####

      for (i in 1:last.possible.index){
        # store the posterior estimate (initial or previous loop)
        index_i_i = (i-1)*(n.step.ahead+1) + 1
        df.out[index_i_i, private$state.names] = x0
        df.out[index_i_i,disp_names] = as.vector(stats::cov2cor(p0))
        df.out[index_i_i,disp_names[var_bool]] = as.vector(diag(p0))
        if(covariance){
          df.out[index_i_i, disp_names] = as.vector(p0)
        }

        # compute the n.step.ahead ahead priors (predictions): x{i+k|i}, k = 1..n.step.ahead
        for(k in 1:n.step.ahead){
          inputs = datmat[i+k-1,private$input.names]

          # solve moment odes using ode.N time-steps with ode.dt step size.
          for(j in 1:ode.N[i+k-1]){
            A. = dfdx(x0,inputs)
            G. = g(x0,inputs)
            x0 = x0 + f(x0,inputs) * ode.dt[i+k-1]
            p0 = p0 + (A. %*% p0 + p0 %*% t(A.) + G. %*% t(G.)) * ode.dt[i+k-1]
          }

          # store n.step.ahead ahead state x_{i+k|i} and covariance p_{i+k|i}
          index_i_k = (i-1) * (n.step.ahead+1) + k + 1
          df.out[index_i_k,private$state.names] = x0
          df.out[index_i_k, disp_names] = as.vector(p0)
        }

        # recover state x{i+1|i} and covariance p{i+1|i} one-step prior
        index_i_1 = (i-1) * (n.step.ahead+1) + 2
        x0 = unlist(df.out[index_i_1, private$state.names])
        p0 = matrix(unlist(df.out[index_i_1, disp_names]), nrow=private$n, ncol=private$n)

        # data:update to update from prior to posterior
        inputs = datmat[i,private$input.names]
        diag(V.obs) = obs_f(x0,inputs)
        y = datmat[i+1,private$obs.names]
        nas = !is.na(y)
        s = sum(nas)
        if (s > 0){
          ynew = y[nas]
          E. = diag(private$n)[nas,,drop=FALSE]
          e. = ynew - E. %*% h(x0,inputs)
          C. = E. %*% dhdx(x0,inputs)
          V. = E. %*% V.obs %*% t(E.)
          R. = C. %*% p0 %*% t(C.) + V.
          Ri. = solve(R.)
          K. = p0 %*% t(C.) %*% Ri.
          # update states
          x0  = x0 + K. %*% e.
          p0  = (I.n - K. %*% C.) %*% p0 %*% t(I.n - K. %*% C.) + K. %*% V. %*% t(K.)
        }
      }
      # convert to correlations instead of covariances
      if(!covariance & return.state.dispersion){
        # convert covariance to correlation
        df.out[,disp_names] = t(apply(df.out[,disp_names],
                                      MARGIN=1,
                                      FUN=function(x){
                                        cov = matrix(x,nrow=2,ncol=2)
                                        cor = stats::cov2cor(cov)
                                        diag(cor) = diag(cov)
                                        return(as.vector(cor))
                                      }
        ))
      }

      # should we return covariance/correlation stats?
      if(!return.state.dispersion){
        bool = !stringr::str_detect(names(df.out),"cor") & !stringr::str_detect(names(df.out),"cov")
        df.out = df.out[,bool]
      }

      # should 1:n.step.ahead be provided, or only n.step.ahead?
      if(give.only.n.step.ahead){
        df.out = df.out[df.out[["n.step.ahead"]]==n.step.ahead,]
      }

      # return data.frame
      return(pred=df.out)
    }

  }


  ############################################################
  ############################################################
  ##################### TMB-STYLE ############################
  ############################################################
  ############################################################

  if (private$method == "tmb") {
    message("'predict' is not available for method 'tmb' yet.")
    return(NULL)
  }



}

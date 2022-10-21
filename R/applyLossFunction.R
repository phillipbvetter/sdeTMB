applyLossFunction = function(data,model,control){
  
  if(control$loss_function=="identity"){
    data$loss_int = 0
    data$huber_c = control$c
    data$tukeypars = rep(1,4)
  }
  if(control$loss_function=="tukey"){
    data$loss_int = 1
    data$huber_c = control$c
    # Compute tukey coefficients
    rtukey = seq(0,100,by=1e-2)
    ctukey = control$c
    funtukey = function(r){
      ifelse(r^2 <= ctukey^2,
             ctukey^2/6 * (1-(1-(r/ctukey)^2)^3),
             ctukey^2/6
      )
    }
    tukeyloss = function(pars){
      res = sum((funtukey(rtukey) - pars[4]*(sigmoid(rtukey,a=pars[1],b=pars[2])+pars[3]))^2)
    }
    tukeyopt = nlminb(start=rep(1,4),objective=tukeyloss)
    data$tukeypars = tukeyopt$par
  }
  if(control$loss_function=="huber"){
    data$loss_int = 2
    data$huber_c = control$c
    data$tukeypars = rep(1,4)
  }
  
  return(data)
}

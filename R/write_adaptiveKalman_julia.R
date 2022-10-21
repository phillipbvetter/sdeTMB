write_adaptiveKalman_julia = function(model,data,silent) {
  
  #Get info
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
  }
  
  # Find time-dep inputs
  timedepInput = sort(unique(unlist(lapply(sdeEq,all.vars))))
  timedepInput = timedepInput[!timedepInput %in% diff.processes]
  timedepInput = timedepInput[!timedepInput %in% paste("d",state,sep="")]
  timedepInput = timedepInput[!timedepInput %in% names(c(data$pars, data$constants))]
  timedepInput = timedepInput[!timedepInput %in% state]
  
  #loading necessary julia packages
  julia_eval("using DifferentialEquations, Zygote, LinearAlgebra, Optimization, Optim, OptimizationOptimJL")
  
  # CONSTRUCTING FIRST MOMENT ODES
  ########################################
  ############ F-VECTOR ##################
  ########################################
  f_u = c()
  for(i in 1:n){
    f_u[i] = deparse(diff.terms[[i]]$dt)
    for(j in 1:n){
      f_u[i] = gsub(pattern=sprintf("\\b%s\\b",state[j]),replacement=sprintf("u[%s]",j),f_u[i])
    }
  }
  
  f_lines = c()
  f_lines[1] = "function compute_f(u,p,t)"
  f_lines[2] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  for(i in 1:n){
    f_lines[2+i] = sprintf("f%s = %s",i,f_u[i])
  }
  f_lines[n+3] = paste("return(vcat(",paste(sprintf(rep("f%s",n),1:n),collapse=","),"))",sep="")
  f_lines[n+4] = "end"
  f_complete = paste(f_lines,collapse=";")
  julia_eval(f_complete)
  
  ########################################
  ############ A-MATRIX ##################
  ########################################
  
  dfdx = matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    terms = lapply(state,function(x) Simplify(D(diff.terms[[i]]$dt,x)))
    for(j in 1:n){
      dfdx[i,j] = paste(deparse(terms[[j]]),collapse="")
      for(k in 1:n){
        dfdx[i,j] = gsub(pattern=sprintf("\\b%s\\b",state[k]),replacement=sprintf("u[%s]",k),dfdx[i,j])
      }
    }
  }
  # write function as strings
  dfdx_lines = c()
  dfdx_lines[1] = "function compute_A(u,p,t)"
  dfdx_lines[length(dfdx_lines)+1] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  count = 1
  for(i in 1:n){
    for(j in 1:n){
      dfdx_lines[length(dfdx_lines)+1] = sprintf("P%s = %s",count,dfdx[j,i])
      count = count + 1
    }
  }
  dfdx_lines[length(dfdx_lines)+1] = sprintf("P = reshape([%s],%s,%s)",paste(sprintf(rep("P%s",count-1),1:(count-1)),collapse=","),n,n)
  dfdx_lines[length(dfdx_lines)+1] = "return(P)"
  dfdx_lines[length(dfdx_lines)+1] = "end"
  dfdx_complete = paste(dfdx_lines,collapse=";")
  julia_eval(dfdx_complete)
  
  ########################################
  ############ G-MATRIX ##################
  ########################################
  
  g = matrix(0,nrow=n,ncol=ng)
  for(i in 1:n){
    for(j in 1:ng){
      g[i,j] = paste(deparse(diff.terms[[i]][[j+1]]),collapse = "")
      for(k in 1:n){
        g[i,j] = gsub(pattern=sprintf("\\b%s\\b",state[k]),replacement=sprintf("u[%s]",k),g[i,j])
      }
    }
  }
  gmat_lines = c()
  gmat_lines[1] = "function compute_G(u,p,t)"
  gmat_lines[length(gmat_lines)+1] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  count = 1
  for(i in 1:ng){
    for(j in 1:n){
      gmat_lines[length(gmat_lines)+1] = sprintf("G%s = %s",count,g[j,i])
      count = count + 1
    }
  }
  gmat_lines[length(gmat_lines)+1] = sprintf("G = reshape([%s],%s,%s)",paste(sprintf(rep("G%s",count-1),1:(count-1)),collapse=","),n,ng)
  gmat_lines[length(gmat_lines)+1] = "return(G)"
  gmat_lines[length(gmat_lines)+1] = "end"
  gmat_complete = paste(gmat_lines,collapse=";")
  julia_eval(gmat_complete)
  
  ########################################
  ############ ODE FUNCTION ##############
  ########################################
  
  ode_lines = c()
  ode_lines[1] ="function kalman_ode(u,p,t)"
  ode_lines[2] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  ode_lines[3] ="du_f = compute_f(u,p,t)"
  ode_lines[4] ="A = compute_A(u,p,t)"
  ode_lines[5] ="G = compute_G(u,p,t)"
  ode_lines[6] = sprintf("P = reshape(u[(%s+1):end],%s,%s)",n,n,n)
  ode_lines[7] ="du_P = vec(A*P + P*transpose(A) + G * transpose(G))"
  ode_lines[8] ="return(vcat(du_f,du_P))"
  ode_lines[9] ="end"
  ode_lines_complete = paste(ode_lines,collapse=";")
  julia_eval(ode_lines_complete)
  
  ########################################
  ############ OBSERVATION FUNCTION ######
  ########################################
  
  h = c()
  for(i in 1:m){
    h[i] = paste(deparse(obsEq[[i]][[1]][[3]]),collapse = "")
    for(j in 1:n){
      h[i] = gsub(pattern=sprintf("\\b%s\\b",state[j]),replacement=sprintf("u[%s]",j),h[i])
    }
  }
  h_lines = c()
  h_lines[1] = "function h(u,p)"
  h_lines[length(h_lines)+1] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  for(i in 1:m){
    h_lines[length(h_lines)+1] = sprintf("h%s = %s",i,h[i])
  }
  h_lines[length(h_lines)+1] = sprintf("h = [%s]",paste(sprintf(rep("h%s",m),1:m),collapse=","))
  h_lines[length(h_lines)+1] = "return(h)"
  h_lines[length(h_lines)+1] = "end"
  h_complete = paste(h_lines,collapse=";")
  julia_eval(h_complete)
  
  ########################################
  ############ JAC OF OBS FUNC ###########
  ########################################
  dhdx = matrix(NA,nrow=m,ncol=n)
  for(i in 1:m){
    terms = lapply(state, function(x) Simplify(D(obsEq[[i]][[1]][[3]],x)))
    for(j in 1:n){
      dhdx[i,j] = paste(deparse(terms[[j]]),collapse = "")
      # dhdx[i,j] = paste(format(deparse(terms[[j]]),nsmall=1),collapse = "")
      for(k in 1:n){
        dhdx[i,j] = gsub(pattern=sprintf("\\b%s\\b",state[k]),replacement=sprintf("u[%s]",k),dhdx[i,j])
      }
    }
  }
  # write function as strings
  dhdx_lines = c()
  dhdx_lines[1] = "function compute_C(u,p)"
  dhdx_lines[length(dhdx_lines)+1] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  count = 1
  for(i in 1:m){
    for(j in 1:n){
      dhdx_lines[length(dhdx_lines)+1] = sprintf("C%s = %s",count,dhdx[i,j])
      count = count + 1
    }
  }
  dhdx_lines[length(dhdx_lines)+1] = sprintf("C = reshape([%s],%s,%s)",paste(sprintf(rep("C%s",count-1),1:(count-1)),collapse=","),m,n)
  dhdx_lines[length(dhdx_lines)+1] = "return(C)"
  dhdx_lines[length(dhdx_lines)+1] = "end"
  dhdx_complete = paste(dhdx_lines,collapse=";")
  julia_eval(dhdx_complete)
  
  
  ########################################
  ############ MISSING OBS MATRIX ########
  ########################################
  
  E_lines = c()
  E_lines[1] = "function construct_E(s,m,nans)"
  E_lines[2] = "E = zeros(s,m)"
  E_lines[3] = "j = 1"
  E_lines[4] = "for i in 1:m"
  E_lines[5] = "if nans[i]"
  E_lines[6] = "E[j,i] = 1.0"
  E_lines[7] = "j += 1"
  E_lines[8] = "end"
  E_lines[9] = "end"
  E_lines[10] = "return(E)"
  E_lines[11] ="end"
  E_complete = paste(E_lines,collapse=";")
  julia_eval(E_complete)
  
  ########################################
  ############ KALMAN FUNCTION ###########
  ########################################
  
  varobs = paste(unlist(lapply(model$obsVar,function(x) deparse(x[[1]]))),collapse=", ")
  
  kalman_lines = c()
  kalman_lines[1] = "function kalman_nll(p,hyperp)"
  kalman_lines[length(kalman_lines)+1] = paste(paste(names(data$pars),collapse=","),"=p",sep="")
  kalman_lines[length(kalman_lines)+1] = "u0 = hyperp[1]"
  kalman_lines[length(kalman_lines)+1] = "tvec = hyperp[2]"
  kalman_lines[length(kalman_lines)+1] = "data = hyperp[3]"
  kalman_lines[length(kalman_lines)+1] = sprintf("V = diagm([%s])",varobs)
  kalman_lines[length(kalman_lines)+1] = "nll = 0"
  kalman_lines[length(kalman_lines)+1] = "prob = ODEProblem(kalman_ode,u0,[tvec[1],tvec[2]],p)"
  kalman_lines[length(kalman_lines)+1] = "for i in 1:(length(tvec)-1)"
  kalman_lines[length(kalman_lines)+1]  = "tspan = [tvec[i],tvec[i+1]]"
  kalman_lines[length(kalman_lines)+1]  = "prob = remake(prob,u0=u0,tspan=tspan,saveat=[tspan[2]])"
  kalman_lines[length(kalman_lines)+1]  = "sol = solve(prob,AutoTsit5(Rodas5()))"
  kalman_lines[length(kalman_lines)+1]  = sprintf("x0 = sol[1][1:%s]",n)
  kalman_lines[length(kalman_lines)+1]  = sprintf("p0 = reshape(sol[1][(%s+1):end],%s,%s)",n,n,n)
  kalman_lines[length(kalman_lines)+1]  = "y = data[i,:]"
  kalman_lines[length(kalman_lines)+1]  = "nans = .!isnan.(y)"
  kalman_lines[length(kalman_lines)+1]  = "s = sum(nans)"
  kalman_lines[length(kalman_lines)+1]  = "if s > 0.5"
  kalman_lines[length(kalman_lines)+1]  = sprintf("E = construct_E(s,%s,nans)",m)
  kalman_lines[length(kalman_lines)+1]  = "e = y[nans] - h(x0,p)[nans]"
  kalman_lines[length(kalman_lines)+1]  = "C = E * compute_C(x0,p)"
  kalman_lines[length(kalman_lines)+1]  = "V2 = E * V * E"
  kalman_lines[length(kalman_lines)+1]  = "R = C * p0 * transpose(C) + V2"
  kalman_lines[length(kalman_lines)+1]  = "Ri = inv(R)"
  kalman_lines[length(kalman_lines)+1]  = "K = p0 * transpose(C) * Ri"
  kalman_lines[length(kalman_lines)+1]  = "x0 = x0 + K * e"
  kalman_lines[length(kalman_lines)+1]  = "p0 = (I-K*C) * p0 * transpose(I-K*C) + K*V2*transpose(K)"
  kalman_lines[length(kalman_lines)+1]  = "nll += 0.5*logdet(R) + 0.5*transpose(e)*Ri*e + 0.5*log(2*pi)*s"
  kalman_lines[length(kalman_lines)+1]  = "end"
  kalman_lines[length(kalman_lines)+1]  = "u0 = vec(vcat(x0,reshape(p0,:,1)))"
  kalman_lines[length(kalman_lines)+1]  = "end"
  kalman_lines[length(kalman_lines)+1]  = "return(nll)"
  kalman_lines[length(kalman_lines)+1]  = "end"
  kalman_complete=paste(kalman_lines,collapse=";")
  julia_eval(kalman_complete)
  
  callback_lines = c()
  callback_lines[1] = "function cb(pars,loss)"
  callback_lines[length(callback_lines)+1] = "@show loss"
  callback_lines[length(callback_lines)+1] = "return(false)"
  callback_lines[length(callback_lines)+1] = "end"
  callback_complete = paste(callback_lines,collapse=";")
  julia_eval(callback_complete)
  
  opt_lines = c()
  opt_lines[1] = sprintf("function %s(p,hyperp)",model$modelname)
  opt_lines[length(opt_lines)+1] = "fun = OptimizationFunction(kalman_nll,Optimization.AutoForwardDiff());"
  opt_lines[length(opt_lines)+1] = "prob = OptimizationProblem(fun,p,hyperp);"
  if(silent){
    opt_lines[length(opt_lines)+1] = "opt = solve(prob,Optim.NewtonTrustRegion())"  
  } else {
    opt_lines[length(opt_lines)+1] = "opt = solve(prob,Optim.NewtonTrustRegion(),callback=cb)"
  }
  opt_lines[length(opt_lines)+1] = "return([opt.u, opt.minimum, opt.retcode, opt.original.iterations, opt.original.f_calls, opt.original.g_calls, opt.original.h_calls])"
  opt_lines[length(opt_lines)+1] = "end"
  opt_complete = paste(opt_lines,collapse=";")
  julia_eval(opt_complete)
}

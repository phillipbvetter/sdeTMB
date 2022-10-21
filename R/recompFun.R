# Function for recompilation
recompFun = function(compile,method,model,data,control,pathdir){
  modelname_with_extension = paste(model$modelname2,".cpp",sep="")
  if(compile){
    switch(method,
           tmb = write_TMB(model,data,control),
           ekf = write_ExtendedKalman(model,data,control),
           ukf = write_UnscentedKalman(model,data,control),
           tmb_exact = write_TMBexact(model,data,control)
    )
    compile(modelname_with_extension)
    #reload the library
    try(dyn.unload(dynlib(model$modelname2)),silent=T)
    try(dyn.load(dynlib(model$modelname2)),silent=T)
  }
  if(!compile & !file.exists(modelname_with_extension)){
    print("You asked me not to compile, but the model C++ file doesn't exist so I will compile anyway")
    recompFun(TRUE,method,model,data)
  }
  #load the library
  try(dyn.load(dynlib(model$modelname2)),silent=T)
}

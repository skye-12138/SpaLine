.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
  env<-reticulate::conda_list()
  if("SpaLine_conda" %in% env$name){
    reticulate::use_condaenv("SpaLine_conda", required = TRUE)
  }else{
    reticulate::conda_create("SpaLine_conda")
    reticulate::use_condaenv("SpaLine_conda", required = TRUE)
  }
  reticulate::py_config()
}

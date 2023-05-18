.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
  reticulate::conda_create("SpaLine_conda")
  reticulate::use_condaenv("SpaLine_conda", required = TRUE)
  reticulate::py_config()
}

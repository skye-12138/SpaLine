.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
  reticulate::conda_create("SpaLine_conda",conda ="/Users/skye/miniforge3/bin/conda")
  reticulate::use_condaenv("SpaLine_conda",conda = "/Users/skye/miniforge3/bin/conda")
  reticulate::py_config()
}

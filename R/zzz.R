.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
  reticulate::conda_create("SpaLine_conda")
  reticulate::use_condaenv("SpaLine_conda", required = TRUE)
  reticulate::py_config()
  reticulate::py_run_file(system.file("python","st_nucleus.py",package = "SpaLine"))
  reticulate::py_run_file(system.file("python","st_segment.py",package = "SpaLine"))
}

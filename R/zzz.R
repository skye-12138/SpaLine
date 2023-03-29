.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
  reticulate::py_run_file(system.file("python","st_nucleus.py",package = "SpaLine"))
  reticulate::py_run_file(system.file("python","st_segment.py",package = "SpaLine"))
}

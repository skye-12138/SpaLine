.onAttach<-function(libname,pkgname){
  reticulate::configure_environment(pkgname)
}

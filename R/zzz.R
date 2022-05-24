# Rcpp::loadModule("double_cpp", TRUE)
.onUnload <- function (libpath) {
  library.dynam.unload("ergm.sign", libpath)
}
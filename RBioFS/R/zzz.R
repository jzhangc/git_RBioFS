.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.

  To cite in publication:
  Zhang J, Hadj-Moussa H, Storey KB. 2016. Current progress of high-throughput microRNA differential
  expression analysis and random forest gene selection for model and non-model systems: an R implementation.
  J Integr Bioinform. 13: 306.

  For more details, please visit: http://kenstoreylab.com/?page_id=2542")
  return(TRUE)
}

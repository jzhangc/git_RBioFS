#' @title uni_PreProc
#'
#' @description a legacy function borrowed from RBioArray. Data pre-processing function for the microarary data, which is now used for univariate reduction during svm nested CV.
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perform a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransParallelComputing If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @details The function does not use design matrix for array weight calculation.
#'          Therefore, the DE analysis based on the output from this function will yield slightly different results from the \code{\link{rbioarray_filter_combine}}.
#'
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values.
#'         If \code{logTrans = TRUE}, the function also outputs a \code{csv} file containing the log transformed data.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @export
uni_PreProc <- function(rawlist, logTrans = FALSE, logTransMethod = "log2",
                        logTransObjT = "data", logTransParallelComputing = FALSE,
                        bgMethod = "auto", normMethod = "quantile", ...){
  if (class(rawlist) == "list"){
    ## log transform  or not
    if (logTrans){
      if (!logTransParallelComputing){
        # log transform
        mtx <- apply(rawlist$E, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
      } else {
        # parallel computing
        # set up cpu cluster
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores, type = "PSOCK")
        on.exit(stopCluster(cl)) # close connect when exiting the function
        # log transform
        mtx <- foreach(i = rawlist$E) %dopar% {
          out <- apply(i, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
        }
      }
      tmpdata <- list(E = mtx, genes = rawlist$genes, targets = rawlist$targets)
      # store and export log transformed data into a csv file
      logTransOut <- data.frame(rawlist$genes, mtx)
      write.csv(logTransOut, file = paste(logTransObjT, "_log_transformed.csv", sep = ""), row.names = FALSE)
    } else {
      tmpdata <- rawlist
    }

    ## normalization
    BgC <- backgroundCorrect.matrix(tmpdata$E, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm) # array weight
    output <- list(E = Norm, genes = rawlist$genes, targets = rawlist$targets, ArrayWeight = Wgt)
  } else {
    ## normalization
    BgC <- backgroundCorrect(rawlist, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm)
    Norm$ArrayWeight <- Wgt
    output <- Norm
  }
  return(output)
}


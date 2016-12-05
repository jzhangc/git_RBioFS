#' @title singleIMP
#'
#' @description data imputation for a vector of values.
#' @param tgtVar Input single target variable.
#' @param methd Imputation method. Options are \code{"mean"}, \code{"random"} and \code{"regression"}. Default is \code{"mean"}.
#' @return Outputs a \code{vector} object with all the missing value imputed.
#' @examples
#' \dontrun{
#' mydfm <- singleIMP(myVector, methd = "mean")
#' }
#' @export
singleIMP <- function(tgtVar, methd = "mean"){
  missing <- is.na(tgtVar)
  nMissing <- sum(missing)
  completeObs <- tgtVar[!missing]
  imp <- tgtVar

  if (methd == "mean"){
    imp[missing] <- mean(completeObs)
  } else if (methd == "random"){
    imp[missing] <- sample(completeObs, size = nMissing, replace = TRUE)
  } else if (methd == "regression"){
    return("to be completed")
  } else {
    stop("please choose the proper imputing methods")
  }

  return(imp)
}


#' @title rbioIMP
#'
#' @description data imputation for raw data based on the factor.
#' @param dfm Input dataframe. Case sensitive and be sure to type with quotation marks.
#' @param fct Input factor variable, based on the which the missing data is imputed.
#' @param annot Input variable or character string used for variable name for the transposed dataframe.
#' @param method Imputation method. Options are \code{"mean"}, \code{"random"} and \code{"rf"}. Default is \code{"mean"}. When \code{rf}, it uses the data impute function from \code{randomForest} pacakge.
#' @param transpo If the output file is transposed. Default is \code{TRUE}.
#' @param ... See \code{rfImpute()} funtion in \code{randomForest} package for details.
#' @return Outputs a \code{dataframe} object with all the missing value imputed.
#' @details Make sure to make the input the data with only the tagert variable columns, meaning no annotation or index variables. The format would be target variables only (column). \code{transpo} arugment is particularly useful for the upcoming data normalization and random forest operations.
#' @importFrom randomForest rfImpute
#' @examples
#' \dontrun{
#' mydfm <- rbioIMP(raw[-c(1:2)], raw$Conditions, raw$sampleName, method = "rf", transpo = TRUE, iter = 10, ntree = 501) # make sure no annotation variable present in the data file.
#' }
#' @export
rbioIMP <- function(dfm, fct, annot, method = "mean", transpo = FALSE, ...){

  out <- dfm
  if (method == "rf"){
    out <- rfImpute(dfm, fct, ...)
    out <- out[, -1, drop = FALSE]

  } else {
    # ave() function applies a function to a dataframe by factor. Make sure to use write out FUN
    # by some R magical property, use [] preserves the dataframe format for object when using lapply()
    out[] <- lapply(out,
                    function(i)ave(i, fct,
                                   FUN = function(j)singleIMP(j, methd = method))
    )
  }

  out <- as.data.frame(out, check.names = FALSE)

  if (transpo == TRUE){
    out <- t(out)
    colnames(out) <- annot
  }

  return(out)
}

#' @title rbioNorm
#'
#' @description Data normalization. This is a shell function using limma package.
#' @param RawData Inpule dataframe, which can be the the results from the function \code{\link{rbioIMP}} with transposition.
#' @param NormMtd Normalization method. Options are \code{"none"}, \code{"scale"}, \code{"quantile"} and \code{"cyclicloess"}. Default is \code{"quantile"}. See \code{limma} package for details.
#' @param correctBG If to correct background. Default is \code{TRUE}.
#' @param BgMtd Background correction method if \code{correctBG} is TRUE. Options are \code{"auto"}, \code{"none"}, \code{"subtract"}, \code{"half"}, \code{"minimum"}, \code{"movingmin"}, \code{"edwards"} and \code{"normexp"}. Default is \code{"normexp"}. See \code{limma} package for details.
#' @param BgOffst Background offset if \code{correctBG} is TRUE. Default is \code{0}. Briefly, this is the constant value to add to the data. See \code{limma} package for details.
#' @return Outputs a \code{dataframe} with processed data.
#' @importFrom limma backgroundCorrect normalizeBetweenArrays
#' @examples
#' \dontrun{
#' normdata <- rbioNorm(rawT, correctBG = FLASE)
#' }
#' @export
rbioNorm<-function(RawData, NormMtd = "quantile",
                      correctBG = TRUE,
                      BgMtd = "normexp", BgOffst = 0){

  if (correctBG == TRUE){
    BgC <- backgroundCorrect(RawData,method = BgMtd, offset = BgOffst) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = NormMtd) # quantile normalization
  } else if (correctBG == FALSE){
    Norm<-normalizeBetweenArrays(RawData, method = NormMtd) # quantile normalization
  } else {
    stop("Please check the arguments")
  }

  out <- as.data.frame(Norm, check.names = FALSE)

  return(out)
}


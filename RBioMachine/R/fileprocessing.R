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
#' @description data imputation for raw data.
#' @param fileName Inpule file name. Case sensitive and be sure to type with quotation marks. Currently only takes \code{.csv} files.
#' @param method Imputation method. Options are \code{"mean"}, \code{"random"} and \code{"regression"}. Default is \code{"mean"}.
#' @param transpo If the output file is transposed. Default is \code{TRUE}.
#' @return Outputs a \code{dataframe} object with all the missing value imputed.
#' @details Make sure to make the input the data with only the tagert variable columns, meaning no annotation or index variables. The format would be target variables only (column). \code{transpo} arugment is particularly useful for the upcoming data normalization and random forest operations.
#' @examples
#' \dontrun{
#' mydfm <- rbioIMP("mydata.csv", method = "mean") # make sure no aannotation variable present in the data file.
#' }
#' @export
rbioIMP <- function(fileName, method = "mean", transpo = TRUE){
  raw <- read.csv(file = fileName, na.strings = " ", stringsAsFactors = FALSE)
  impdata <- as.data.frame(sapply(raw, singleImp, methd = method)) # impute the missing values
  return(impdata)
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

  out <- as.data.frame(Norm)

  return(out)
}


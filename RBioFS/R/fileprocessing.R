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
rbioNorm <- function(RawData, NormMtd = "quantile",
                      correctBG = TRUE,
                      BgMtd = "normexp", BgOffst = 0){
  if (correctBG == TRUE){
    BgC <- backgroundCorrect(RawData,method = BgMtd, offset = BgOffst) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = NormMtd) # quantile normalization
  } else if (correctBG == FALSE){
    Norm <- normalizeBetweenArrays(RawData, method = NormMtd) # quantile normalization
  } else {
    stop("Please check the arguments")
  }

  out <- as.data.frame(Norm, check.names = FALSE)
  return(out)
}


#' @title center_scale
#'
#' @description data centering function with scale opiton
#' @param x Input matrix. Make sure it is a matrix. Typically, row is samples, and column the features.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @return Outputs a list containing centered matrix.
#' @details This function is needed to pre-process data when conducting plsda analysis.
#'          And the function can also be used for other purposes when needed.
#'          With the colSds function from matrixStats pacakge to calcualte column standard deviation,
#'          this function is faster than the native \code{scale} function when \code{scale = TRUE}.
#' @importFrom matrixStats colSds
#' @examples
#' \dontrun{
#' centerX <- center_scale(data, scale = TRUE)
#' }
#' @export
center_scale <- function(x, scale = TRUE){
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (class(x) == "data.frame" | is.vector(x)){
    x <- as.matrix(sapply(x, as.numeric))
  }

  col.mean <- colMeans(x, na.rm = TRUE)
  if (scale){
    col.sd <- colSds(x, center = col.mean, na.rm = TRUE) # matrixStats::colSds
    # check zero values
    if (any(abs(col.sd) < .Machine$double.eps^0.5)) warning("Scaling with (near) zero standard deviation")
    temp <- t((t(x) - col.mean) / col.sd)  # centre + scale
  } else {
    col.sd <- NULL
    temp <- t((t(x) - col.mean)) # centre
  }
  out <- list(centerX = temp, meanX = col.mean, scale = ifelse(scale, TRUE, FALSE), columnSD = col.sd)
  return(out)
}


#' @title fs_csv_generator
#'
#' @description Creat CSV file from microarray or RNAseq experiments ready for RBioFS
#' @param objTitle Prefix for the output file name. Default is \code{"data"}.
#' @param fltdata Filtered expression data from microarray or RNAseq experiments.
#' @param rmControl Only functional for Agilent platform microarray datasets, to remove control probes. For other datasets, set to \code{FALSE}. Default is \code{FALSE}.
#' @param anno Gene annotation dataframe.
#' @param geneSymbol.var Name for the gene symbol variable from \code{anno}. Default is \code{"GeneSymbol"}.
#' @param tgt Sample annotation data frame.
#' @param group.var Name for the sample group variable from from \code{tgt}. Default is \code{"Group"}.
#' @param de.dfm Differetinal expression resutls dataframe.
#' @param probeID.var Name for the probe ID variable from \code{de.dfm}. Default is \code{"ProbeName"}.
#' @param fc.var Name for the fold change variable from \code{de.dfm}. Default is \code{"logFC"}.
#' @param p.var Name for the raw p value variable from \code{de.dfm}. Defaults is \code{"P.Value"}.
#' @param adjP.var Name for the adjusted value variable from \code{de.dmf}. Default is \code{"adj.P.Val"}.
#' @param q.value P value threshold for subsetting the filtered data. Default is \code{0.05}.
#' @param fc.value Fold change threshold for subsetting filtered data. Default is \code{1.5}.
#' @param fc.log.trans Wether the input fold change is log2 transformed. Default is \code{TRUE}.
#' @details The function works best with the results from RBioArray.
#' @return Outputs a \code{csv} file to working directory for RBioFS. If not \code{Elist} class from limma, Make sure that \code{fltdata} is a list object with \code{genes} (name) dataframe for gene annotation and \code{E} (name) as expression matrix/dataframe. For \code{tgt}, make sure to have sample name as rownames.
#' @examples
#' \dontrun{
#' fs_csv_generator(objTitle = "test", geneSymbol.var = "Gene.Symbol", fltdata = fltdata, anno = annot,
#' de.dfm = diana_medaka_DE$HTvC, tgt = Tgt, group.var = "Treatment")
#' }
#' @export
fs_csv_generator <- function(objTitle = "data",
                             fltdata = NULL, rmControl = FALSE,
                             anno = NULL, geneSymbol.var = "GeneSymbol",
                             tgt = NULL, group.var = "Group",
                             de.dfm = NULL, probeID.var = "ProbeName", fc.var = "logFC", p.var = "P.Value", adjP.var = "adj.P.Val",
                             q.value = 0.05, fc.value = 1.5, fc.log.trans = TRUE){

  ## check all the variables
  if (is.null(fltdata)){
    stop("Filtered data not found. Function aborted.")
  }
  if (is.null(anno)){
    stop("Annotation dataframe not found. Function aborted.")
  }
  if (is.null(de.dfm)){
    stop("DE data not found. Function aborted.")
  }
  if (is.null(tgt)){
    stop("Target index dataframe not found. Function aborted.")
  }
  if (!geneSymbol.var %in% names(anno)){
    stop("Gene symbol variable not found. Function aborted.")
  }
  if (!probeID.var %in% names(anno)){
    stop("Probe ID variable not found in the annotation data. Function aborted.")
  }
  if (!probeID.var %in% names(de.dfm)){
    stop("Probe ID variable not found in the DE data. Function aborted.")
  }
  if (FALSE %in% (c(p.var, adjP.var,fc.var) %in% names(de.dfm))){
    stop("One of the required variables not found in the DE data. Function aborted.")
  }
  if (!group.var %in% names(tgt)){
    stop("Group variable not found in the target index. Function aborted.")
  }

  ## start
  dfm <- data.frame(fltdata$genes, fltdata$E)
  if (rmControl){
    dfm <- dfm[dfm$ControlType == 0, ]
  }
  if (!probeID.var %in% names(dfm)){ # check the variable for dfm
    stop("Probe ID variable not found in the filtered data. Function aborted.")
  }

  # threosholds
  p <- max(de.dfm[de.dfm[, adjP.var] < 0.05, p.var])

  if (fc.log.trans){
    fc <- log2(fc.value)
  } else {
    fc <- fc.value
  }

  # function core
  tmp <- dfm[dfm[, probeID.var] %in% de.dfm[de.dfm[, p.var] < p & abs(de.dfm[, fc.var]) > fc, probeID.var], ]
  tmp <- merge(anno, tmp, by = probeID.var, all.y = TRUE)
  tmp <- tmp[complete.cases(tmp[, geneSymbol.var]), ]

  out <- tmp[, - c(1:(ncol(tmp) - ncol(fltdata$E)))]
  out <- t(out)
  colnames(out) <- tmp[, geneSymbol.var]
  out <- data.frame(Sample = rownames(out), Group = Tgt[rownames(out) %in% rownames(tgt), group.var], out)

  ## output
  write.csv(out, paste(objTitle, ".fs_ready.csv", sep = ""), row.names = FALSE)
}


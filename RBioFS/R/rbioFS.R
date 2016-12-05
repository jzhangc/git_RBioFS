#' @title rbioFS
#'
#' @description Recursive nested random froest variable importance (vi) and OOB error rate computation in a sequential forward selection (SFS) manner.
#' @param objTitle The title for the output data frame. Default is \code{"x_vs_tgt"}
#' @param file Input file. Only takes \code{.csv} format.
#' @param impute Wether to use data imputation functionality to impute missing data points. Default is \code{FALSE}.
#' @param imputeMethod The method for data imputation, only set if \code{impute = TRUE}. Default is \code{"rf"}. See \code{\link{rbioIMP}} for details.
#' @param imputeIter RF imputation iteration value. Default is \code{10}. See \code{\link{rbioIMP}} for details.
#' @param imputeNtree RF imputation ntree value. Default is \code{501}. See \code{\link{rbioIMP}} for details.
#' @param quantileNorm Wehter ot use quantile nomalization on the raw data. Default is \code{FALSE}.
#' @param nTimes Number of random forest for both initial FS and SFS-like FS. Default is \code{50} times.
#' @param nTree Number of trees generated for each random forest run for both initial FS and SFS-like FS. Default is \code{1001} trees.
#' @param SFS_mTry Number of randomly selected featurs for constructing trees for SFS-like FS step. When \code{"recur_default"}, it'll be based on \code{p / 3}; when \code{"rf_default"}, it will use the default setting in \code{randomForest} package. Default is \code{"recur_default"}.
#' @param multicore If to use parallel computing. Default is \code{TRUE}.
#' @param plot If to plot results, bar graph for initial FS and joint-point curve for SFS-like FS. Default is \code{TRUE}
#' @param initalFS_Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param initalFS_xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param initalFS_yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param initalFS_n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param initalFS_errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param initalFS_errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param initalFS_xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param initalFS_yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param initalFS_plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param initalFS_plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @param SFS_xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param SFS_yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param SFS_n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param SFS_errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param SFS_errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param SFS_symbolSize Size of the symbol. Default is \code{2}.
#' @param SFS_xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param SFS_yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param SFS_plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param SFS_plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs two \code{list} objects and \code{.csv} for the two FS steps. And, if \code{plot = TRUE}, a bargraph for inital FS and joint-point curve for SFS-like FS step. A \code{.csv} file with imputed data is also generated if \code{impute = TRIE}.
#' @details Make sure to arrange input \code{.csv} file with first two columns for smaple ID and conditions, and the rest for features (e.g., genes).
#' @examples
#' \dontrun{
#' rbioFS(file = "test.csv", impute = TRUE, imputeIter = 50, quantileNorm = TRUE)
#' }
#' @export
rbioFS <- function(file, impute = FALSE, imputeMethod = "rf", imputeIter = 10, imputeNtree = 501,
                   quantileNorm = FALSE,
                   nTimes = 50, nTree = 1001, SFS_mTry = "recur_default",
                   multicore = TRUE,
                   plot = TRUE,
                   initalFS_n = "all",
                   initalFS_Title = NULL, initalFS_xLabel = "Mean Decrease in Accuracy", initalFS_yLabel = NULL,
                   initalFS_errorbar = "SEM", initalFS_errorbarWidth = 0.2,
                   initalFS_xTxtSize = 10, initalFS_yTxtSize =10,
                   initalFS_plotWidth = 170, initalFS_plotHeight = 150,
                   SFS_n = "all",
                   SFS_Title = NULL, SFS_xLabel = NULL, SFS_yLabel = NULL,
                   SFS_errorbar = "SEM", SFS_errorbarWidth = 0.2,
                   SFS_symbolSize = 2, SFS_xTxtSize = 10, SFS_yTxtSize =10,
                   SFS_plotWidth = 170, SFS_plotHeight = 150
                   ){
  raw <- read.csv(file = file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, check.names = FALSE)
  tgt <- factor(raw[[2]], levels = unique(raw[[2]])) # target variable

  if (impute){
    input <- RBioFS::rbioIMP(dfm = raw[-c(1:2)], method = imputeMethod,
                             iter = imputeIter, ntree = imputeNtree,
                             fct = tgt, annot = raw[[1]], transpo = FALSE)
    input <- data.frame(input, check.names = FALSE)

    out <- data.frame(raw[, c(1:2)], input)
    write.csv(out, file = paste(substr(noquote(file), 1, nchar(file) - 4), "_imputed.csv", sep = ""), row.names = FALSE)
  } else {
    input <- data.frame(raw[, -c(1:2)], check.names = FALSE)
  }

  if (quantileNorm){ # normalization
    input <- t(input)
    input <- RBioFS::rbioNorm(input, correctBG = FALSE)
    input <- t(input)
    input <- data.frame(input, check.names = FALSE)
  }

  ## FS
  if (plot){

    RBioFS::rbioRF_initialFS(objTitle = substr(noquote(file), 1, nchar(file) - 4), x = input, targetVar = tgt,
                             nTimes = nTimes, nTree = nTree,
                             plot = TRUE, n = initalFS_n,
                             Title = initalFS_Title,
                             errorbar = initalFS_errorbar, errorbarWidth = initalFS_errorbarWidth,
                             xLabel = initalFS_xLabel, yLabel = initalFS_yLabel,
                             xTxtSize = initalFS_xTxtSize, yTxtSize = initalFS_yTxtSize,
                             plotWidth = initalFS_plotWidth, plotHeight = initalFS_plotHeight) # initial FS


    RBioFS::rbioRF_SFS(objTitle = substr(noquote(file), 1, nchar(file) - 4),
                       x = get(paste(substr(noquote(file), 1, nchar(file) - 4), "_inital_FS", sep = ""))$matrix_initial_FS,
                       targetVar = tgt, nTimes = nTimes, mTry = SFS_mTry,
                       plot = TRUE,
                       n = SFS_n,
                       Title = SFS_Title, xLabel = SFS_xLabel, yLabel = SFS_yLabel,
                       errorbar = SFS_errorbar, errorbarWidth = SFS_errorbarWidth,
                       symbolSize = SFS_symbolSize, xTxtSize = SFS_xTxtSize, yTxtSize = SFS_yTxtSize,
                       plotWidth = SFS_plotWidth, plotHeight = SFS_plotHeight) # SFS



  } else {

    RBioFS::rbioRF_initialFS(objTitle = substr(noquote(file), 1, nchar(file) - 4), x = input, targetVar = tgt,
                             nTimes = nTimes, nTree = nTree,
                             plot = FALSE) # initial FS


    RBioFS::rbioRF_SFS(objTitle = substr(noquote(file), 1, nchar(file) - 4),
                       x = get(paste(substr(noquote(file), 1, nchar(file) - 4), "_inital_FS", sep = ""))$matrix_initial_FS,
                       targetVar = tgt, nTimes = nTimes, mTry = SFS_mTry,
                       plot = FALSE) # SFS

  }

}

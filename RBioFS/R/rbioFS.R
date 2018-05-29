#' @title rbioFS
#'
#' @description Recursive nested random forest variable importance (vi) and OOB error rate computation in a sequential forward selection (SFS) manner.
#' @param objTitle The title for the output data frame. Default is \code{"data"}
#' @param file Input file. Only takes \code{.csv} format. Set either this or \code{input}, but not both.
#' @param input Input data frame. They type should be \code{data frame} or \code{matrix}. Set either this or \code{file}, but not both.
#' @param sampleIDVar Sample variable name. It's a character string.
#' @param groupIDVar Group variable name. It's a character string.
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
#' @param initialFS_Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param initialFS_xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param initialFS_yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param initialFS_n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param initialFS_errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param initialFS_errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param initialFS_xTickLblSize Font size for the x-axis text. Default is \code{10}.
#' @param initialFS_yTickLblSize Font size for the y-axis text. Default is \code{10}.
#' @param initialFS_plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param initialFS_plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @param SFS_xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param SFS_yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param SFS_n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param SFS_errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param SFS_errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param SFS_symbolSize Size of the symbol. Default is \code{2}.
#' @param SFS_xTickLblSize Font size for the x-axis text. Default is \code{10}.
#' @param SFS_yTickLblSize Font size for the y-axis text. Default is \code{10}.
#' @param SFS_plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param SFS_plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs two \code{list} objects and \code{.csv} for the two FS steps. And, if \code{plot = TRUE}, a bargraph for initial FS and joint-point curve for SFS-like FS step. A \code{.csv} file with imputed data is also generated if \code{impute = TRIE}.
#' @details Make sure to arrange input \code{.csv} file with first two columns for smaple ID and conditions, and the rest for features (e.g., genes).
#' @examples
#' \dontrun{
#' rbioFS(file = "test.csv", impute = TRUE, imputeIter = 50, quantileNorm = TRUE)
#' }
#' @export
rbioFS <- function(objTitle = "data", file = NULL, input = NULL, sampleIDVar = NULL, groupIDVar = NULL,
                   impute = FALSE, imputeMethod = "rf", imputeIter = 10, imputeNtree = 501,
                   quantileNorm = FALSE,
                   nTimes = 50, nTree = 1001, SFS_mTry = "recur_default",
                   multicore = TRUE,
                   plot = TRUE,
                   initialFS_n = "all",
                   initialFS_Title = NULL, initialFS_xLabel = "Mean Decrease in Accuracy", initialFS_yLabel = NULL,
                   initialFS_errorbar = "SEM", initialFS_errorbarWidth = 0.2,
                   initialFS_xTickLblSize = 10, initialFS_yTickLblSize =10,
                   initialFS_plotWidth = 170, initialFS_plotHeight = 150,
                   SFS_n = "all",
                   SFS_Title = NULL, SFS_xLabel = NULL, SFS_yLabel = NULL,
                   SFS_errorbar = "SEM", SFS_errorbarWidth = 0.2,
                   SFS_symbolSize = 2, SFS_xTickLblSize = 10, SFS_yTickLblSize =10,
                   SFS_plotWidth = 170, SFS_plotHeight = 150
                   ){
  ## argument check
  if (!is.null(file) & !is.null(input)) stop("set only one of the \"file\" and \"input\"")
  if (is.null(file) & is.null(input)) stop("set one of the \"file\" and \"input\"")
  if (is.null(sampleIDVar) | is.null(groupIDVar)) stop("set both sampleIDVar and groupIDVar")

  ## input constuction
  if (is.null(input)){
    raw <- read.csv(file = file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    raw <- data.frame(input, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!all(c(sampleIDVar, groupIDVar) %in% names(raw))) stop("sampleIDvar and/or groupIDvar not found in the input dataframe.")

  tgt <- factor(raw[, groupIDVar], levels = unique(raw[, groupIDVar])) # target variable

  ## imputation
  if (impute){
    imp_data <- RBioFS::rbioIMP(dfm = raw[, !names(raw) %in% c(sampleIDVar, groupIDVar)], method = imputeMethod,
                             iter = imputeIter, ntree = imputeNtree,
                             fct = tgt, annot = raw[, sampleIDVar], transpo = FALSE)
    imp_data <- data.frame(imp_data, check.names = FALSE)
    fs_data <- imp_data

    # export imputation results
    out <- data.frame(raw[, c(sampleIDVar, groupIDVar)], imp_data, check.names = FALSE)
    write.csv(out, file = paste(objTitle, "_imputed.csv", sep = ""), row.names = FALSE)
  } else {
    fs_data <- data.frame(raw[, !names(raw) %in% c(sampleIDVar, groupIDVar)], check.names = FALSE)
  }

  if (quantileNorm){ # normalization
    fs_data <- t(fs_data)
    fs_data <- RBioFS::rbioNorm(fs_data, correctBG = FALSE)
    fs_data <- t(fs_data)
    fs_data <- data.frame(fs_data, check.names = FALSE)
  }

  ## FS
  if (plot){
    cat(paste("Initial selection with plotting...", sep = ""))  # initial message
    RBioFS::rbioRF_initialFS(objTitle = objTitle, x = fs_data, targetVar = tgt,
                             nTimes = nTimes, nTree = nTree,
                             plot = TRUE, n = initialFS_n,
                             plot.title = initialFS_Title,
                             plot.errorbar = initialFS_errorbar, plot.errorbarWidth = initialFS_errorbarWidth,
                             plot.xLabel = initialFS_xLabel, plot.yLabel = initialFS_yLabel,
                             plot.xTickLblSize = initialFS_xTickLblSize, plot.yTickLblSize = initialFS_yTickLblSize,
                             plot.Width = initialFS_plotWidth, plot.Height = initialFS_plotHeight) # initial FS
    cat(paste("Done!\n", sep = ""))  # final message

    cat(paste("Sequential forward selection with plotting...", sep = ""))  # initial message
    RBioFS::rbioRF_SFS(objTitle = objTitle,
                       x = get(paste(objTitle, "_initial_FS", sep = ""))$matrix_initial_FS,
                       targetVar = tgt, nTimes = nTimes, mTry = SFS_mTry,
                       plot = TRUE,
                       n = SFS_n,
                       plot.title = SFS_Title, plot.xLabel = SFS_xLabel, plot.yLabel = SFS_yLabel,
                       plot.errorbar = SFS_errorbar, plot.errorbarWidth = SFS_errorbarWidth,
                       plot.symbolSize = SFS_symbolSize, plot.xTickLblSize = SFS_xTickLblSize, plot.yTickLblSize = SFS_yTickLblSize,
                       plot.Width = SFS_plotWidth, plot.Height = SFS_plotHeight) # SFS
    cat(paste("Done!\n", sep = ""))  # final message


  } else {
    cat(paste("Initial selection without plotting...", sep = ""))  # initial message
    RBioFS::rbioRF_initialFS(objTitle = objTitle, x = fs_data, targetVar = tgt,
                             nTimes = nTimes, nTree = nTree,
                             plot = FALSE) # initial FS
    cat(paste("Done!\n", sep = ""))  # final message

    cat(paste("Sequential forward selection without plotting...", sep = ""))  # initial message
    RBioFS::rbioRF_SFS(objTitle = objTitle,
                       x = get(paste(objTitle, "_initial_FS", sep = ""))$matrix_initial_FS,
                       targetVar = tgt, nTimes = nTimes, mTry = SFS_mTry,
                       plot = FALSE) # SFS
    cat(paste("Done!\n", sep = ""))  # final message
  }
}

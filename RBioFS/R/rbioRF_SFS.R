#' @title rbioRF_SFS
#'
#' @description Recursive nested random froest variable importance (vi) and OOB error rate computation in a sequential forward selection (SFS) manner.
#' @param objTitle The title for the output data frame. Default is \code{"x_vs_tgt"}
#' @param x Input dataframe or matrix. Make sure to arrange the data with features as column names.
#' @param targetVar The target variable for random forest feature selection. This is a factor object.
#' @param nTimes Number of random forest vi computation runs. Default is \code{50} times.
#' @param nTree Number of trees generated for each random forest run. Default is \code{1001} trees.
#' @param mTry Number of randomly selected featurs for constructing trees. When \code{"recur_default"}, it'll be based on \code{p / 3}; when \code{"rf_default"}, it will use the default setting in \code{randomForest} package. Default is \code{"recur_default"}.
#' @param multicore If to use parallel computing. Default is \code{TRUE}.
#' @param plot If to plot a bargraph to visualize vi and the ranking. Default is \code{TRUE}
#' @param Title Figure title. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param xLabel X-axis label. Make sure to use quotation marks. Use \code{NULL} to hide. Default is \code{NULL}.
#' @param yLabel Y-axis label. Make sure to use quotatio marks. Use \code{NULL} to hide. Default is \code{"Mean Decrease in Accurac"}
#' @param n Number of features to show. Takes integer numbers. Default is \code{"all"} (make sure to include quotation marks).
#' @param errorbar The type of errorbar in the graph. Options are \code{"SEM"} (standard error of the mean) or \code{"SD"} (standard deviation). Default is \code{"SEM"}.
#' @param errorbarWidth The width of the errorbar. Default is \code{0.2}.
#' @param symbolSize Size of the symbol. Default is \code{2}.
#' @param xTxtSize Font size for the x-axis text. Default is \code{10}.
#' @param yTxtSize Font size for the y-axis text. Default is \code{10}.
#' @param plotWidth The width of the figure for the final output figure file. Default is \code{170}.
#' @param plotHeight The height of the figure for the final output figure file. Default is \code{150}.
#' @return Outputs a \code{list} object with  OOB error rate summary, and a joint-point curve in \code{csv} format.
#' @details Make sure to arrange data (dfm) with feature (e.g., gene) as variables (i.e., columns), and rownames as sample names.
#' @import ggplot2
#' @import foreach
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom randomForest randomForest importance
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @examples
#' \dontrun{
#' rbioRF_SFS(training_HCvTC, tgtVar_HCvTC, multicore = TRUE)
#' }
#' @export
rbioRF_SFS <- function(objTitle = "x_vs_tgt",
                       x, targetVar, nTimes = 50, nTree = 1001, mTry = "recur_default",
                       multicore = TRUE,
                       plot = TRUE, n = "all",
                       Title = NULL, xLabel = NULL, yLabel = NULL,
                       errorbar = "SEM", errorbarWidth = 0.2,
                       symbolSize = 2, xTxtSize = 10, yTxtSize =10,
                       plotWidth = 170, plotHeight = 150
){
  ## prepare the dataframe
  training <- data.frame(x, check.names = FALSE)

  ### pepare the target variable
  tgt <- factor(as.character(targetVar), levels = unique(targetVar))

  ### prepare draw size. this uses down-sampling if the samples are unbalanced
  nlvl <- length(levels(tgt))
  size <- min(as.vector(table(tgt))) # down-sampling
  drawSize <- rep(size, nlvl)

  ## prepare blank tree OOB error matrics
  singleerrmtx <- matrix(nrow = 1, ncol = nTimes) # for the recursive OOB error rates from a single tree
  ooberrmtx <- matrix(nrow = ncol(training), ncol = nTimes) # for the recursive OOB error rates from all trees.

  if (!multicore){

    ## signle core computing: recursive structure
    tmpFunc <- function(n, m, tmperrmtx, tmpTraining, tmpTgt,
                        tmpTree, tmpTry, tmpSize){

      if (n == 0){
        return(tmperrmtx)

      } else {
        if (tmpTry == "recur_default"){

          if (ncol(tmpTraining) < 4){
            rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, importance = TRUE,
                               proximity = TRUE, drawSize = tmpSize)
          } else {
            rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree, mtry = max(floor(ncol(tmpTraining) / 3), 2),
                               importance = TRUE,
                               proximity = TRUE, drawSize = tmpSize)
          }

        } else if (tmpTry == "rf_default"){
          rf <- randomForest(x = tmpTraining, y = tmpTgt, ntree = tmpTree,
                             importance = TRUE,
                             proximity = TRUE, drawSize = tmpSize)
        } else {
          stop("Please select a proper mtry setting")
        }

        tmperrmtx[, m] <- tail(rf$err.rate[, 1], n = 1) # fill the OOB error rate
        tmpFunc(n - 1, m + 1, tmperrmtx, tmpTraining, tmpTgt,
                tmpTree, tmpTry, tmpSize)
      }
    }

    tmpFunc2 <- function(i, j, tmp2mtx, ...){
      if (i == 0){
        rownames(tmp2mtx) <- seq(j - 1)
        colnames(tmp2mtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

        return(tmp2mtx)
      } else {
        tmp2mtx[j, ] <- tmpFunc(n = nTimes, m = 1, tmperrmtx = singleerrmtx,
                                tmpTraining = training[, 1:j, drop = FALSE], tmpTgt = tgt, tmpTree = nTree, tmpTry = mTry,
                                tmpSize = drawSize)
        tmpFunc2(i - 1, j + 1, tmp2mtx, ...)
      }
    }

    mtxforfunc2 <- ooberrmtx

    ooberrmtx <- tmpFunc2(i = ncol(training), j = 1, tmp2mtx = mtxforfunc2) # j is the tree index

  } else { ## parallel computing
    # set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # recursive RF using par-apply functions
    tmpfunc3 <- function(j){
      if (mTry == "recur_default"){
        if (j < 4){
          rf <- randomForest::randomForest(x = training[, 1:j, drop = FALSE], y = tgt, ntree = nTree, importance = TRUE,
                                           proximity = TRUE, drawSize = drawSize)

        } else {
          rf <- randomForest::randomForest(x = training[, 1:j, drop = FALSE], y = tgt, ntree = nTree, mtry = max(floor(ncol(training[1:j]) / 3), 2),
                                           importance = TRUE,
                                           proximity = TRUE, drawSize = drawSize)
        }
      } else if (mTry == "rf_default"){
        rf <- randomForest::randomForest(x = training[, 1:j, drop = FALSE], y = tgt, ntree = nTree, importance = TRUE,
                                         proximity = TRUE, drawSize = drawSize)
      } else {
        stop("Please select a proper mtry setting")
      }

      tmperrmtx <- tail(rf$err.rate[, 1], n = 1) # compute the OOB error rate
      lst <- list(tmperrmtx = tmperrmtx)
    }

    l <- foreach(i = 1:ncol(training), .packages = c("foreach")) %dopar% {
      tmp <- foreach(j = 1:nTimes) %dopar% tmpfunc3(i)
      errmtx <- foreach(i = 1:nTimes, .combine = cbind) %dopar% tmp[[i]]$tmperrmtx
      lst <- list(errmtx = errmtx)
    }
    ooberrmtx <- foreach(j = 1:ncol(training), .combine = rbind) %dopar% l[[j]]$errmtx

    rownames(ooberrmtx) <- seq(ncol(training))
    colnames(ooberrmtx) <- c(paste("OOB_error_tree_rep", seq(nTimes), sep = "_"))

  }

  ## perpare the summary dataframe for OOB error rates
  ooberrnames <- rownames(ooberrmtx)
  ooberrmean <- rowMeans(ooberrmtx)
  ooberrSD <- apply(ooberrmtx, 1, sd)
  ooberrSEM <- sapply(ooberrSD, function(x)x / sqrt(ncol(ooberrmtx)))
  ooberrsummary <- data.frame(Features = ooberrnames, Mean = ooberrmean, SD = ooberrSD,
                              SEM = ooberrSEM, stringsAsFactors = FALSE)
  ooberrsummary$Features <- factor(ooberrsummary$Features, levels = unique(ooberrsummary$Features))

  ## output
  mean_min_idx <- which.min(ooberrsummary$Mean)  # index for the minimum mean oob feature group
  sd_min <- ooberrsummary$SD[mean_min_idx]  # oob SD for the feature group above
  minerrsd <- with(ooberrsummary, which(Mean <= (Mean[mean_min_idx] + sd_min)))  # 1sd minimum selection

  minfeatures <- colnames(training)[1:min(minerrsd)]
  sfsmatrix <- training[, 1:min(minerrsd), drop = FALSE]

  outlst <- list(selected_features = minfeatures,
                 feature_subsets_with_min_OOBerror_plus_1SD = minerrsd,
                 OOB_error_rate_summary = ooberrsummary,
                 SFS_matrix = sfsmatrix)

  sink(file = paste(objTitle,".SFS.txt",sep = ""), append = FALSE) # dump the results to a file
  print(outlst)
  sink() # end dump

  ## plot
  if (plot){
    # check the feature number
    if (nrow(ooberrsummary) == 1){
      ## print msg
      print("Only single feature subset detected. No need to plot.")
      ## output to env
      return(assign(paste(objTitle, "_SFS", sep = ""), outlst, envir = .GlobalEnv))

    } else {
      loclEnv <- environment()
      # prepare plotting dataframe
      if (n != "all"){
        pltdfm <- head(ooberrsummary, n)
      } else {
        pltdfm <- ooberrsummary
      }

      # plotting
      baseplt <- ggplot(pltdfm, aes(x = Features, y = Mean, group = 1), environment = loclEnv) +
        geom_line() +
        geom_point(size = symbolSize) +
        scale_x_discrete(expand = c(0.01, 0)) +
        ggtitle(Title) +
        xlab(xLabel) + # the arguments for x and y labls are switched as the figure is rotated
        ylab(yLabel) + # the arguments for x and y labls are switched as the figure is rotated
        geom_vline(xintercept = min(minerrsd), linetype = "dashed") +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom",
              legend.title = element_blank(),
              axis.text.x = element_text(size = xTxtSize),
              axis.text.y = element_text(size = yTxtSize, hjust = 0.5))

      if (errorbar == "SEM"){
        plt <- baseplt +
          geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = errorbarWidth, position = position_dodge(0.9)) +
          scale_y_continuous(expand = c(0, 0),
                             limits = c(with(pltdfm, min(Mean - SEM) * 0.6),
                                        with(pltdfm, max(Mean + SEM) * 1.2)))
      } else if (errorbar == "SD") {
        plt <- baseplt +
          geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbarWidth, position = position_dodge(0.9)) +
          scale_y_continuous(expand = c(0, 0),
                             limits = c(with(pltdfm, min(Mean - SD) * 0.6),
                                        with(pltdfm, max(Mean + SD) * 1.2)))
      }

      grid.newpage()

      # extract gtable
      pltgtb <- ggplot_gtable(ggplot_build(plt))

      # add the right side y axis
      Aa <- which(pltgtb$layout$name == "axis-l")
      pltgtb_a <- pltgtb$grobs[[Aa]]
      axs <- pltgtb_a$children[[2]]
      axs$widths <- rev(axs$widths)
      axs$grobs <- rev(axs$grobs)
      axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
      Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
      pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
      pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

      # export the file and draw a preview
      ggsave(filename = paste(objTitle,".OOB.plot.pdf", sep = ""), plot = pltgtb,
             width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(x)) converts object name into a character string
      grid.draw(pltgtb) # preview
    }
  }

  ## output to env
  return(assign(paste(objTitle, "_SFS", sep = ""), outlst, envir = .GlobalEnv))
}
